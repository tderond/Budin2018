from __future__ import print_function
from __future__ import division
import numpy
import time
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import argparse




parser = argparse.ArgumentParser(description="""
    Note that the DIFFUSION_COEFFICIENT is in real time units, while the other two movement arguments are in terms of timesteps
    """)
parser.add_argument("--particles",type=int,#nargs=1,
                    required=True,help="total size of the quinone pool")

timesteps=parser.add_mutually_exclusive_group(required=True)
timesteps.add_argument("--num_timesteps",type=int,#nargs=1,
              help="use explicit number of timesteps.")
timesteps.add_argument("--total_time",type=int,#nargs=1,
              help="use a total simulation time. in unit of seconds. number of timesteps end up being total_time/TIMESTEP_LENGTH")
timesteps.add_argument("--over_variance",type=float,#nargs=1,
              help="number of timesteps becomes this number divided by the movement step variance (so it scales inversely with the diffusion coefficient)")


speed_or_diffconst=parser.add_mutually_exclusive_group(required=True)
speed_or_diffconst.add_argument("--diffusion_coefficient",type=float,#nargs=1,
                                help="like a self-diffusion coefficient, in units of um^2/sec. The standard deviation of the gaussian taken every step will be sqrt(2*DIFFUSION_COEFFICIENT*TIMESTEP_LENGTH)")
speed_or_diffconst.add_argument("--speed",type=float,#nargs=1,
                                help="average displacement in um every timestep. To achieve this, the standard deviation of the gaussian taken every step will be SPEED*sqrt(pi/2). In units of um/TIMESTEP_LENGTH (so if the timestep length is 1ms, the unit is um/ms or mm/sec)")
speed_or_diffconst.add_argument("--movement_standard_deviation",type=float,#nargs=1,
                                help="directly set the standard deviation of movement during each timestep. In units of um/TIMESTEP_LENGTH (so if the timestep length is 1ms, the unit is um/ms or mm/sec)")
parser.add_argument("--timestep_length",type=float,nargs="?",default=0.001,
                    help="how much time each simulation step represents, in seconds. Default is 0.001s or 1ms")

parser.add_argument("--lengthscale",type=float,nargs="?",default=1,
                    help="the simulation plane is a square. This number represents how much length in um each of the sides represents. Default = 1um")



parser.add_argument("--kcats",type=float,nargs=2,required=True,
                    help="first argument is kcat of terminal oxidase, second argument is kcat of reducase. can't be more than 1/timestep_length")



parser.add_argument("--goal_slowdown",type=float,nargs="?",default=4.0,
                    help="Make the goals this much slower than the particles (i.e. divide their movement standard deviation by this number, (default=4.0)")

parser.add_argument("--goal_radius",type=float,nargs="+",
                    required=True,help="radius of the enzymes. If one argument, number of each type of goal, if two arguments, set number of terminal oxidases and reductases separately")
                    
parser.add_argument("--number_of_goals",type=int,nargs="+",
                    required=True,help="if one argument, number of each type of goal, if two arguments, set number of terminal oxidases and reductases separately")

parser.add_argument("--tethering",type=float,nargs=2,
                    help="Tethers akin to a spring. First argument is the spring's ground state length and the second argument is the spring constant. If the number of oxidizing and reducing goals is not the same, it will only tether until it runs out of one of them.")

parser.add_argument("--spawn_tethered",action='store_true',help="sets the initial locations of the reducing goals to be next to their tethered oxidizing goal. only possible if the number of oxidizing goals equal the number of reducing goals")

parser.add_argument("--start_measuring_stats",type=float,nargs="?",default=0.5,
                    help="Start measuring statistics this far through the runs (default=0.5, which means the statistics are based on only the second half of timesteps)")

parser.add_argument("--calculate_average_distribution",action='store_true')

parser.add_argument("--verbose",action='store_true')

parser.add_argument("--plot_history",action='store_true')

parser.add_argument("--realtime_output",action='store_true',
                    help="output stats in a machine-readable format to stdout every step")

parser.add_argument("--plot_final_frame",type=str,#nargs=1,
             help="plots the final frame of the simulation as png and svg, argument is filename, without extension")
parser.add_argument("--plot_all_frames",type=str,#nargs=1,
             help="plots every frame of the simulation as png, with an svg every 100 frames, argument is filename prefix, will get appended with frame number and extension")


parser.add_argument("--emphasis",type=int,nargs="?",default=0,const=2)

parser.add_argument("--name",type=str,nargs="?",default="no_name_given",
                    help="Just a name to help you keep track of the output")


#arguments="--verbose --lengthscale 2.5 --name tethered --number_of_goals 800 900 --kcats 1611 990 --diffusion_coefficient 0.1 --particles 10000 --goal_radius 0.004 --goal_slowdown 2.3 --num_timesteps 5000 --timestep_length 0.0005 --tethering 0.01 0.1"
#namespace=parser.parse_args(arguments.split())
namespace=parser.parse_args()




time0=time.time()


timestep_length=namespace.timestep_length
scale=namespace.lengthscale
number_of_particles=namespace.particles



if(namespace.diffusion_coefficient):
    movement_stdev=numpy.sqrt(namespace.diffusion_coefficient*2.0*timestep_length)/scale
elif(namespace.speed):
    movement_stdev=namespace.speed*numpy.sqrt(0.5*numpy.pi)/scale
elif(namespace.movement_standard_deviation):
    movement_stdev=namespace.movement_standard_deviation/scale
else:
    sys.exit(1)
    
if(namespace.total_time):
    number_of_timesteps=namespace.total_time/timestep_length
elif(namespace.over_variance):
    number_of_timesteps=int(namespace.over_variance/(movement_stdev*movement_stdev))
elif(namespace.num_timesteps):
    number_of_timesteps=namespace.num_timesteps
else:
    sys.exit(1)    
    
number_of_oxidizing_goals=namespace.number_of_goals[0]
if(len(namespace.number_of_goals)==1):
    number_of_reducing_goals=number_of_oxidizing_goals
elif(len(namespace.number_of_goals)==2):
    number_of_reducing_goals=namespace.number_of_goals[1]
else:
    sys.exit(1)
    
    
oxidizing_goal_radius=namespace.goal_radius[0]/scale
if(len(namespace.goal_radius)==1):
    reducing_goal_radius=oxidizing_goal_radius
elif(len(namespace.goal_radius)==2):
    reducing_goal_radius=namespace.goal_radius[1]/scale
else:
    sys.exit(1)

    
max_kcat=1/timestep_length
if((namespace.kcats[0]>max_kcat)or(namespace.kcats[1]>max_kcat)):
    print("kcats cannot be larger than {:1.1f}".format(max_kcat),file=sys.stderr)
    sys.exit(1)
oxidation_probability=namespace.kcats[0]*timestep_length
reduction_probability=namespace.kcats[1]*timestep_length
    


if(namespace.plot_all_frames):
    fig_prefix=namespace.plot_all_frames
else:
    fig_prefix=""
    
if(namespace.plot_final_frame):
    final_fig_prefix=namespace.plot_all_frames
else:
    final_fig_prefix=""

emphasis_time=namespace.emphasis

if(namespace.verbose):
    verbose=True
else:
    verbose=False
    
if(namespace.realtime_output):
    realtime_output=True
else:
    realtime_output=False

if(namespace.tethering):
    tethering=True
    tethering_distance=namespace.tethering[0]/scale
    tethering_strength=namespace.tethering[1]
    number_to_tether=min(number_of_oxidizing_goals,number_of_reducing_goals)
    if(number_of_oxidizing_goals!=number_of_reducing_goals):
        if(verbose):
            print("@{:1.1f}s: Number of oxidizing and reducing goals not the same, only tethering the first {:d}".format(
                time.time()-time0,number_to_tether),file=sys.stderr)

else:
    tethering=False


goalmovement_stdev=movement_stdev/namespace.goal_slowdown #make goals move slower


#plotting parameters
fig_every_so_many=1
oxidizing_goal_color="LightBlue"
oxidizing_goal_plot_radius=max(oxidizing_goal_radius,0)
oxidizing_goal_alpha=0.7
reducing_goal_color="Pink"
reducing_goal_plot_radius=max(reducing_goal_radius,0)
reducing_goal_alpha=0.7
particle_size=5 # (not a radius. radius scales by square root)

# clustering_analysis_radius=oxidizing_goal_radius*5
# output_clustering_every=10
output_oxrate_every=1

densityradii_stepsize=0.005
density_radii=numpy.arange(0.01,0.26,densityradii_stepsize) #used for calculating the densities near goals


expected_wellmixed_rate=0.005*number_of_particles #used for plotting


#simulation setup

particles=numpy.random.rand(number_of_particles,2) #initualize the quinones

expected_reduced=number_of_particles/2
if (verbose):
    print("setting half the particles to initially reduced, which is {}".format(expected_reduced),file=sys.stderr)

particle_states=numpy.zeros(number_of_particles,dtype='int')
particle_states[0:int(expected_reduced)]=1 #set half the quinones to reduced to begin with
particle_state_ages=numpy.zeros(number_of_particles,dtype='int')

oxidizing_goals=numpy.random.rand(number_of_oxidizing_goals,2)

if(namespace.spawn_tethered):
    reducing_goals=oxidizing_goals+numpy.random.normal(scale=tethering_distance*1.25,size=(number_of_oxidizing_goals,2)) #this also forces the number of reducing goals to be the same as the oxidizing goals
else:
    reducing_goals=numpy.random.rand(number_of_reducing_goals,2)

    
#stats kept track of during the simulation
total_winners=0
clustering_total=0
near_oxidized_goal_total=0
total_reduced=0
oxrate_tally=0

#for plotting various graphs
oxrate_history=[]
conc_reduced_history=[]
conc_near_history=[]

#for plotting the movement history of a particle/quinone
particle0_xhistory=[]
particle0_yhistory=[]
particle0_redoxhistory=[]

cumulative_totals_within_radii=numpy.zeros(len(density_radii))

calculated_dic={"movement_stdev":movement_stdev,
     "number_of_timesteps":number_of_timesteps
     }


def tesselate(array):
    #there's probably a smarter way to do this but it works.
    x1=numpy.array([1,0],dtype="float")
    y1=numpy.array([0,1],dtype="float")
    r=numpy.append(
        numpy.append(
            numpy.append(
                numpy.append(
                    numpy.append(
                        numpy.append(
                            numpy.append(
                                numpy.append(
                                    array,
                                    array-y1,axis=0),
                                array+y1,axis=0),
                            array-x1,axis=0),
                        array+x1,axis=0),
                    array+x1+y1,axis=0),
                array+x1-y1,axis=0),
            array-x1+y1,axis=0),
        array-x1-y1,axis=0)
    return r



def do_tethering(oxidizing_goals,reducing_goals):
    if tethering==False:
        return (oxidizing_goals,reducing_goals)
    if(verbose):
        print("@{:1.1f}s: Tethering goals for timestep {:d}".format(time.time()-time0,t),file=sys.stderr)
    new_oxidizing_goals=oxidizing_goals.copy()
    new_reducing_goals=reducing_goals.copy()
    for i in xrange(number_to_tether):
        og=oxidizing_goals[i]
        rg=reducing_goals[i]
        tesselated_rg=tesselate(rg[numpy.newaxis,:])
        distance=numpy.linalg.norm(og-tesselated_rg,axis=1)
        closest_index=numpy.argmin(distance)
        #because of how the field wraps, an OG's buddy RG may be a tesselated version of the true RG. pick the closest one. 
        movement=(og-tesselated_rg[closest_index])/distance[closest_index]*(tethering_distance-distance[closest_index])*tethering_strength
        #you can think of the "left" side of this as calculating the unit vector encoding the direction of this movement vector
        #while the "right" side (after the first multiplication) calculates the magnitude of this movement vector.
        new_oxidizing_goals[i]+=movement*.5
        new_reducing_goals[i]-=movement*.5
        #the tethering movement gets doled out equally between the two buddy goals
    return (new_oxidizing_goals,new_reducing_goals)
            
            

def push_goals_apart(oxidizing_goals,reducing_goals)
    if(verbose):
        print("@{:1.1f}s: Preventing goal overlap for timestep {:d}".format(time.time()-time0,t),file=sys.stderr)
    #all_goals=oxidizing_goals.append(reducing_goals,axis=1)
    tesselated_oxidizing_goals=tesselate(oxidizing_goals)
    tesselated_reducing_goals=tesselate(reducing_goals)
    new_oxidizing_goals=oxidizing_goals.copy()
    new_reducing_goals=reducing_goals.copy()
    for og in tesselated_oxidizing_goals:
        #"og" is pushing
        oxidizing_comparison_array=oxidizing_goals-og
        oxidizing_distances=numpy.linalg.norm(oxidizing_comparison_array,axis=1)
        test=numpy.logical_and(oxidizing_distances<(oxidizing_goal_radius*2), #for those that are overlapping
                               numpy.all(oxidizing_goals!=og,axis=1)) #make sure you don't push yourself away.
        new_oxidizing_goals[test]+=(oxidizing_comparison_array[test]/oxidizing_distances[test,numpy.newaxis]*(2*oxidizing_goal_radius-oxidizing_distances[test,numpy.newaxis]))\
                                *0.5        
        reducing_comparison_array=reducing_goals-og
        reducing_distances=numpy.linalg.norm(reducing_comparison_array,axis=1)
        test=reducing_distances<(oxidizing_goal_radius+reducing_goal_plot_radius) #no need to check for self here.
        new_reducing_goals[test]+=reducing_comparison_array[test]/reducing_distances[test,numpy.newaxis]*((oxidizing_goal_radius+reducing_goal_radius)-reducing_distances[test,numpy.newaxis])\
                                *(oxidizing_goal_radius/(oxidizing_goal_radius+reducing_goal_radius))
            
            
    for rg in tesselated_reducing_goals:
        #"rg" is pushing
        oxidizing_comparison_array=oxidizing_goals-rg
        oxidizing_distances=numpy.linalg.norm(oxidizing_comparison_array,axis=1)
        test=oxidizing_distances<(oxidizing_goal_radius+reducing_goal_radius) #no need to check for self.
        new_oxidizing_goals[test]+=(oxidizing_comparison_array[test]/oxidizing_distances[test,numpy.newaxis]*(reducing_goal_radius+oxidizing_goal_radius-oxidizing_distances[test,numpy.newaxis]))\
                                *(reducing_goal_radius/(oxidizing_goal_radius+reducing_goal_radius))
            
        reducing_comparison_array=reducing_goals-rg
        reducing_distances=numpy.linalg.norm(reducing_comparison_array,axis=1)
        test=numpy.logical_and(reducing_distances<(reducing_goal_radius*2), #for those that are close
                               numpy.all(reducing_goals!=rg,axis=1)) #make sure you don't push yourself away.
        new_reducing_goals[test]+=reducing_comparison_array[test]/reducing_distances[test,numpy.newaxis]*((2*reducing_goal_radius)-reducing_distances[test,numpy.newaxis])\
                                *0.5
        
    return (new_oxidizing_goals,new_reducing_goals)

def plot_image(filename,final=False):
        fig=matplotlib.pyplot.figure(figsize=(12,9))
        gs=matplotlib.gridspec.GridSpec(3,4)
        ax =fig.add_subplot(gs[:,1:])
        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(r"${:0.1f} \mu m \times \, {:0.1f} \mu m$".format(scale,scale))

        #plot goals
        for og in oxidizing_goals:
            ax.add_patch(matplotlib.patches.Circle((og[0],og[1]),radius=oxidizing_goal_plot_radius
                                                   ,linewidth=0,alpha=oxidizing_goal_alpha,color=oxidizing_goal_color))
            ax.add_patch(matplotlib.patches.Circle((og[0]+1,og[1]),radius=oxidizing_goal_plot_radius
                                                   ,linewidth=0,alpha=oxidizing_goal_alpha,color=oxidizing_goal_color))
            ax.add_patch(matplotlib.patches.Circle((og[0],og[1]+1),radius=oxidizing_goal_plot_radius
                                                   ,linewidth=0,alpha=oxidizing_goal_alpha,color=oxidizing_goal_color))
            ax.add_patch(matplotlib.patches.Circle((og[0]-1,og[1]),radius=oxidizing_goal_plot_radius
                                                   ,linewidth=0,alpha=oxidizing_goal_alpha,color=oxidizing_goal_color))
            ax.add_patch(matplotlib.patches.Circle((og[0],og[1]-1),radius=oxidizing_goal_plot_radius
                                                   ,linewidth=0,alpha=oxidizing_goal_alpha,color=oxidizing_goal_color))

        for rg in reducing_goals:
            ax.add_patch(matplotlib.patches.Circle((rg[0],rg[1]),radius=reducing_goal_plot_radius
                                                   ,linewidth=0,alpha=reducing_goal_alpha,color=reducing_goal_color))
            ax.add_patch(matplotlib.patches.Circle((rg[0]+1,rg[1]),radius=reducing_goal_plot_radius
                                                   ,linewidth=0,alpha=reducing_goal_alpha,color=reducing_goal_color))
            ax.add_patch(matplotlib.patches.Circle((rg[0],rg[1]+1),radius=reducing_goal_plot_radius
                                                   ,linewidth=0,alpha=reducing_goal_alpha,color=reducing_goal_color))
            ax.add_patch(matplotlib.patches.Circle((rg[0]-1,rg[1]),radius=reducing_goal_plot_radius
                                                   ,linewidth=0,alpha=reducing_goal_alpha,color=reducing_goal_color))
            ax.add_patch(matplotlib.patches.Circle((rg[0],rg[1]-1),radius=reducing_goal_plot_radius
                                                   ,linewidth=0,alpha=reducing_goal_alpha,color=reducing_goal_color))


                
        
        
        #plot tethers
        if tethering:
            for i in xrange(number_to_tether):
                og=oxidizing_goals[i]
                rg=reducing_goals[i]
                if((og[0]-rg[0])*(og[0]-rg[0])+(og[1]-rg[1])*(og[1]-rg[1])<0.1):
                    ax.add_line(matplotlib.lines.Line2D([og[0],rg[0]],[og[1],rg[1]],color="gold",linewidth=2))
                
        if(final):
            fig.savefig("{}_noparticles.png".format(filename))
            fig.savefig("{}_noparticles.svg".format(filename))
        
        
        #plot movement history path for particle 0
        if(namespace.plot_history):    
            number_of_lines=len(particle0_xhistory)-1
            for i in range(number_of_lines):
                heavyness=(i/number_of_lines)*0.7+0.3
                ax.add_line(matplotlib.lines.Line2D([particle0_xhistory[i],particle0_xhistory[i+1]],
                                        [particle0_yhistory[i],particle0_yhistory[i+1]],
                                        linewidth=1,
                                        color=[particle0_redoxhistory[i]*0.7*heavyness+(1-heavyness),
                                               (1-heavyness),
                                               (1-particle0_redoxhistory[i])*0.7*heavyness+(1-heavyness)]))

        
        
        #plot particles/quinones
        colors=numpy.full(number_of_particles,"Blue",dtype="str")
        colors[particle_states==1]="Red"
        newly_changed=particle_state_ages<0
        emphatic_scatter=ax.scatter(particles[newly_changed,0],particles[newly_changed,1],
                                    s=10*pow((2.0*emphasis_time-particle_state_ages[newly_changed])/emphasis_time,2),
                                    c="Yellow",linewidths=0) #supposed to emphasize newly-changed particles
        
        scatter=ax.scatter(particles[:,0],particles[:,1],s=particle_size,linewidths=0,c=colors) #actual plotting of quinones
        
        
        
        #plot the stats
        total_reduced_ax=fig.add_subplot(gs[1,0])
        total_reduced_ax.set_xlim(1,number_of_timesteps)

        total_reduced_ax.set_ylim(0,1)
        total_reduced_ax.set_yticklabels(['{:2.0f}%'.format(x*100) for x in total_reduced_ax.get_yticks()])
        total_reduced_ax.set_ylabel("Proportion of\n$Q$ pool reduced")
        total_reduced_ax.set_xlabel("Timestep",labelpad=0)
        total_reduced_ax.plot(range(1,len(conc_reduced_history)+1),conc_reduced_history)
        


        reduced_within_radii_ax=fig.add_subplot(gs[2,0])
        reduced_within_radii_ax.set_xlim(min(density_radii),max(density_radii))
        reduced_within_radii_ax.set_ylim(0,number_of_particles/(scale*scale))
        reduced_within_radii_ax.set_ylabel("Average $QH_2$ density\n(particles per $\mu m^2$)")
        reduced_within_radii_ax.set_xlabel("Distance from center of $TR$")
        reduced_within_radii_ax.bar(density_radii,totals_within_radii/(numpy.pi*density_radii*density_radii*number_of_oxidizing_goals),linewidth=0,width=densityradii_stepsize)
        reduced_within_radii_ax.set_xticklabels(['{:2.2f}'.format(x) for x in reduced_within_radii_ax.get_xticks()],size="xx-small")
        xax=reduced_within_radii_ax.get_xaxis()
        
        
        ox_rate_ax=fig.add_subplot(gs[0,0])
        ox_rate_ax.set_xlim(1,number_of_timesteps)
        ox_rate_ax.set_xlabel("Timestep",labelpad=0)
        ox_rate_ax.set_ylim(0,expected_wellmixed_rate*number_of_timesteps)
        ox_rate_ax.set_ylabel("Cumulative\noxidation events")
        ox_rate_ax.plot(range(output_oxrate_every-1+1,len(oxrate_history)*output_oxrate_every+1,output_oxrate_every),
                        numpy.cumsum(oxrate_history))
        
        
        
        ax.set_title("""Timestep #{:05d}: {:03d} oxidation events. {:d} out of {:d} quinones reduced.""".format(
            t+1,total_winners,numpy.count_nonzero(particle_states),number_of_particles)
            ,loc="left",size=12)
        if(final):
            fig.savefig("{}.png".format(filename))
            fig.savefig("{}.svg".format(filename))
        else:
            fig.savefig("{}_t{:06d}.png".format(filename,t))
            if(t%100==0):
                fig.savefig("{}_t{:06d}.svg".format(filename,t))
        matplotlib.pyplot.close()
        
        
#plotting setup

# colors=numpy.full(number_of_particles,"Blue",dtype="str")
# colors[particle_states==1]="Red"


winners_this_round=0
cumulative_steps=0
cumulative_winners=0
cumulative_reduced=0



if(realtime_output):
    print("REALTIME_OUTPUT",
          "runtime:","{:0.1f}".format(time.time()-time0),
          "name:",namespace.name,
          "timestep:",0,
          "elapsed_time (sec):",0,
          "oxidation_events_this_timestep:","N/A",
          "cumulative_oxidation_events:",0,
          "since_measuring_stats:",0,
          "num_reduced_quinones:",numpy.count_nonzero(particle_states),
          "program_name:",sys.argv[0],
          "parameters:",namespace,"calculated:",calculated_dic,
          sep="\t")



#START OF THE MAIN SIMULATION LOOP
for t in xrange(number_of_timesteps):

    
    #Apply Wiener process to quinones and enzymes
    particles+=numpy.random.normal(scale=movement_stdev,size=(number_of_particles,2))
    particles=1.0-numpy.modf(1.0-numpy.modf(particles)[0])[0] #wraps less than 0 and more than 1
    
    oxidizing_goals+=numpy.random.normal(scale=goalmovement_stdev,size=(number_of_oxidizing_goals,2))
    oxidizing_goals=1.0-numpy.modf(1.0-numpy.modf(oxidizing_goals)[0])[0] #wraps less than 0 and more than 1
    
    reducing_goals+=numpy.random.normal(scale=goalmovement_stdev,size=(number_of_reducing_goals,2))
    reducing_goals=1.0-numpy.modf(1.0-numpy.modf(reducing_goals)[0])[0] #wraps less than 0 and more than 1
    
    #invoke functions defined above to apply tethering and avoid goal overlap
    (oxidizing_goals,reducing_goals)=do_tethering(oxidizing_goals,reducing_goals)
    
    (oxidizing_goals,reducing_goals)=push_goals_apart(oxidizing_goals,reducing_goals)# prevent bumping into each other
    
    
    winners_this_round=0
    particle_state_ages+=1
    
    totals_within_radii=numpy.zeros(len(density_radii))
    
    
    new_particle_states=particle_states.copy()
    
    if(verbose):
        print("@{:1.1f}s: Checking particle distances to goals for timestep {:d}".format(time.time()-time0,t),file=sys.stderr)
    
    #this is used to implement the kcats
    luckyness=numpy.random.rand(number_of_particles)
    lucky_enough_to_get_oxidized=luckyness<oxidation_probability
    lucky_enough_to_get_reduced=luckyness<reduction_probability
    
    #test for particles in goals and density around goals
    for og in tesselate(oxidizing_goals):
        comparison_array=particles-og
        distances=numpy.linalg.norm(comparison_array,axis=1) #that's the same as
                #numpy.sqrt(numpy.power(comparison_array[:,0],2)+numpy.power(comparison_array[:,1],2))
                #I checked it gives the same results, but numpy.linalg.norm is about twice as fast
                #Using numpy.power(comparison_array[:,0],2)+numpy.power(comparison_array[:,1],2) or
                #numpy.sum(numpy.square(comparison_array),axis=1) and comparing to the squared radii
                #would also work, but is the same speed as numpy.linalg.norm
                
        #set the ones inside the goal to oxidized
        test=numpy.logical_and(
            numpy.logical_and(
                distances<=oxidizing_goal_radius,
                particle_states==1),
                lucky_enough_to_get_oxidized
            )
        number_of_winners=numpy.count_nonzero(test)
        new_particle_states[test]=0
        particle_state_ages[test]=-emphasis_time
        
        
        #calculate how many reduced are within various radii from the goals (using the "old" particle states before this step)
        if(namespace.calculate_average_distribution):
            for i in xrange(len(density_radii)):
                totals_within_radii[i]+=numpy.count_nonzero(numpy.logical_and(distances<=density_radii[i]/scale,new_particle_states==1))
                #count_nonzero is WAY faster than sum() here
        
        total_winners+=number_of_winners
        winners_this_round+=number_of_winners
                        
        
        
    for rg in tesselate(reducing_goals):
        comparison_array=particles-rg
        distances=numpy.linalg.norm(comparison_array,axis=1) #that's the same as
                #numpy.sqrt(numpy.power(comparison_array[:,0],2)+numpy.power(comparison_array[:,1],2))
                #I checked it gives the same results, but numpy.linalg.norm is about twice as fast
        test=numpy.logical_and(
            numpy.logical_and(
                distances<=reducing_goal_radius,
                particle_states==0),
                lucky_enough_to_get_reduced
            )
        new_particle_states[test]=1
        particle_state_ages[test]=-emphasis_time
    
    particle_states=new_particle_states
    total_reduced+=numpy.count_nonzero(particle_states)
    
#     if (t%100==0):
#         print("timestep #{:d}: {:d} out of {:d} are now reduced. {:d} oxidation events have happened".format(t,sum(particle_states),len(particle_states),total_winners),file=sys.stderr)
#         print(t,sum(particle_states),total_winners)

#     if(t==(number_of_timesteps/50-1)):
#         print("{:0.1f}".format(time.time()-time0),number_of_particles,number_of_timesteps,oxidizing_goal_radius,
#           speed,goalspeed,number_of_oxidizing_goals,number_of_reducing_goals,t,total_winners,sum(particle_states))
#         total_winners=0


#only start keeping track of these stats after a certain fraction of the simulation steps have been completed
    if(t>=number_of_timesteps*namespace.start_measuring_stats):
        cumulative_steps+=1
        cumulative_totals_within_radii+=totals_within_radii
        cumulative_winners+=winners_this_round
        cumulative_reduced+=numpy.count_nonzero(particle_states)
        
    if(namespace.plot_history):    
        particle0_xhistory.append(particles[0,0])
        particle0_yhistory.append(particles[0,1])
        particle0_redoxhistory.append(particle_states[0])
        
    if(realtime_output):
        print("REALTIME_OUTPUT",
              "runtime:","{:0.1f}".format(time.time()-time0),
              "name:",namespace.name,
              "timestep:",t+1,
              "elapsed_time (sec):",(t+1)*timestep_length,
              "oxidation_events_this_timestep:",winners_this_round,
              "cumulative_oxidation_events:",total_winners,
              "since_measuring_stats:",cumulative_winners,
              "num_reduced_quinones:",numpy.count_nonzero(particle_states),
              "program_name:",sys.argv[0],
              "parameters:",namespace,"calculated:",calculated_dic,
              sep="\t")
        
        
        
    conc_reduced_history.append(numpy.count_nonzero(particle_states)/len(particle_states))
    if(t%output_oxrate_every==output_oxrate_every-1):
        oxrate_history.append(oxrate_tally)
        oxrate_tally=0
    oxrate_tally+=winners_this_round
    
        
    if ((fig_prefix!="") and (t%fig_every_so_many==0)):
        
        if(verbose):
            print("@{:1.1f}s: Plotting timestep {:d}".format(time.time()-time0,t),file=sys.stderr)
        
        plot_image(fig_prefix) #function defined above
        
        
#END OF MAIN SIMULATION LOOP        
        
if(namespace.plot_final_frame):
        plot_image(namespace.plot_final_frame,final=True)

if(verbose):
    print("In {:d} timesteps, {:d} reached a goal (took {:0.1f} sec to calculate)".format(
        number_of_timesteps,total_winners,time.time()-time0),
         file=sys.stderr)

print("FINAL_STATS",
      "runtime:","{:0.1f}".format(time.time()-time0),
      "name:",namespace.name,
      "diffusioncoefficient (um^2/sec):",namespace.diffusion_coefficient,
      "oxidation_rate (/sec):",cumulative_winners/(cumulative_steps*timestep_length),
      "reduced_qs_proportion (reduced/total):",cumulative_reduced/(cumulative_steps*number_of_particles),
      "stats_over_number_of_steps:",cumulative_steps,
      "number_of_oxidizing_goals:", number_of_oxidizing_goals,
      "program_name:",sys.argv[0],
      "parameters:",namespace,"calculated:",calculated_dic,
      sep="\t")
if(namespace.calculate_average_distribution):
    densities=totals_within_radii/(numpy.pi*density_radii*density_radii*number_of_oxidizing_goals)
    average_densities=cumulative_totals_within_radii/(cumulative_steps*numpy.pi*density_radii*density_radii*number_of_oxidizing_goals)
    for i in xrange(len(density_radii)):
        print("QUINONE_DISTRIBUTION",
            "name:",namespace.name,
            "diffusion_coefficient (cm^2/s):",namespace.diffusion_coefficient,
            "distance_from_goal (um):",density_radii[i],
            "particle_density_last_step (particles/um^2):",densities[i],
            "particle_density_average (particles/um^2):",average_densities[i],
             sep="\t")
