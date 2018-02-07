# Here are some of the commands used to calculate some of the simulations shown in the paper. GNU Parallel (https://www.gnu.org/software/parallel/) greatly simplifies invoking the script repeatedly with varying parameters, while handling CPU load. To test parallel commands without actually running them, use --dryrun as the first flag.

# calculate oxidation rates when the oxidase concentration is varied from 500 to 3000 enzymes per cell (25 to 125 20-mers) (`seq 25 5 150`) at two diffusion coefficients (0.45 450), Each calculation is invoked in pentuplicate (`seq 5`), running 8 calculations concurrently (-j 8).
parallel -j 8 python sim11dec-c.py --lengthscale 2.5 --name SingleDH-20 --number_of_goals {3} 126 --kcats 1604 300 --diffusion_coefficient {2} --particles 10000 --goal_radius 0.017888544 --goal_slowdown 13.3 --num_timesteps 5000 --timestep_length 0.0005 >> vary-TO-outputfile ::: `seq 5` ::: 0.45 450 ::: `seq 25 5 150`

# refers to the conditions.tsv file to calculate oxidation rates over a range of diffusion coefficients (`seq 0 0.02 0.6`). Each calculation is invoked in pentuplicate (`seq 5`), running 16 calculations concurrently (-j 16).
parallel --colsep '\t' -j 16 python ETCsimulation.py --lengthscale 2.5 --name {1} --number_of_goals {2} {3} --kcats {4} {5} --diffusion_coefficient {9} --particles 10000 --goal_radius {6} --goal_slowdown {7} --num_timesteps 5000 --timestep_length 0.0005 >> outputfile :::: conditions.tsv ::: `seq 5` ::: `seq 0 0.02 0.6`

# plot the frames for a movie with tethered goals, and also output the cumulative quinol density curve.
python ETCsimulation.py --calculate_average_distribution --lengthscale 2.5 --number_of_goals 37 --spawn_tethered --diffusion_coefficient 0.1 --kcats 1611 1676 --goal_radius 0.017888543819998 --particles 50000 --goal_slowdown 7.4 --timestep_length 0.0005 --num_timesteps 1000 --name NADH_domains20_tethered --verbose --plot_all_frames imageoutputfile --tethering 0.02 0.5 > distributionoutputfile
