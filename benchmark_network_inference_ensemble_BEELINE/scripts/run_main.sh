# this is to plot out the time usage for each of the GRN algorithms
Rscript plot_time_usage.R

# preprocess the networks 
# convert the threshold networks into concerete networks via the F1 maximization
Rscript preprocess_networks.R

# calculate for precision, recall and F1 for various networks 
Rscript precision_recall.R

# run the boolean simulation to find steady states 
Rscript Boolean_steady_states.R

# run simulation for iQcell 
# jupyter nbconvert --execute iQcell_simulation_steady_states.ipynb

# get the gold standard steady states 
Rscript find_golden_steady_states.R

# check the steady states of the boolean figures 
Rscript check_steady_states_boolean_sim.R

# plot out the OneSC manual curated data 
# jupyter nbconvert --execute save_cluster_G.ipynb
# jupyter nbconvert --execute plot_BEELINE_requirements.ipynb


