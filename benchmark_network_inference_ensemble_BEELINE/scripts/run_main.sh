# this is to plot out the time usage for each of the algorithms
Rscript plot_time_usage.R

# preprocess the networks 
Rscript preprocess_networks.R

# check for precision and recall 
Rscript precision_recall.R

# run the AUPRC 
Rscript AUPRC.R

# run the boolean simulation to find steady states 
Rscript Boolean_steady_states.R

# run simulation for iQcell 
# jupyter nbconvert --execute iQcell_simulation_steady_states.ipynb

# check the steady states of the boolean figures 
Rscript check_steady_states_boolean_sim.R


