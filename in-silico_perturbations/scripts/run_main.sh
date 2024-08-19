# move all the required files for simulation of the perturbation 
python move_networks.py 

# run perturbation simulation -- takes a long time 
python perturb_simulations.py 

# count the proportion of end states 
python count_proportion.py 

# plot out the proportions of steady states after perturbation 
Rscript plot_proportion.R

# curate the data needed for UMAP -- takes a long time 
python curate_UMAP_data.py 

# compute the correlational similarity 
python compute_correlation_similarity.py 

# plot out the UMAP simulations 
Rscript plot_UMAP_simulation.R

# plot out the heatmap simulations 
# Rscript plot_heatmap_simulation.R
