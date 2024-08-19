Rscript preprocess_networks.R

python simulation_OneCC.py 

# extract the different states 
python extract_states.py

# compile steady states 
python compile_steady_states.py

# plot out the networks 
python plot_networks.py

# plot out steady states information
# Rscript plot_steady_states.R

# plot out all the steady states
Rscript plot_steady_states_profiles.R

# rerun this script for better plotting 
python curate_UMAP_data.py

# run correlation similarity between simulated data and real data 
python compute_correlation_similarity.py

# make pretty UMAPs of simulated profiles 
Rscript plot_UMAP_simulation.R

# measure cross correlation 
python measure_cross_correlation.py

# compile the cross correlation 
python compiled_cross_correlation.py

# compile dtw 
python compiled_dtw.py 

# plot out the cross-correlation 
Rscript plot_cross_correlation.R

# plot out the steady states similarity 
Rscript ss_similarity.R