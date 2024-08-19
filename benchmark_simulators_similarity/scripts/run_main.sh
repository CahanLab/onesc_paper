# run boolODE to extend out the Linear dynamics 
cd ../BoolODE/
python boolode.py --config config-files/dyn-LI_extended.yaml

cd ../scripts
# run the simulation of the networks 
python simulate_networks.py 

# cluster the BEELINE data 
python clustering_BEELINE_data.py 

# compile sim profiles 
python compile_sim_profiles.py 

# check UMAP structure 
python check_UMAP_overlap.py

# check cross correlation 
python check_cross_correlation.py 

# plot out the cross correlation plots 
Rscript plot_cross_correlation.R

# perform pearson and spearman correlation to assess similarities 
python run_correlation_sim.py

# average out the correlation 
python average_correlation_sim.py

# plot out average correlation 
Rscript plot_average_correlation.R 

# plot out the original UMAP + cluster id 
# plot out the query UMAP + cluster id 
Rscript plot_correlation.R

##### these scripts are not used ######

# check UMAP structure using query 
python check_UMAP_overlap_query_UMAP.py

# TODO turns out that pySCN is pretty useless in this 
# perform pySCN to assess the similarities 
python run_pySCN.py

# plot out original UMAP + cluster id SCN classification 
# plot out query UMAP + cluster id
Rscript plot_classification.R