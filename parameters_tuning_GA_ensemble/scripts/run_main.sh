
# start investigating the effect of clustering info on performance 
python screen_cluster_info.py 

# start investgating the effect of generation info on the performance 
python screen_generation_time.py 

python extract_generation_time_runtime.py

# investigate the effect of generation on the performance 
python screen_generation_time_finer.py 


# start investigating the effect of network density paramter on the performance 
python screen_network_density.py 

##### the below are to make plots #####
Rscript plot_cluster_info.R

Rscript plot_generation_time.R

Rscript plot_network_density.R