import pandas as pd 
import numpy as np 
import onesc
import os
import seaborn as sns
import pickle 

# check for trajectory mapping 
# see if the trajectory of the clusters match the appropriate 

# Here are some metric suggesting 
# overlap OneCC and BoolODE together on the same UMAP 
# x axis == pseudotime; y axis == GRN similarities (euclidean distance) -- for individiual trajectory 
# x axis == pseudotime; y axis == pySCN scores -- individual trajectory 
    # cluster network plot. the length of the edge represents the pseudotime it requires to reach another cluster 
# cross correlation -- probably not needed 


# we perform the sampling of the data 

# meta table columns [cell id, simulation time, lineage]

if 'compiled_simulation' not in os.listdir('../output'):
    os.makedirs("../output/compiled_simulation")

def cluster_simulation(sim_folder_path): 
    all_files = os.listdir(sim_folder_path)
    all_files = [x for x in all_files if 'simulated_exp.csv' in x]

    steady_states_dict = dict()
    for exp_file in all_files: 
        sim_exp = pd.read_csv(sim_folder_path + exp_file, index_col=0)
        ss_profile = sim_exp.iloc[:, -1]
        ss_profile = np.array(ss_profile >= 1)
        ss_profile = ss_profile.astype(int)
        ss_profile = ss_profile.astype(str)
        ss_profile = "_".join(ss_profile)
        if ss_profile not in steady_states_dict.keys(): 
            #steady_states_dict[ss_profile] = [sim_exp]
            steady_states_dict[ss_profile] = [exp_file]
        else:
            #steady_states_dict[ss_profile].append(sim_exp)
            steady_states_dict[ss_profile].append(exp_file)
    return steady_states_dict

def sample_data(data_type):
    # cluster the simulation trajectory in terms of 
    steady_states_dict = cluster_simulation("../output/OneSC_sim/" + data_type + "/")
    samp_tab = pd.DataFrame()
    exp_tab = pd.DataFrame()
    for steady_states in steady_states_dict.keys():
        print(steady_states)
        for exp_file in steady_states_dict[steady_states]:
            cell_id = exp_file.split("_")[0]
            temp_exp = pd.read_csv("../output/OneSC_sim/" + data_type + "/" + exp_file, index_col = 0)
            temp_st = pd.DataFrame()
            temp_st['sim_time'] = temp_exp.columns
            temp_exp.columns = cell_id + "_" + temp_exp.columns
            temp_st['cell_id'] = temp_exp.columns
            temp_st.index = temp_exp.columns
            temp_st['steady_states'] = steady_states
            samp_tab = pd.concat([samp_tab, temp_st])
            exp_tab = pd.concat([exp_tab, temp_exp], axis = 1)
    exp_tab.to_csv('../output/compiled_simulation/' + data_type + "_exp.csv")
    samp_tab.to_csv('../output/compiled_simulation/' + data_type + "_st.csv")

# compile the simulated expression profiles 
for data_type in ['dyn-BF', 'dyn-BFC', 'dyn-TF', 'dyn-LI', 'dyn-CY']:
    sample_data(data_type)

