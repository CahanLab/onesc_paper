from operator import index
import pandas as pd
import numpy as np 
import os
import scanpy as sc 
import pickle

data_types = os.listdir("../output/run_correlation/")
data_types = [x for x in data_types if 'DS_Store' not in x]

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
            steady_states_dict[ss_profile] = [exp_file.removesuffix("_exp.csv")]
        else:
            #steady_states_dict[ss_profile].append(sim_exp)
            steady_states_dict[ss_profile].append(exp_file.removesuffix("_exp.csv"))
    return steady_states_dict

for data_type in data_types: 
    print(data_type)
    query_cluster_corr = pd.read_csv("../output/run_correlation/" + data_type + "/correlation_mat.csv", index_col = 0)
    query_cluster_st = pd.read_csv("../output/run_correlation/" + data_type + "/sample_tab.csv", index_col = 0)
    sim_folder_path = "../output/simulation_OneCC/" + data_type + "/run_OneCC_pyEpoch/density/simulation/"
    sim_cluster = cluster_simulation(sim_folder_path)
    inv_map = dict()
    for ss in sim_cluster.keys():
        for temp_name in sim_cluster[ss]:
            inv_map[temp_name] = ss
    
    query_cluster_st['trajectory_id'] = None
    query_cluster_st['steady_states'] = None
    for temp_index in query_cluster_st.index:
        query_cluster_st.loc[temp_index, 'trajectory_id'] = "_".join(temp_index.split("_", 2)[:2])
        query_cluster_st.loc[temp_index, 'steady_states'] = inv_map["_".join(temp_index.split("_", 2)[:2])]
    
    trajectory_id = list(query_cluster_st['trajectory_id'])
    trajectory_id = [x.split("_")[0] for x in trajectory_id]
    query_cluster_st['trajectory_id'] = trajectory_id


    for ss in np.unique(query_cluster_st['steady_states']):
        sub_traj_st = query_cluster_st.loc[query_cluster_st['steady_states'] == ss, :]
        df_list = list()
        for trajectory in np.unique(sub_traj_st['trajectory_id']):
            small_traj_st = sub_traj_st.loc[sub_traj_st['trajectory_id'] == trajectory, :]
            small_corr = query_cluster_corr.loc[:, small_traj_st.index]
            small_corr.columns = [x.split("_")[2] for x in small_corr.columns]
            df_list.append(small_corr)
        mean_df = pd.concat(df_list).groupby(level=0).mean()
        min_df = pd.concat(df_list).groupby(level=0).min()
        max_df = pd.concat(df_list).groupby(level=0).max()

        mean_df.to_csv("../output/run_correlation/" + data_type + "/" + ss + "_mean_corr.csv")
        min_df.to_csv("../output/run_correlation/" + data_type + "/" + ss + "_min_corr.csv")
        max_df.to_csv("../output/run_correlation/" + data_type + "/" + ss + "_max_corr.csv")


