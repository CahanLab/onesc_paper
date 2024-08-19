from operator import index
import pandas as pd
import numpy as np 
import os
import scanpy as sc 
import pickle

data_types = os.listdir("../output/run_correlation/")
data_types = [x for x in data_types if 'dyn' in x]

for data_type in data_types: 
    print(data_type)
    query_cluster_corr = pd.read_csv("../output/run_correlation/" + data_type + "/correlation_mat.csv", index_col = 0)
    query_cluster_st = pd.read_csv("../output/run_correlation/" + data_type + "/sample_tab.csv", index_col = 0)
    trajectory_id = list(query_cluster_st['cell_id'])
    trajectory_id = [x.split("_")[0] for x in trajectory_id]
    query_cluster_st['trajectory_id'] = trajectory_id

    if data_type == 'dyn-CY':
        df_list = list()
        for trajectory in np.unique(query_cluster_st['trajectory_id']):
            small_traj_st = query_cluster_st.loc[query_cluster_st['trajectory_id'] == trajectory, :]
            small_corr = query_cluster_corr.loc[:, small_traj_st.index]
            small_corr.columns = [x.split("_")[1] for x in small_corr.columns]
            df_list.append(small_corr)
        mean_df = pd.concat(df_list).groupby(level=0).mean()
        min_df = pd.concat(df_list).groupby(level=0).min()
        max_df = pd.concat(df_list).groupby(level=0).max()

        mean_df.to_csv("../output/run_correlation/" + data_type + "/" + "mean_corr.csv")
        min_df.to_csv("../output/run_correlation/" + data_type + "/" + "min_corr.csv")
        max_df.to_csv("../output/run_correlation/" + data_type + "/" + "max_corr.csv")
    else:
        for ss in np.unique(query_cluster_st['steady_states']):
            sub_traj_st = query_cluster_st.loc[query_cluster_st['steady_states'] == ss, :]
            df_list = list()
            for trajectory in np.unique(sub_traj_st['trajectory_id']):
                small_traj_st = sub_traj_st.loc[sub_traj_st['trajectory_id'] == trajectory, :]
                small_corr = query_cluster_corr.loc[:, small_traj_st.index]
                small_corr.columns = [x.split("_")[1] for x in small_corr.columns]
                df_list.append(small_corr)
            mean_df = pd.concat(df_list).groupby(level=0).mean()
            min_df = pd.concat(df_list).groupby(level=0).min()
            max_df = pd.concat(df_list).groupby(level=0).max()

            mean_df.to_csv("../output/run_correlation/" + data_type + "/" + ss + "_mean_corr.csv")
            min_df.to_csv("../output/run_correlation/" + data_type + "/" + ss + "_min_corr.csv")
            max_df.to_csv("../output/run_correlation/" + data_type + "/" + ss + "_max_corr.csv")


