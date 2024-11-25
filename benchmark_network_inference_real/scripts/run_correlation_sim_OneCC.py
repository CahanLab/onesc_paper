import pandas as pd
import numpy as np 
import os
import scanpy as sc 
import pickle

if "run_correlation" not in os.listdir("../output/"):
    os.mkdir("../output/run_correlation")

data_type_list = ['t_cell_diff', 'HSC', 'day_4_EB']

def create_avg_profile(data_type):
    train_exp = pd.read_csv("../Beeline_benchmark/inputs/real_data/" + data_type + "/train_exp.csv", index_col=0)
    train_st = pd.read_csv("../Beeline_benchmark/inputs/real_data/" + data_type + "/samp_tab.csv", index_col=0)

    avg_df = pd.DataFrame(data = None, index = train_exp.index, columns = np.unique(train_st['cell_types']))
    for cluster_id in avg_df.columns:
        temp_exp = train_exp.loc[:, train_st['cell_types'] == cluster_id]
        avg_df.loc[:, cluster_id] = temp_exp.mean(axis = 1)
    avg_df = avg_df.T
    avg_df = avg_df - avg_df.min()
    avg_df = avg_df / avg_df.max()
    avg_df = avg_df * 2
    return avg_df.T

def correlation_similarity(avg_df, query_exp):
    corr_mat = pd.DataFrame(data = None, index = avg_df.columns, columns = query_exp.columns)
    for cluster_id in corr_mat.index:
        print(cluster_id)
        corr_mat.loc[cluster_id, :] = query_exp.corrwith(avg_df[cluster_id])
    return corr_mat 

for data_type in data_type_list: 
    if data_type not in os.listdir("../output/run_correlation"):
        os.mkdir("../output/run_correlation/" + data_type)
        
    target_dir = "../output/run_correlation/" + data_type + '/'

    avg_df = create_avg_profile(data_type)
    
    if data_type == 'dyn-LI_Extended':
        data_type = 'dyn-LI'
    query_exp = pd.read_csv("../output/simulation_OneCC/" + data_type + "/run_OneCC_pyEpoch/density/big_sim_df.csv", index_col = 0)
    query_st = pd.read_csv("../output/simulation_OneCC/" + data_type + "/run_OneCC_pyEpoch/density/UMAP_coord.csv", index_col = 0)

    avg_df = avg_df.loc[query_exp.index, :]
    corr_mat = correlation_similarity(avg_df, query_exp)

    corr_mat.to_csv(target_dir + "correlation_mat.csv")
    avg_df.to_csv(target_dir + "avg_df.csv")

for data_type in data_type_list:
    print(data_type)
    target_dir = "../output/run_correlation/" + data_type + '/'
    corr_mat = pd.read_csv(target_dir + "correlation_mat.csv", index_col=0)

    if data_type == 'dyn-LI_Extended':
        data_type = 'dyn-LI'
    query_st = pd.read_csv("../output/simulation_OneCC/" + data_type + "/run_OneCC_pyEpoch/density/UMAP_coord.csv", index_col = 0)
    query_st = query_st.loc[corr_mat.columns, :]
    query_st['cluster_id'] = corr_mat.idxmax()
    query_st.to_csv(target_dir + 'sample_tab.csv')
