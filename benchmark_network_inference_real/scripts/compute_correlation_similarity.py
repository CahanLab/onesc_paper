import pandas as pd
import numpy as np 
import os
import scanpy as sc 
import pickle

data_type_list = ['myeloid_progenitors_full']

def create_avg_profile(data_type):
    train_exp = pd.read_csv("../Beeline_benchmark/inputs/real_data/" + data_type + "/train_exp.csv", index_col=0)
    train_st = pd.read_csv("../Beeline_benchmark/inputs/real_data/" + data_type + "/samp_tab.csv", index_col=0)

    avg_df = pd.DataFrame(data = None, index = train_exp.index, columns = np.unique(train_st['cell_types']))
    for cluster_id in avg_df.columns:
        temp_exp = train_exp.loc[:, train_st['cell_types'] == cluster_id]
        avg_df.loc[:, cluster_id] = temp_exp.mean(axis = 1)
    avg_df = avg_df.T
    avg_df = avg_df - avg_df.min()
    avg_df = avg_df / (avg_df.max() - avg_df.min())
    avg_df = avg_df * 2
    return avg_df.T

def correlation_similarity(avg_df, query_exp):
    corr_mat = pd.DataFrame(data = None, index = avg_df.columns, columns = query_exp.columns)
    for cluster_id in corr_mat.index:
        print(cluster_id)
        corr_mat.loc[cluster_id, :] = query_exp.corrwith(avg_df[cluster_id], method = 'pearson')
    return corr_mat 

def create_avg_bool_profile(data_type):
    input_path = os.path.join("../output/extract_states", data_type)
    combined_df = pd.DataFrame()
    for temp_file in os.listdir(input_path):
        temp_df = pd.read_csv(os.path.join(input_path, temp_file), index_col = 0)
        combined_df = pd.concat([combined_df, temp_df], axis = 1)
        print(temp_file)
    combined_df = combined_df.T.drop_duplicates().T
    return combined_df 

for data_type in data_type_list:  
    target_dir = "../output/plot_UMAPs/" + data_type + '/'
    #avg_df = create_avg_profile(data_type)
    avg_df = create_avg_bool_profile(data_type)
    query_exp = pd.read_csv("../output/plot_UMAPs/" + data_type + "/big_sim_df.csv", index_col = 0)
    query_st = pd.read_csv("../output/plot_UMAPs/" + data_type + "/UMAP_coord.csv", index_col = 0)
    avg_df = avg_df.loc[query_exp.index, :]
    corr_mat = correlation_similarity(avg_df, query_exp)
    corr_mat.to_csv(target_dir + "correlation_mat.csv")
    avg_df.to_csv(target_dir + "avg_df.csv")


