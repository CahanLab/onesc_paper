import pandas as pd
import numpy as np 
import os
import scanpy as sc 
import pickle
from scipy.spatial.distance import cdist

def correlation_similarity(avg_df, query_exp):
    corr_mat = pd.DataFrame(data = None, index = avg_df.columns, columns = query_exp.columns)
    for cluster_id in corr_mat.index:
        print(cluster_id)
        corr_mat.loc[cluster_id, :] = query_exp.corrwith(avg_df[cluster_id])
    return corr_mat 

def distance_similarity(avg_df, query_exp):
    dist_mat = pd.DataFrame(data = None, index = avg_df.columns, columns = query_exp.columns)
    for cluster_id in dist_mat.index: 
        dist_mat.loc[cluster_id, :] = cdist(query_exp.T, [np.array(avg_df[cluster_id]) * 2], metric = 'euclidean')[:, 0]
    return dist_mat

def distance_bool_similarity(avg_df, query_exp):
    dist_mat = pd.DataFrame(data = None, index = avg_df.columns, columns = query_exp.columns)
    bool_query_exp = query_exp.copy()
    bool_query_exp[bool_query_exp <= 1] = 0 
    bool_query_exp[bool_query_exp > 1] = 2
    for cluster_id in dist_mat.index: 
        bool_dist = np.abs(bool_query_exp.T - np.array(avg_df[cluster_id]) * 2).sum(axis = 1) / (avg_df.shape[0] * 2)
        dist_mat.loc[cluster_id, :] = bool_dist
    return dist_mat

def distance_bool_similarity_other(avg_df, query_exp):
    dist_mat = pd.DataFrame(data = None, index = avg_df.columns, columns = query_exp.columns)
    bool_query_exp = query_exp.copy()
    bool_query_exp[bool_query_exp <= 1] = 0 
    bool_query_exp[bool_query_exp > 1] = 2
    for cluster_id in dist_mat.index: 
        dist_mat.loc[cluster_id, :] = cdist(bool_query_exp.T, [np.array(avg_df[cluster_id]) * 2], metric = 'euclidean')[:, 0]
    return dist_mat

def create_avg_bool_profile(data_type):
    input_path = os.path.join("../input/extract_states", data_type)
    combined_df = pd.DataFrame()
    for temp_file in os.listdir(input_path):
        temp_df = pd.read_csv(os.path.join(input_path, temp_file), index_col = 0)
        combined_df = pd.concat([combined_df, temp_df], axis = 1)
        print(temp_file)
    combined_df = combined_df.T.drop_duplicates().T
    return combined_df 

data_type = 'myeloid_progenitors_full'
    
target_dir = "../output/test_zero_floor_OneSC/myeloid_progenitors_full/wildtype_UMAPs/"
#avg_df = create_avg_profile(data_type)
avg_df = create_avg_bool_profile(data_type)
query_exp = pd.read_csv("../output/test_zero_floor_OneSC/myeloid_progenitors_full/wildtype_UMAPs/big_sim_df.csv", index_col = 0)
query_st = pd.read_csv("../output/test_zero_floor_OneSC/myeloid_progenitors_full/wildtype_UMAPs/UMAP_coord.csv", index_col = 0)
avg_df = avg_df.loc[query_exp.index, :]
dist_mat = distance_similarity(avg_df, query_exp)
dist_mat.to_csv(target_dir + "dist_mat.csv")
bool_dist_mat = distance_bool_similarity(avg_df, query_exp)
bool_dist_mat.to_csv(target_dir + "bool_dist_mat.csv")


