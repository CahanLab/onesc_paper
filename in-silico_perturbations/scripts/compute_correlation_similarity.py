import pandas as pd
import numpy as np 
import os
import scanpy as sc 
import pickle

def correlation_similarity(avg_df, query_exp):
    corr_mat = pd.DataFrame(data = None, index = avg_df.columns, columns = query_exp.columns)
    for cluster_id in corr_mat.index:
        print(cluster_id)
        corr_mat.loc[cluster_id, :] = query_exp.corrwith(avg_df[cluster_id])
    return corr_mat 

def create_avg_bool_profile():
    input_path = os.path.join("../output/files_needed/myeloid_progenitors_full/extract_states")
    combined_df = pd.DataFrame()
    for temp_file in os.listdir(input_path):
        temp_df = pd.read_csv(os.path.join(input_path, temp_file), index_col = 0)
        combined_df = pd.concat([combined_df, temp_df], axis = 1)
        print(temp_file)
    combined_df = combined_df.T.drop_duplicates().T
    return combined_df 

root_output = "../output/plot_UMAPs/"
perturb_types = os.listdir(root_output)
perturb_types = [x for x in perturb_types if 'DS_Store' not in x]

avg_df = create_avg_bool_profile()
for perturb_type in perturb_types: 
    genes_list = os.listdir(os.path.join(root_output, perturb_type))
    genes_list = [x for x in genes_list if "DS_Store" not in x]
    for temp_gene in genes_list:
        output_path = os.path.join(root_output, perturb_type, temp_gene)
        query_exp = pd.read_csv(os.path.join(output_path, "big_sim_df.csv"), index_col = 0)
        query_st = pd.read_csv(os.path.join(output_path, "UMAP_coord.csv"), index_col = 0)
        avg_df = avg_df.loc[query_exp.index, :]
        corr_mat = correlation_similarity(avg_df, query_exp)
        corr_mat.to_csv(os.path.join(output_path, "correlation_mat.csv"))
        avg_df.to_csv(os.path.join(output_path, "avg_df.csv"))



