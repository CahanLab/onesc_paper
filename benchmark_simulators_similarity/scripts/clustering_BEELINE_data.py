import pandas as pd 
import numpy as np 
import OneCC_bool
import os
import scanpy as sc
import string

if "BEELINE_clusters" not in os.listdir("../output/"):
    os.mkdir("../output/BEELINE_clusters")

# perform clustering of the BEELINE data

def cluster_adata_synthetic(data_type, resolution = 0.3):

    file_path = "../BEELINE-data/inputs/Synthetic/" + data_type + "/" + data_type + "-5000-1/"
    if data_type == 'dyn-LI_Extended':
        file_path = '../BoolODE/Debug/dyn-linear-1_extended/dyn-linear-1_extended-5000-1/'

    train_exp = pd.read_csv(file_path + "ExpressionData.csv", index_col=0)
    pt_st = pd.read_csv(file_path + "PseudoTime.csv", index_col = 0)

    pt_st = pt_st.fillna(0)
    pt_st['pseudoTime'] = pt_st.sum(axis = 1)
    pt_st = pt_st.loc[train_exp.columns, :]

    adata = sc.AnnData(train_exp.T)
    sc.tl.pca(adata, svd_solver='arpack')

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=9)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)
    sc.pl.umap(adata, color='leiden')
    adata.obs['pseudoTime'] = pt_st['pseudoTime']

    return adata


def assign_alphabet(adata):
    cluster_list = np.unique(adata.obs['leiden'])
    pt_list = list()
    for unique_cluster in cluster_list:
        pt_list.append(np.mean(adata.obs.loc[adata.obs['leiden'] == unique_cluster, 'pseudoTime']))
    assign_df = pd.DataFrame()
    assign_df['leiden'] = cluster_list
    assign_df['mean_pseudoTime'] = pt_list
    assign_df = assign_df.sort_values("mean_pseudoTime")
    alphabet_string = string.ascii_uppercase
    alphabet_list = list(alphabet_string)
    alphabet_list = ['cluster_' + x for x in alphabet_list]
    assign_df['cluster_id'] = alphabet_list[0:len(cluster_list)]
    assign_df.index = assign_df['leiden']

    adata.obs['cluster_id'] = None
    for unique_cluster in assign_df.index:
        adata.obs.loc[adata.obs['leiden'] == str(unique_cluster), 'cluster_id'] = assign_df.loc[unique_cluster, 'cluster_id']
    return adata
# dyn-BF 
data_type = 'dyn-BF'
adata = cluster_adata_synthetic(data_type, resolution = 0.3)
adata = assign_alphabet(adata)
sc.pl.umap(adata, color = 'cluster_id')

if data_type not in os.listdir("../output/BEELINE_clusters/"):
    os.mkdir("../output/BEELINE_clusters/" + data_type)
target_dir = "../output/BEELINE_clusters/" + data_type + "/"
adata.write_h5ad(target_dir + "redefined_adata.h5ad")
adata.obs.to_csv(target_dir + "compiled_sampTab_pt.csv")
adata.to_df().to_csv(target_dir + "compiled_expTab.csv")

# dyn-BFC 
data_type = 'dyn-BFC'
adata = cluster_adata_synthetic(data_type, resolution = 0.3)
adata = assign_alphabet(adata)
sc.pl.umap(adata, color = 'cluster_id')

if data_type not in os.listdir("../output/BEELINE_clusters/"):
    os.mkdir("../output/BEELINE_clusters/" + data_type)
target_dir = "../output/BEELINE_clusters/" + data_type + "/"
adata.write_h5ad(target_dir + "redefined_adata.h5ad")
adata.obs.to_csv(target_dir + "compiled_sampTab_pt.csv")
adata.to_df().to_csv(target_dir + "compiled_expTab.csv")

# dyn-LI
data_type = 'dyn-LI'
adata = cluster_adata_synthetic(data_type, resolution = 0.3)
adata = assign_alphabet(adata)
sc.pl.umap(adata, color = 'cluster_id')

if data_type not in os.listdir("../output/BEELINE_clusters/"):
    os.mkdir("../output/BEELINE_clusters/" + data_type)
target_dir = "../output/BEELINE_clusters/" + data_type + "/"
adata.write_h5ad(target_dir + "redefined_adata.h5ad")
adata.obs.to_csv(target_dir + "compiled_sampTab_pt.csv")
adata.to_df().to_csv(target_dir + "compiled_expTab.csv")

# dyn-TF
data_type = 'dyn-TF'
adata = cluster_adata_synthetic(data_type, resolution = 0.3)
adata.obs.loc[adata.obs['leiden'] == '8', 'leiden'] = '6'
adata = assign_alphabet(adata)
sc.pl.umap(adata, color = 'cluster_id')

if data_type not in os.listdir("../output/BEELINE_clusters/"):
    os.mkdir("../output/BEELINE_clusters/" + data_type)
target_dir = "../output/BEELINE_clusters/" + data_type + "/"
adata.write_h5ad(target_dir + "redefined_adata.h5ad")
adata.obs.to_csv(target_dir + "compiled_sampTab_pt.csv")
adata.to_df().to_csv(target_dir + "compiled_expTab.csv")

# dyn-CY
data_type = 'dyn-CY'
adata = cluster_adata_synthetic(data_type, resolution = 0.1)
adata = assign_alphabet(adata)
sc.pl.umap(adata, color = 'cluster_id')
sc.pl.umap(adata, color = 'pseudoTime')
# switch cluster A and cluster B
old_A_index = adata.obs['cluster_id'] == 'cluster_A'
old_B_index = adata.obs['cluster_id'] == 'cluster_B'
adata.obs.loc[old_A_index, 'cluster_id'] = 'cluster_B'
adata.obs.loc[old_B_index, 'cluster_id'] = 'cluster_A'
sc.pl.umap(adata, color = 'cluster_id')

if data_type not in os.listdir("../output/BEELINE_clusters/"):
    os.mkdir("../output/BEELINE_clusters/" + data_type)
target_dir = "../output/BEELINE_clusters/" + data_type + "/"
adata.write_h5ad(target_dir + "redefined_adata.h5ad")
adata.obs.to_csv(target_dir + "compiled_sampTab_pt.csv")
adata.to_df().to_csv(target_dir + "compiled_expTab.csv")

# dyn-LI_Extended
data_type = 'dyn-LI_Extended'
adata = cluster_adata_synthetic(data_type, resolution = 0.2)
adata = assign_alphabet(adata)
sc.pl.umap(adata, color = 'cluster_id')
sc.pl.umap(adata, color = 'pseudoTime')

if data_type not in os.listdir("../output/BEELINE_clusters/"):
    os.mkdir("../output/BEELINE_clusters/" + data_type)
target_dir = "../output/BEELINE_clusters/" + data_type + "/"
adata.write_h5ad(target_dir + "redefined_adata.h5ad")
adata.obs.to_csv(target_dir + "compiled_sampTab_pt.csv")
adata.to_df().to_csv(target_dir + "compiled_expTab.csv")

