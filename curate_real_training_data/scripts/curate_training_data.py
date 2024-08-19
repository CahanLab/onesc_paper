import pandas as pd
import numpy as np
import scanpy as sc
import onesc
from os import path, mkdir, makedirs
import pickle
import seaborn
import networkx as nx

mmTFs = pd.read_csv("../input/mouse_TFs/Mus_musculus_TF.txt", sep = '\t')
mmTFs_list = list(mmTFs['Symbol'])

# this is to get the myeloid stuff 2 
test_adata = sc.read_h5ad("../input/myeloid_progenitors/filtered_adata.h5ad")
samp_tab = test_adata.obs

sc.tl.paga(test_adata, groups="cell_types")
sc.pl.paga(test_adata, color=["cell_types"])

exp_tab = test_adata.raw.to_adata().to_df()
exp_tab = exp_tab.T
pca_coord = test_adata.obsm['X_pca']
pca_coord = pd.DataFrame(pca_coord, index = samp_tab.index).T

trajectory_dict = dict()
trajectory_dict['T1'] = ['CMP', 'MEP', 'Erythrocytes']
trajectory_dict['T2'] = ['CMP', 'GMP', 'Granulocytes']
trajectory_dict['T3'] = ['CMP', 'GMP', 'Monocytes']
trajectory_dict['T4'] = ['CMP', 'MK']

cluster_col = 'cell_types'
pt_col = 'dpt_pseudotime'

initial_clusters = ['CMP']
end_clusters = ['Erythrocytes', 'Granulocytes', 'Monocytes', 'MK']
clusters_G = onesc.construct_cluster_network(exp_tab, samp_tab, 
                                             initial_clusters = initial_clusters, 
                                             terminal_clusters = end_clusters, 
                                             cluster_col = cluster_col, 
                                             pseudo_col = pt_col, 
                                             pt_margin = 0.12)
nx.draw(clusters_G, with_labels = True)

my_df = onesc.suggest_dynamic_genes(exp_tab.loc[exp_tab.index.isin(mmTFs_list), :].copy(), samp_tab, trajectory_dict, cluster_col, pt_col, adj_p_cutoff = 0.05, log2_change_cutoff = 3, min_exp_cutoff = 0.4)
interesting_TFs = np.unique(my_df.index)

curated_exp_tab = exp_tab.loc[interesting_TFs, :]
if path.isdir('../output/myeloid_progenitors_full') == False:
    makedirs('../output/myeloid_progenitors_full')
curated_exp_tab.to_csv('../output/myeloid_progenitors_full/train_exp.csv')
samp_tab.to_csv('../output/myeloid_progenitors_full/samp_tab.csv')
pt_tab = samp_tab.loc[:, ['dpt_pseudotime']]
pt_tab.columns = ['PseudoTime']
pt_tab.to_csv('../output/myeloid_progenitors_full/pseudotime.csv')