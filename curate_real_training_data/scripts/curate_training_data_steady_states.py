import pandas as pd
import numpy as np
import scanpy as sc
import onesc
from os import path, mkdir, makedirs
import os
import networkx as nx

data_type = 'myeloid_progenitors_full'
initial_clusters = ['CMP']
end_clusters = ['Erythrocytes', 'Granulocytes', 'Monocytes', 'MK']

train_exp = pd.read_csv("../output/" + data_type + "/train_exp.csv", index_col = 0)
samp_tab = pd.read_csv("../output/" + data_type + "/samp_tab.csv", index_col = 0)

pt_col = 'dpt_pseudotime'
cluster_col = 'cell_types'

samp_tab.index = [str(x) for x in samp_tab.index]
train_exp.columns = [str(x) for x in train_exp.columns]

clusters_G = onesc.construct_cluster_network(train_exp, samp_tab, initial_clusters = initial_clusters, terminal_clusters = end_clusters, cluster_col = cluster_col, pseudo_col = pt_col)

lineage_cluster = onesc.extract_trajectory(clusters_G,initial_clusters, end_clusters)
vector_thresh = onesc.find_threshold_vector(train_exp, samp_tab, cluster_col = "leiden",  cutoff_percentage=0.45)

state_dict = onesc.define_states(train_exp, samp_tab, lineage_cluster, vector_thresh, cluster_col)
transition_dict = onesc.define_transition(state_dict)

for temp_lineage in state_dict.keys():
    state_df = state_dict[temp_lineage]
    state_df.to_csv("../output/" + data_type + "/" + temp_lineage + "_states.csv")
