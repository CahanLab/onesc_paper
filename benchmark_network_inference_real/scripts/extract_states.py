import pandas as pd
import numpy as np 
import onesc
import os

data_types = ['myeloid_progenitors_full']

output_path = '../output/extract_states/'
for data_type in data_types:
    if data_type == 'day_4_EB':
        initial_clusters = ['Primed']
        end_clusters = ['PGC', 'NeurEct', 'PrimStr']
    elif data_type == 't_cell_diff':
        initial_clusters = ['ETP']
        end_clusters = ['DN3']
    elif data_type == 'myeloid_progenitors_full':
        initial_clusters = ['CMP']
        end_clusters = ['Erythrocytes', 'Granulocytes', 'Monocytes', 'MK']
    train_exp = pd.read_csv("../Beeline_benchmark/inputs/real_data/" + data_type + "/train_exp.csv", index_col = 0)
    samp_tab = pd.read_csv("../Beeline_benchmark/inputs/real_data/" + data_type + "/samp_tab.csv", index_col = 0)
    pt_col = 'dpt_pseudotime'
    cluster_col = 'cell_types'

    clusters_G = onesc.construct_cluster_network(train_exp, samp_tab, initial_clusters = initial_clusters, terminal_clusters = end_clusters, cluster_col = cluster_col, pseudo_col = pt_col)

    lineage_cluster = onesc.extract_trajectory(clusters_G,initial_clusters, end_clusters)

    vector_thresh = onesc.find_threshold_vector(train_exp, samp_tab, cluster_col = "cell_types", cutoff_percentage=0.4)
    state_dict = onesc.define_states(train_exp, samp_tab, lineage_cluster, vector_thresh, cluster_col, percent_exp = 0.3)

    if os.path.isdir(os.path.join(output_path, data_type)) == False:
        os.makedirs(os.path.join(output_path, data_type))
    for lineage_name in state_dict.keys():
        temp_state = state_dict[lineage_name]
        temp_state.to_csv(os.path.join(output_path, data_type) + "/" + lineage_name + "_states.csv")
