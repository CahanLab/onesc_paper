import scanpy as sc
import pandas as pd
import numpy as np 
import onesc

import networkx as nx
import seaborn as sns
import pickle
import matplotlib.pyplot as plt

import time

import os
import math

import warnings
warnings.filterwarnings("ignore")

data_type = 'myeloid_progenitors_full'

initial_clusters = ['CMP']
end_clusters = ['Erythrocytes', 'Granulocytes', 'Monocytes', 'MK']

train_exp = pd.read_csv("../input/curated_myeloid_data" + "/train_exp.csv", index_col = 0)
samp_tab = pd.read_csv("../input/curated_myeloid_data" + "/samp_tab.csv", index_col = 0)

pt_col = 'dpt_pseudotime'
cluster_col = 'cell_types'

output_path = "../output/alternate_path/"
os.makedirs(output_path, exist_ok=True)

edge_list = [("MEP", "MK"), ("CMP", "MEP"), ("MEP", "Erythrocytes"), ("CMP", "GMP"), ("GMP", "Granulocytes"), ("GMP", "Monocytes")]
clusters_G = nx.DiGraph(edge_list)

f = plt.figure()
nx.draw(clusters_G, with_labels = True)
f.savefig(output_path + "/graph.png")

lineage_cluster = onesc.extract_trajectory(clusters_G,initial_clusters, end_clusters)

vector_thresh = onesc.find_threshold_vector(train_exp, samp_tab, cluster_col = "cell_types", cutoff_percentage=0.4)
lineage_time_change_dict = onesc.find_gene_change_trajectory(train_exp, samp_tab, lineage_cluster, cluster_col, pt_col, vector_thresh, pseudoTime_bin=0.01) # this pseudotime bin could be 

state_dict = onesc.define_states(train_exp, samp_tab, lineage_cluster, vector_thresh, cluster_col, percent_exp = 0.3)
transition_dict = onesc.define_transition(state_dict)
pickle.dump(state_dict, open(output_path + "/state_dict.pickle", "wb"))

training_data = onesc.curate_training_data(state_dict, transition_dict, lineage_time_change_dict, samp_tab, cluster_id = cluster_col, pt_id = pt_col,act_tolerance = 0.04, show_conflicts=True)

corr_mat = onesc.calc_corr(train_exp)

#weight_dict = onesc.define_weight(state_dict)

start_time = time.time() 
networks_ensemble = onesc.create_network_ensemble(training_data, 
                                    corr_mat, 
                                    ideal_edges = round(0.4 * corr_mat.shape[1]), 
                                    num_generations = 300, 
                                    max_iter = 30, 
                                    num_parents_mating = 4, 
                                    sol_per_pop = 30, 
                                    reduce_auto_reg = True, 
                                    GA_seed_list = [1, 2, 3, 4, 5], 
                                    init_pop_seed_list = [21, 22, 23, 24, 25])
time_lapse = time.time() - start_time

print(data_type)
print(time_lapse)
networks_ensemble[0].to_csv(output_path + "/OneSC_network.csv")