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

import argparse

import igraph as ig
from igraph import Graph

parser = argparse.ArgumentParser(description='machine type')
parser.add_argument('--machine', type=str, help = 'machine type')
args = parser.parse_args()

file_name = args.machine 

data_type = 'myeloid_progenitors_full'
initial_clusters = ['CMP']
end_clusters = ['Erythrocytes', 'Granulocytes', 'Monocytes', 'MK']

train_exp = pd.read_csv("../inputs/real_data/" + data_type + "/train_exp.csv", index_col = 0)
samp_tab = pd.read_csv("../inputs/real_data/" + data_type + "/samp_tab.csv", index_col = 0)

pt_col = 'dpt_pseudotime'
cluster_col = 'cell_types'

output_path = '../output/'

clusters_G = onesc.construct_cluster_network(train_exp, samp_tab, initial_clusters = initial_clusters, terminal_clusters = end_clusters, cluster_col = cluster_col, pseudo_col = pt_col)
onesc.plot_state_graph(clusters_G)

lineage_cluster = onesc.extract_trajectory(clusters_G,initial_clusters, end_clusters)

vector_thresh = onesc.find_threshold_vector(train_exp, samp_tab, cluster_col = "cell_types", cutoff_percentage=0.4)
lineage_time_change_dict = onesc.find_gene_change_trajectory(train_exp, samp_tab, lineage_cluster, cluster_col, pt_col, vector_thresh, pseudoTime_bin=0.01) # this pseudotime bin could be 

state_dict = onesc.define_states(train_exp, samp_tab, lineage_cluster, vector_thresh, cluster_col, percent_exp = 0.3)
transition_dict = onesc.define_transition(state_dict)
training_data = onesc.curate_training_data(state_dict, transition_dict, lineage_time_change_dict, samp_tab, cluster_id = cluster_col, pt_id = pt_col,act_tolerance = 0.04)
corr_mat = onesc.calc_corr(train_exp)

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
pickle.dump(time_lapse, open(os.path.join(output_path, file_name + '_time_lapse.pickle'), 'wb'))
networks_ensemble[0].to_csv(os.path.join(output_path, file_name + "_OneSC_network.csv"))