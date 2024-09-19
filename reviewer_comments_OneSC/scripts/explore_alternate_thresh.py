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
from scipy import stats

warnings.filterwarnings("ignore")

data_type = 'myeloid_progenitors_full'

initial_clusters = ['CMP']
end_clusters = ['Erythrocytes', 'Granulocytes', 'Monocytes', 'MK']

train_exp = pd.read_csv("../input/curated_myeloid_data" + "/train_exp.csv", index_col = 0)
samp_tab = pd.read_csv("../input/curated_myeloid_data" + "/samp_tab.csv", index_col = 0)

pt_col = 'dpt_pseudotime'
cluster_col = 'cell_types'

def percentile_threshold(train_exp, percentile_cut = 0.25): 
    cutoff_dict = dict()
    for tmp_gene in train_exp.index: 
        cutoff_dict[tmp_gene] = np.quantile(train_exp.loc[tmp_gene, :][train_exp.loc[tmp_gene, :] != 0], percentile_cut)
    return pd.Series(cutoff_dict)

quant_thresh = percentile_threshold(train_exp)

#########################
output_path = "../output/alternate_thresh/"
os.makedirs(output_path, exist_ok=True)

edge_list = [("MEP", "MK"), ("CMP", "MEP"), ("MEP", "Erythrocytes"), ("CMP", "GMP"), ("GMP", "Granulocytes"), ("GMP", "Monocytes")]
clusters_G = nx.DiGraph(edge_list)

f = plt.figure()
nx.draw(clusters_G, with_labels = True)

lineage_cluster = onesc.extract_trajectory(clusters_G,initial_clusters, end_clusters)

vector_thresh = onesc.find_threshold_vector(train_exp, samp_tab, cluster_col = "cell_types", cutoff_percentage=0.4)

####### plot out the correlation between the two methods 
slope, intercept, r_value, p_value, std_err = stats.linregress(vector_thresh, quant_thresh)

# Regression line values
regression_line = slope * vector_thresh + intercept

# Create scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(vector_thresh, quant_thresh, color='blue', label='Data points')
plt.plot(vector_thresh, regression_line, color='red', label=f'Linear fit: y = {slope:.2f}x + {intercept:.2f}')

# Add labels and title
plt.xlabel('Current thresholding scheme')
plt.ylabel('Percentile thresholding scheme')
plt.title('Percentile thresh vs current thresh')
plt.legend()
plt.savefig(os.path.join(output_path, 'quant_vs_curr_thresh.png'))

####### check the Boolean profiles ######
lineage_time_change_dict = onesc.find_gene_change_trajectory(train_exp, samp_tab, lineage_cluster, cluster_col, pt_col, vector_thresh, pseudoTime_bin=0.01) # this pseudotime bin could be 

state_dict = onesc.define_states(train_exp, samp_tab, lineage_cluster, vector_thresh, cluster_col, percent_exp = 0.3)
state_dict_new = onesc.define_states(train_exp, samp_tab, lineage_cluster, quant_thresh, cluster_col, percent_exp = 0.3)

for tmp_traj in state_dict.keys():
    print(state_dict[tmp_traj])
    print(state_dict_new[tmp_traj])
    print(state_dict[tmp_traj] == state_dict_new[tmp_traj])
