import networkx as nx
import pandas as pd
import numpy as np
import random
import os
import math

def make_network(nx_object):
    network_df = pd.DataFrame(np.empty((len(nx_object.edges),3),dtype=pd.Timestamp),columns=['TF','TG','Type'])
    i = 0 
    for temp_edge in nx_object.edges:
        network_df.loc[i, 'TF'] = "g" + str(temp_edge[0])
        network_df.loc[i, 'TG'] = "g" + str(temp_edge[1])
        rand_type = random.randint(0, 1)
        if rand_type == 0:
            network_df.loc[i, 'Type'] = "+"
        else:
            network_df.loc[i, 'Type'] = "-"
        i = i + 1
    
    self_reg_df = pd.DataFrame(np.empty((len(nx_object.nodes),3),dtype=pd.Timestamp),columns=['TF','TG','Type'])
    all_genes = list(nx_object.nodes)
    all_genes = ["g"+str(x) for x in all_genes]
    self_reg_df['TF'] = all_genes
    self_reg_df['TG'] = all_genes
    self_reg_df['Type'] = "+"

    network_df = pd.concat([network_df, self_reg_df])
    return network_df

def output_string(reg_df, type = "+"):
    if reg_df.shape[0] == 0:
        return ""
    TF_list = list(reg_df['TF'])
    TF_list = [str(x) for x in TF_list]
    reg_string = " or ".join(TF_list)
    reg_string = "( " + reg_string + " )"
    if type == "+":
        return reg_string
    else:
        reg_string = "not" + reg_string
        return reg_string

def booleanize_network(network_df):
    booleanized_network_df = pd.DataFrame(np.empty((len(np.unique(network_df['TG'])),2),dtype=pd.Timestamp),columns=['Gene','Rule'])
    booleanized_network_df['Gene'] = np.unique(network_df['TG'])
    for target_gene in np.unique(network_df['TG']):
        sub_network_df = network_df.loc[network_df['TG'] == target_gene, :]

        pos_reg_df = sub_network_df.loc[sub_network_df['Type'] == "+", :]
        neg_reg_df = sub_network_df.loc[sub_network_df['Type'] == '-', :]

        pos_string = output_string(pos_reg_df, type = "+")
        neg_string = output_string(neg_reg_df, type = "-")

        if pos_string != "" and neg_string != "":
            combined_string = pos_string + " and " + neg_string
        elif pos_string == "" and neg_string != "":
            combined_string = neg_string
        elif pos_string != "" and neg_string == "":
            combined_string = pos_string
        booleanized_network_df.loc[booleanized_network_df['Gene'] == target_gene, 'Rule'] = combined_string
    return booleanized_network_df

# 0.2, 0.4, 0.6, 0.8, 1
# 20, 40, 60, 80, 100

if "random_networks" not in os.listdir("../output"):
    os.mkdir("../output/random_networks")

trial_runs = 10
proportion_list = [0.2, 0.4, 0.6, 0.8, 1]
num_node_list = [5, 10, 15, 20]

for trial in range(0, trial_runs):
    for proportion in proportion_list: 
        for num_node in num_node_list:
            num_edges = math.floor((num_node * num_node) * proportion - num_node)
            nx_object = nx.generators.random_graphs.gnm_random_graph(num_node, num_edges)
            network_df = make_network(nx_object)

            folder_name = "nodes_" + str(num_node) + "_edges_" + str(num_edges)
            my_path = "../output/random_networks/" + "trial_" + str(trial) + "/" + folder_name

            if not os.path.exists(my_path):
                os.makedirs(my_path)

            network_df.to_csv(my_path + "/OneCC_network.txt", sep = '\t', index = False)
            booleanized_network_df = booleanize_network(network_df)

            if not os.path.exists(my_path + "/BoolODE_original"):
                os.makedirs(my_path + "/BoolODE_original")
            
            if not os.path.exists(my_path + "/BoolODE_heaviside"):
                os.makedirs(my_path + "/BoolODE_heaviside")
            booleanized_network_df.to_csv(my_path + "/BoolODE_original/BoolODE_network.txt", sep = '\t', index = False)
            booleanized_network_df.to_csv(my_path + "/BoolODE_heaviside/BoolODE_network.txt", sep = '\t', index = False)


