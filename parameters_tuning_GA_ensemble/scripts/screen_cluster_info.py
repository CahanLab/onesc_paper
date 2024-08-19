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
import pickle
import scanpy 

prng = np.random.default_rng(2023)

input_path = "../input/preprocessed_data"
subfolders = [ f.path for f in os.scandir(input_path) if f.is_dir() ]
subfolders = [x.strip("./") for x in subfolders]
subfolders = [x for x in subfolders if "ipynb" not in x]
subfolders = [x for x in subfolders if 'dyn-' in x]

def merge_clusters(samp_tab, clusters_G, cluster_col = 'leiden', num_merge = 2):
    i = 0
    while True:
        selected_cluster = prng.choice(np.unique(samp_tab[cluster_col]))
        print(selected_cluster)
        cur_out_edges = list(clusters_G.out_edges(selected_cluster))
        cur_in_edges = list(clusters_G.in_edges(selected_cluster))
        if len(cur_out_edges) == 0: 
            continue
        else:
            cur_edge_tuples = cur_out_edges[0] # if there are multiple out edges, only going to mix the first 1
            cur_node = cur_edge_tuples[0]
            next_node = cur_edge_tuples[1]
            next_out_edges = list(clusters_G.out_edges(next_node))
            clusters_G.add_node(cur_node + "+" + next_node)
            if len(cur_in_edges) > 0:
                for temp_tuple in cur_in_edges:
                    clusters_G.remove_edges_from([temp_tuple])
                    clusters_G.add_edges_from([(temp_tuple[0], cur_node + "+" + next_node)])
            if len(next_out_edges) > 0:
                for temp_tuple in next_out_edges:
                    clusters_G.remove_edges_from([temp_tuple])
                    clusters_G.add_edges_from([(cur_node + "+" + next_node, temp_tuple[1])])
            if len(cur_out_edges) > 1:
                for temp_tuple in cur_out_edges:
                    new_next_node = temp_tuple[1]
                    clusters_G.remove_edges_from([temp_tuple])
                    if new_next_node != next_node: 
                        clusters_G.add_edges_from([(cur_node + "+" + next_node, new_next_node)])
            clusters_G.remove_node(cur_node) # remove the edge 
            clusters_G.remove_node(next_node)
            i = i + 1
            cluster_info = np.array(samp_tab[cluster_col])
            cluster_info[np.array(samp_tab[cluster_col].isin([cur_node, next_node]))] = cur_node + "+" + next_node
            samp_tab[cluster_col] = cluster_info
        if i == num_merge:
            break 
    return [samp_tab, clusters_G] 

def split_clusters(samp_tab, clusters_G, cluster_col = 'leiden', pt_col = 'pseudoTime', num_merge = 2):
    i = 0
    while True: 
        selected_cluster = prng.choice(np.unique(samp_tab[cluster_col]))
        if "+part" in selected_cluster: 
            continue
        mid_pt = np.median(samp_tab.loc[samp_tab[cluster_col] == selected_cluster, pt_col])
        cur_out_edges = list(clusters_G.out_edges(selected_cluster))
        cur_in_edges = list(clusters_G.in_edges(selected_cluster))
        clusters_G.add_node(selected_cluster + "+part1")
        clusters_G.add_node(selected_cluster + "+part2")
        clusters_G.add_edges_from([(selected_cluster + "+part1", selected_cluster + "+part2")])
        if len(cur_in_edges) > 0: 
            for temp_tuple in cur_in_edges:
                clusters_G.remove_edges_from([temp_tuple])
                clusters_G.add_edges_from([(temp_tuple[0], selected_cluster + "+part1")])
        if len(cur_out_edges) > 0: 
            for temp_tuple in cur_out_edges:
                clusters_G.remove_edges_from([temp_tuple])
                clusters_G.add_edges_from([(selected_cluster + "+part2", temp_tuple[1])])
        clusters_G.remove_node(selected_cluster)
        i = i + 1
        cluster_info = np.array(samp_tab[cluster_col])
        cluster_info[list(np.logical_and(samp_tab[cluster_col] == selected_cluster, samp_tab[pt_col] > mid_pt))] = selected_cluster + "+part2"
        cluster_info[list(np.logical_and(samp_tab[cluster_col] == selected_cluster, samp_tab[pt_col] <= mid_pt))] = selected_cluster + "+part1"
        samp_tab[cluster_col] = cluster_info
        if i == num_merge: 
            break
    return [samp_tab, clusters_G]  

method_list = ['split', 'merge']
degree_list = [1, 2]
repeats_list = [1, 2, 3, 4, 5]
output_path = "../output/screen_cluster_assignment"
if os.path.exists(output_path) == False:
    os.makedirs(output_path)

for sub_folder in subfolders:
    for temp_method in method_list: 
        for temp_degree in degree_list: 
            for temp_trial in repeats_list:
                print(sub_folder)
                adata = sc.read_h5ad(os.path.join("../", sub_folder, "redefined_adata.h5ad"))
                exp_tab = adata.to_df().T
                samp_tab = adata.obs
                pt_col = 'pseudoTime'
                cluster_col = 'leiden'
                clusters_G = pickle.load(open(os.path.join("../", sub_folder, 'clusters_G.pickle'), "rb"))
                if temp_method == 'split':
                    [new_samp_tab, new_clusters_g] = split_clusters(samp_tab, 
                                                                    clusters_G.copy(), 
                                                                    num_merge=temp_degree)
                else: 
                    [new_samp_tab, new_clusters_g] = merge_clusters(samp_tab, 
                                                                    clusters_G.copy(), 
                                                                    num_merge=temp_degree)
                data_type = sub_folder.split("/")[2]
                final_output = os.path.join(output_path, data_type, "method_" + temp_method + "_degree_" + str(temp_degree), str(temp_trial))
                if os.path.exists(final_output) == False: 
                    os.makedirs(final_output)
                new_samp_tab.to_csv(os.path.join(final_output, "new_samp_tab.csv"))
                pickle.dump(new_clusters_g, open(os.path.join(final_output, "new_cluster_g.pickle"), "wb"))
                end_clusters = [x for x in new_clusters_g.nodes() if new_clusters_g.out_degree(x)==0]
                initial_clusters = [x for x in new_clusters_g.nodes() if new_clusters_g.in_degree(x)==0]
                lineage_cluster = onesc.extract_trajectory(new_clusters_g,initial_clusters , end_clusters)
                vector_thresh = onesc.find_threshold_vector(exp_tab, new_samp_tab, cluster_col = "leiden")
                lineage_time_change_dict = onesc.find_gene_change_trajectory(exp_tab, new_samp_tab, lineage_cluster, cluster_col, pt_col, vector_thresh, pseudoTime_bin=0.01)
                state_dict = onesc.define_states(exp_tab, new_samp_tab, lineage_cluster, vector_thresh, cluster_col)
                transition_dict = onesc.define_transition(state_dict)
                training_dict = onesc.curate_training_data(state_dict, transition_dict, lineage_time_change_dict, new_samp_tab, cluster_id = cluster_col, pt_id = pt_col, act_tolerance = 0.02)
                pickle.dump(training_dict, open(os.path.join(final_output, "training_data.pickle"), "wb"))
                start_time = time.time()
                corr_mat = onesc.calc_corr(exp_tab)
                temp_grn = onesc.create_network_ensemble(training_dict, 
                                                    corr_mat, 
                                                    ideal_edges = round(0.4 * corr_mat.shape[1]), 
                                                    num_generations = 300, 
                                                    max_iter = 30, 
                                                    num_parents_mating = 6, 
                                                    sol_per_pop = 30, 
                                                    reduce_auto_reg = True)
                time_lapse = time.time() - start_time
                pickle.dump(time_lapse, open(os.path.join(final_output, 'time_lapse.pickle'), 'wb'))
                temp_grn[0].to_csv(os.path.join(final_output, "OneSC_network.csv"))



