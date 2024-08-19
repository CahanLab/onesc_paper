import pandas as pd
import numpy as np 
import OneCC_bool
import networkx as nx

import seaborn as sns
import pickle
import matplotlib.pyplot as plt

import os
import scanpy as sc
import multiprocessing as mp

import matplotlib.pyplot as plt

input_path = "../output/simulation_perturbation/myeloid_progenitors_full/"
perturb_types = os.listdir(input_path)
perturb_types = [x for x in perturb_types if 'DS_Store' not in x]
perturb_types = [x for x in perturb_types if 'OneCC_simulator' not in x]

out_root_path = "../output/plot_UMAPs"
if os.path.isdir(out_root_path) == False: 
    os.makedirs(out_root_path)

def find_all_steady_states(sim_folder_path): 
    all_files = os.listdir(sim_folder_path)
    all_files = [x for x in all_files if 'simulated_exp.csv' in x]

    steady_states_df = pd.DataFrame()
    for exp_file in all_files: 
        sample_id = exp_file.split("_")[0]
        sim_exp = pd.read_csv(sim_folder_path + exp_file, index_col=0)
        ss_profile = sim_exp.iloc[:, -1]
        ss_profile = np.array(ss_profile >= 1)
        ss_profile = ss_profile.astype(int)
        steady_states_df[sample_id] = ss_profile
        steady_states_df.index = sim_exp.index

    steady_states_df = steady_states_df.T
    steady_states_df = steady_states_df.drop_duplicates()

    return steady_states_df.T

def cluster_simulation(sim_folder_path): 
    all_files = os.listdir(sim_folder_path)
    all_files = [x for x in all_files if 'simulated_exp.csv' in x]
    steady_states_dict = dict()
    for exp_file in all_files: 
        sim_exp = pd.read_csv(sim_folder_path + exp_file, index_col=0)
        ss_profile = sim_exp.iloc[:, -1]
        ss_profile = np.array(ss_profile >= 1)
        ss_profile = ss_profile.astype(int)
        ss_profile = ss_profile.astype(str)
        ss_profile = "_".join(ss_profile)
        if ss_profile not in steady_states_dict.keys(): 
            #steady_states_dict[ss_profile] = [sim_exp]
            steady_states_dict[ss_profile] = [exp_file]
        else:
            #steady_states_dict[ss_profile].append(sim_exp)
            steady_states_dict[ss_profile].append(exp_file)
    return steady_states_dict

def to_string_pattern(gene_pattern):
    gene_pattern = [str(x) for x in gene_pattern]
    return "_".join(gene_pattern)

for perturb_type in perturb_types: 
    print(perturb_type)
    TFs_list = os.listdir(os.path.join(input_path, perturb_type))
    TFs_list = [x for x in TFs_list if 'DS_Store' not in x]
    for temp_TF in TFs_list: 
        out_path = os.path.join(out_root_path, perturb_type, temp_TF) 
        if os.path.isdir(out_path) == False: 
            os.makedirs(out_path)

        sim_files = [x for x in os.listdir(os.path.join(input_path, perturb_type, temp_TF)) if '.csv' in x]
        raw_cluster_files = cluster_simulation(os.path.join(input_path, perturb_type, temp_TF + "/"))

        states_dict = dict()
        for temp_file in os.listdir(os.path.join("../output/files_needed/myeloid_progenitors_full/extract_states")):
            states_tab = pd.read_csv(os.path.join('../output/files_needed/myeloid_progenitors_full/extract_states', temp_file), index_col = 0)
            string_pattern = to_string_pattern(states_tab.loc[:, states_tab.columns[-1]])
            states_dict[string_pattern] = states_tab.columns[-1]

        # map all simulation samples to a end type 
        sample_ss_dict = dict()
        for temp_string in raw_cluster_files.keys():
            for temp_sample in raw_cluster_files[temp_string]:
                if temp_string in states_dict.keys():
                    sample_ss_dict[temp_sample.replace("_exp.csv", "")] = states_dict[temp_string]
                else: 
                    sample_ss_dict[temp_sample.replace("_exp.csv", "")] = 'other'
                    
        big_sim_df = pd.DataFrame()
        for sim_file in sim_files: 
            experiment_title = sim_file.replace("_exp.csv", "")
            temp_sim = pd.read_csv(os.path.join(input_path, perturb_type, temp_TF + "/") + sim_file, index_col = 0)
            temp_sim = temp_sim[temp_sim.columns[::50]] # probably have to do it every 10 cells 
            temp_sim.columns = experiment_title + "-" + temp_sim.columns
            big_sim_df = pd.concat([big_sim_df, temp_sim], axis = 1)

        train_obj = pickle.load(open("../output/files_needed/myeloid_progenitors_full/train_obj.pickle", "rb"))
        UMAP_coord = OneCC_bool.UMAP_embedding_apply(train_obj, big_sim_df)
        UMAP_coord['sim_time'] = [int(x.split("-")[1]) for x in list(UMAP_coord.index)]
        UMAP_coord['sample'] = [str(x.split("-")[0]) for x in list(UMAP_coord.index)]
        UMAP_coord['steady_states'] = None
        for temp_sample in np.unique(UMAP_coord['sample']):
            UMAP_coord.loc[UMAP_coord['sample'] == temp_sample, "steady_states"] = sample_ss_dict[temp_sample]
        UMAP_coord = pd.concat([UMAP_coord, big_sim_df.T], axis = 1)

        big_sim_df.to_csv(os.path.join(out_path, "big_sim_df.csv"))
        UMAP_coord.to_csv(os.path.join(out_path, 'UMAP_coord.csv'))


