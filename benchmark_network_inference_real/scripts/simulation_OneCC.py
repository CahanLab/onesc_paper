from ast import AsyncFunctionDef
import pandas as pd
import numpy as np 
import onesc
import networkx as nx

import seaborn as sns
import pickle
import matplotlib.pyplot as plt

import os
import scanpy as sc
import multiprocessing as mp

pool = mp.pool.ThreadPool(mp.cpu_count())

if os.path.isdir("../output/simulation_OneSC") == False:
    os.mkdir("../output/simulation_OneSC")

data_types = os.listdir("../output/curated_networks")
data_types = [x for x in data_types if 'DS_Store' not in x]

for data_type in data_types:
    state_dict = pickle.load(open("../Beeline_benchmark/run_OneSC/" + data_type + "/state_dict.pickle", "rb"))
    init_state = state_dict['trajectory_0'].iloc[:, 0]

    exp_dict = dict()
    for gene in init_state.index: 
        if init_state[gene] == 1:
            exp_dict[gene] = 2
        else:
            exp_dict[gene] = 0
    
    if os.path.isdir("../output/simulation_OneSC/" + data_type) == False:
        os.mkdir("../output/simulation_OneSC/" + data_type)

    GRN_methods = os.listdir("../output/curated_networks/" + data_type)
    GRN_methods = [x for x in GRN_methods if 'DS_Store' not in x]
    for GRN_method in GRN_methods:
        print(GRN_method)
        if os.path.isdir("../output/simulation_OneSC/" + data_type + "/" + GRN_method) == False:
            os.mkdir("../output/simulation_OneSC/" + data_type + "/" + GRN_method)
        for network_type in ['density']:
            if os.path.isdir("../output/simulation_OneSC/" + data_type + "/" + GRN_method + "/" + network_type) == False:
                os.mkdir("../output/simulation_OneSC/" + data_type + "/" + GRN_method + "/" + network_type)
            save_folder_path = "../output/simulation_OneSC/" + data_type + "/" + GRN_method + "/" + network_type
            if network_type == 'orphan':
                if GRN_method == 'run_OneSC':
                    continue
                train_grn = pd.read_csv("../output/curated_networks/" + data_type + "/" + GRN_method + "/curated_network_orphan.csv", sep = ',', index_col=0)
            else:
                train_grn = pd.read_csv("../output/curated_networks/" + data_type + "/" + GRN_method + "/curated_network_density.csv", sep = ',', index_col=0)
            MyNetwork = onesc.network_structure()
            MyNetwork.fit_grn(train_grn)
            pickle.dump(MyNetwork, open(save_folder_path + "/OneSC_subnetwork.pickle", "wb"))
            temp_simulator = onesc.OneSC_simulator()
            temp_simulator.add_network_compilation(GRN_method, MyNetwork)
            pickle.dump(temp_simulator, file = open(save_folder_path + "/OneSC_simulator.pickle", "wb"))
            num_samples = 200 # simulate 200 cells 
            num_runs = list(range(0, num_samples))     
            if os.path.isdir(save_folder_path + "/simulation") == False:
                os.mkdir(save_folder_path + "/simulation") 
            def run_parallel(i):
                np.random.seed(i)
                init_exp_dict = dict()
                for temp_gene in exp_dict.keys(): 
                    init_exp_dict[temp_gene] = exp_dict[temp_gene]
                temp_simulator.simulate_exp(init_exp_dict, GRN_method, num_sim = 1800, t_interval = 0.1, noise_amp = 0.5, random_seed = i)
                sim_exp = temp_simulator.sim_exp.copy()
                sim_exp.to_csv(save_folder_path + "/simulation/" + str(i) + "_simulated_exp.csv")   
            results = pool.map(run_parallel, num_runs) 

    print(data_type)



