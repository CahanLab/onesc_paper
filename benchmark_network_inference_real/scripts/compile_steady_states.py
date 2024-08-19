from statistics import mean
import pandas as pd 
import numpy as np 
import os
import scanpy as sc
import statsmodels.api as sm
import pickle
import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
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

datasets = ['myeloid_progenitors_full']   
for dataset in datasets: 
    print(dataset)
    methods_list = os.listdir("../output/simulation_OneSC/" + dataset)
    methods_list = [x for x in methods_list if 'DS_Store' not in x]
    for method in methods_list:
        print(method)
        for grn_style in ['density']:
            if grn_style == 'orphan' and method == 'run_OneSC':
                continue
            load_path = "../output/simulation_OneSC/" + dataset + "/" + method + "/" + grn_style + "/simulation/"
            unique_ss = find_all_steady_states(load_path)
            save_path = "../output/steady_states/" + dataset + "/" + method + "/" + grn_style
            Path(save_path).mkdir( parents=True, exist_ok=True )
            unique_ss.to_csv(save_path + "/" + "unique_steady_states.csv")

for dataset in datasets: 
    print(dataset)
    methods_list = os.listdir("../output/simulation_OneSC/" + dataset)
    methods_list = [x for x in methods_list if 'DS_Store' not in x]
    steady_states = pickle.load(open("../Beeline_benchmark/run_OneSC/" + dataset + "/state_dict.pickle", "rb"))

    ss_list = []
    for temp_lineage in steady_states.keys(): 
        ss_list.append(steady_states[temp_lineage].columns[steady_states[temp_lineage].shape[1] - 1])
    
    benchmark_df = pd.DataFrame(data = 0, index=ss_list + ['other'], columns=methods_list)

    for method in methods_list:
        print(method)

        for grn_style in ['density']:
            if grn_style == 'orphan' and method == 'run_OneSC':
                continue
            load_path = "../output/simulation_OneSC/" + dataset + "/" + method + "/" + grn_style + "/simulation/"
            state_dict = cluster_simulation(load_path)
            all_files = os.listdir(load_path)

            for temp_lineage in steady_states.keys(): 
                ss_string = to_string_pattern(steady_states[temp_lineage].iloc[:, steady_states[temp_lineage].shape[1] - 1])
                ss_name = steady_states[temp_lineage].columns[steady_states[temp_lineage].shape[1] - 1]
                if ss_string in state_dict.keys(): 
                    benchmark_df.loc[ss_name, method] = len(state_dict[ss_string])
            
            benchmark_df.loc['other', method] = len(all_files) - np.sum(benchmark_df[method])

            save_path = "../output/steady_states/" + dataset + "/" + method + "/" + grn_style
            benchmark_df.to_csv(save_path + "/" + 'ss_proportion.csv')


