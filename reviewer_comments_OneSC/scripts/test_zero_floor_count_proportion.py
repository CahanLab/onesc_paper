import numpy as np 
import pandas as pd
import os
import pickle 

out_path = "../output/test_zero_floor_OneSC/myeloid_progenitors_full/ss_prop"
if os.path.isdir(out_path) == False: 
    os.makedirs(out_path)

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

data_type = 'myeloid_progenitors_full'
print(data_type)
states_dict = dict()
all_states_dict = pickle.load(open(os.path.join('../input/files_needed/myeloid_progenitors_full/states_dict.pickle'), "rb"))
ss_states = list()
for temp_index in all_states_dict.keys():
    states_tab = all_states_dict[temp_index]
    for temp_col in range(0, states_tab.shape[1]):
        string_pattern = to_string_pattern(states_tab.loc[:, states_tab.columns[temp_col]])
        if string_pattern in states_dict.keys():
            continue
        states_dict[string_pattern] = states_tab.columns[temp_col]
        ss_states.append(states_tab.columns[temp_col])
ss_states.append("other")

# this is to get the overexpression 
input_path = '../output/test_zero_floor_OneSC'
perturb_df = pd.DataFrame(data = 0, columns = ss_states, index = list(states_tab.index) + ['wild_type'])
perturb_type = 'overexpression'
all_genes = os.listdir(os.path.join(input_path, data_type, perturb_type))
all_genes = [x for x in all_genes if x != '.DS_Store']
all_genes = all_genes + ['wild_type']
for gene in all_genes:
    if gene == 'wild_type':
        sim_folder_path = '../output/test_zero_floor_OneSC/myeloid_progenitors_full/wildtype/'
        sim_files = [x for x in os.listdir(sim_folder_path) if '.csv' in x]
        raw_cluster_files = cluster_simulation(sim_folder_path)
        for temp_index in raw_cluster_files.keys():
            if temp_index in states_dict.keys():
                perturb_df.loc[gene, states_dict[temp_index]] = perturb_df.loc[gene, states_dict[temp_index]] + len(raw_cluster_files[temp_index])
            else:
                perturb_df.loc[gene, 'other'] = perturb_df.loc[gene, 'other'] + len(raw_cluster_files[temp_index])
    else:
        sim_folder_path = '../output/test_zero_floor_OneSC/myeloid_progenitors_full/' + "/" + perturb_type + "/" + gene + "/"
        sim_files = [x for x in os.listdir(sim_folder_path) if '.csv' in x]
        raw_cluster_files = cluster_simulation(sim_folder_path)
        for temp_index in raw_cluster_files.keys():
            if temp_index in states_dict.keys():
                perturb_df.loc[gene, states_dict[temp_index]] = perturb_df.loc[gene, states_dict[temp_index]] + len(raw_cluster_files[temp_index])
            else:
                perturb_df.loc[gene, 'other'] = perturb_df.loc[gene, 'other'] + len(raw_cluster_files[temp_index])
    
perturb_df.to_csv(os.path.join(out_path, data_type + "_" + perturb_type + "_proportion.csv"))

perturb_df = pd.DataFrame(data = 0, columns = ss_states, index = list(states_tab.index) + ['wild_type'])
perturb_type = 'knockout'
all_genes = os.listdir(os.path.join(input_path, data_type, perturb_type))
all_genes = [x for x in all_genes if x != '.DS_Store']
all_genes = all_genes + ['wild_type']
for gene in all_genes:
    if gene == 'wild_type':
        sim_folder_path = '../output/test_zero_floor_OneSC/myeloid_progenitors_full/wildtype/'
        sim_files = [x for x in os.listdir(sim_folder_path) if '.csv' in x]
        raw_cluster_files = cluster_simulation(sim_folder_path)
        for temp_index in raw_cluster_files.keys():
            if temp_index in states_dict.keys():
                perturb_df.loc[gene, states_dict[temp_index]] = perturb_df.loc[gene, states_dict[temp_index]] + len(raw_cluster_files[temp_index])
            else:
                perturb_df.loc[gene, 'other'] = perturb_df.loc[gene, 'other'] + len(raw_cluster_files[temp_index])
    else:
        sim_folder_path = '../output/test_zero_floor_OneSC/myeloid_progenitors_full/' + "/" + perturb_type + "/" + gene + "/"
        sim_files = [x for x in os.listdir(sim_folder_path) if '.csv' in x]
        raw_cluster_files = cluster_simulation(sim_folder_path)
        for temp_index in raw_cluster_files.keys():
            if temp_index in states_dict.keys():
                perturb_df.loc[gene, states_dict[temp_index]] = perturb_df.loc[gene, states_dict[temp_index]] + len(raw_cluster_files[temp_index])
            else:
                perturb_df.loc[gene, 'other'] = perturb_df.loc[gene, 'other'] + len(raw_cluster_files[temp_index])
perturb_df.to_csv(os.path.join(out_path, data_type + "_" + perturb_type + "_proportion.csv"))
