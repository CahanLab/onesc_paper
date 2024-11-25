import pandas as pd
import numpy as np 
import onesc

import pickle

import os
import multiprocessing as mp

pool = mp.pool.ThreadPool(mp.cpu_count() - 1)

output_path = "../output/test_zero_floor_OneSC"
if os.path.isdir(output_path) == False:
    os.mkdir(output_path)

data_types = os.listdir("../input/files_needed")
data_types = [x for x in data_types if os.path.isdir(os.path.join("../output/files_needed", x)) == True]

data_type = 'myeloid_progenitors_full'

state_dict = pickle.load(open("../input/files_needed/" + data_type + "/states_dict.pickle", "rb"))
init_state = state_dict['trajectory_0'].iloc[:, 0]
exp_dict = dict()
for gene in init_state.index: 
    if init_state[gene] == 1:
        exp_dict[gene] = 2
    else:
        exp_dict[gene] = 0

train_grn = pd.read_csv("../input/files_needed/" + data_type +  "/curated_network.csv", sep = ',', index_col=0)

MyNetwork = onesc.network_structure()
MyNetwork.fit_grn(train_grn, min_val = 0)

temp_simulator = onesc.OneSC_simulator()
temp_simulator.add_network_compilation('Myeloid', MyNetwork)
TFs = np.unique(train_grn['TF'])
temp_simulator.TFs = TFs 
if os.path.isdir(os.path.join(output_path, data_type)) == False:
    os.makedirs(os.path.join(output_path, data_type))

pickle.dump(temp_simulator, file = open(os.path.join(output_path, data_type, "OneCC_simulator.pickle"), "wb"))

# this is to generate a wildtype from the data 
save_folder_path = os.path.join(output_path, data_type, 'wildtype')
os.makedirs(save_folder_path, exist_ok=True)
num_samples = 200 
num_runs = list(range(0, num_samples))     
def run_parallel(i):
    np.random.seed(i)
    init_exp_dict = dict()
    for temp_gene in exp_dict.keys(): 
        #init_exp_dict[temp_gene] = max(0.02, exp_dict[temp_gene] + np.random.normal(0, 1))
        init_exp_dict[temp_gene] = exp_dict[temp_gene]
    temp_simulator.simulate_exp(init_exp_dict, "Myeloid", num_sim = 1800, t_interval = 0.1, noise_amp = 0.5, random_seed = i)
    sim_exp = temp_simulator.sim_exp.copy()
    sim_exp.to_csv(save_folder_path + "/" + str(i) + "_simulated_exp.csv")   
results = pool.map(run_parallel, num_runs) 

# this is to generate overexpression 
save_folder_path = os.path.join(output_path, data_type, 'overexpression')
for temp_TF in TFs:
    print(temp_TF)
    if os.path.isdir(save_folder_path + "/" + temp_TF) == False:
        os.makedirs(save_folder_path + "/" + temp_TF) 
    num_samples = 200 # TODO make sure you change this to 50
    num_runs = list(range(0, num_samples))     
    def run_parallel(i):
        np.random.seed(i)
        init_exp_dict = dict()
        for temp_gene in exp_dict.keys(): 
            init_exp_dict[temp_gene] = exp_dict[temp_gene]
        perturb_dict = dict()
        perturb_dict[temp_TF] = 3
        temp_simulator.simulate_exp(init_exp_dict, "Myeloid", perturb_dict = perturb_dict, num_sim = 1800, t_interval = 0.1, noise_amp = 0.5, random_seed = i)
        sim_exp = temp_simulator.sim_exp.copy()
        sim_exp.to_csv(save_folder_path + "/" + temp_TF + "/" + str(i) + "_simulated_exp.csv")   
    results = pool.map(run_parallel, num_runs) 

# this is to generate knockout 
save_folder_path = os.path.join(output_path, data_type, 'knockout')
for temp_TF in TFs:
    print(temp_TF)
    if os.path.isdir(save_folder_path + "/" + temp_TF) == False:
        os.makedirs(save_folder_path + "/" + temp_TF) 
    num_samples = 200 
    num_runs = list(range(0, num_samples))     
    def run_parallel(i):
        np.random.seed(i)
        init_exp_dict = dict()
        for temp_gene in exp_dict.keys(): 
            init_exp_dict[temp_gene] = exp_dict[temp_gene]
        perturb_dict = dict()
        perturb_dict[temp_TF] = -3
        temp_simulator.simulate_exp(init_exp_dict, "Myeloid", perturb_dict = perturb_dict, num_sim = 1800, t_interval = 0.1, noise_amp = 0.5, random_seed = i)
        sim_exp = temp_simulator.sim_exp.copy()
        sim_exp.to_csv(save_folder_path + "/" + temp_TF + "/" + str(i) + "_simulated_exp.csv")   
    results = pool.map(run_parallel, num_runs) 
