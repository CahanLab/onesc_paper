import pandas as pd 
import numpy as np 
import onesc
import os
import multiprocessing as mp

folder_list = os.listdir("../BEELINE-data/inputs/Synthetic")
folder_list = [x for x in folder_list if 'dyn' in x]
folder_list = [x for x in folder_list if 'LL' not in x]
# simulate maybe 100 runs for each trajectory
# probably let it run overnight or something 

init_dict = dict()
#init_dict['dyn-LL'] = ['g1']
init_dict['dyn-LI'] = ['g1']
init_dict['dyn-BFC'] = ['g1']
init_dict['dyn-BF'] = ['g1']
init_dict['dyn-CY'] = ['g1', 'g2', 'g3']
init_dict['dyn-TF'] = ['g1']


if "OneSC_sim" not in os.listdir("../output"):
    os.mkdir("../output/OneSC_sim")

pool = mp.pool.ThreadPool(mp.cpu_count() - 1)
for temp_folder in folder_list: 
    init_cond = init_dict[temp_folder]
    grn_tab = pd.read_csv("../BEELINE-data/inputs/Synthetic/" + temp_folder + "/" + temp_folder + "-2000-1/refNetwork.csv")
    grn_tab.columns = ['TF', 'TG', 'Type']

    exp_tab = pd.read_csv("../BEELINE-data/inputs/Synthetic/" + temp_folder + "/" + temp_folder + "-2000-1/ExpressionData.csv", index_col=0)

    MyNetwork = onesc.network_structure()
    MyNetwork.fit_grn(grn_tab)

    temp_simulator = onesc.OneSC_simulator()
    temp_simulator.add_network_compilation('OneSC', MyNetwork)
    
    exp_dict = dict() 
    for gene in np.unique(exp_tab.index): 
        if gene in init_cond:
            exp_dict[gene] = 2
        else:
            exp_dict[gene] = 0

    TFs = np.unique(grn_tab['TF'])
    temp_simulator.TFs = TFs 

    num_samples = 100
    num_runs = list(range(0, num_samples))

    if temp_folder not in os.listdir("../output/OneSC_sim/"):
        os.mkdir("../output/OneSC_sim/" + temp_folder)
    
    target_dir = "../output/OneSC_sim/" + temp_folder

    def run_parallel(i):
        init_exp_dict = dict()
        for temp_gene in exp_dict.keys(): 
            init_exp_dict[temp_gene] = exp_dict[temp_gene]
        temp_simulator.simulate_exp(init_exp_dict, "OneSC", num_sim = 2500, t_interval = 0.05, noise_amp = 0.1, random_seed = i)
        sim_exp = temp_simulator.sim_exp.copy()
        sim_exp.to_csv("../output/OneSC_sim/" + temp_folder + "/" + str(i) + "_simulated_exp.csv")
        return None 
    
    results = pool.map(run_parallel, num_runs)

pool.close()