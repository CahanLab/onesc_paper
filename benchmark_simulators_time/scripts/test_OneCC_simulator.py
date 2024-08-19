import pandas as pd
import numpy as np 
import onesc
import os
import time

folder_list = os.listdir("../output/random_networks/")
for temp_trial in folder_list: 
    for node_combo in os.listdir("../output/random_networks/" + temp_trial):

        if 'OneSC_time.txt' in os.listdir("../output/random_networks/" + temp_trial + "/" + node_combo):
            continue 

        grn_tab = pd.read_csv("../output/random_networks/" + temp_trial + "/" + node_combo + "/OneCC_network.txt", sep = '\t')
        MyNetwork = onesc.network_structure()
        MyNetwork.fit_grn(grn_tab)
        temp_simulator = onesc.OneSC_simulator()
        temp_simulator.add_network_compilation('OneSC', MyNetwork)
        
        exp_dict = dict() 
        for gene in np.unique(grn_tab['TF']): 
            exp_dict[gene] = 1

        num_samples = 5
        num_runs = list(range(0, num_samples))

        print("../output/random_networks/" + temp_trial + "/" + node_combo)
        if "OneSC_simulation" not in os.listdir("../output/random_networks/" + temp_trial + "/" + node_combo):
            os.mkdir("../output/random_networks/" + temp_trial + "/" + node_combo + "/OneSC_simulation/")

        def run_seq(i):
            np.random.seed(i)
            init_exp_dict = dict()
            for temp_gene in exp_dict.keys(): 
                #init_exp_dict[temp_gene] = max(0.02, exp_dict[temp_gene] + np.random.normal(0, 1))
                init_exp_dict[temp_gene] = exp_dict[temp_gene]
            temp_simulator.simulate_exp(init_exp_dict, "OneSC", num_sim = 5000, t_interval = 0.1, noise_amp = 0.1, random_seed = i)
            sim_exp = temp_simulator.sim_exp.copy()
            sim_exp.to_csv("../output/random_networks/" + temp_trial + "/" + node_combo + "/OneSC_simulation" + "/" + str(i) + "_simulated_exp.csv")

        start = time.time()
        for i in num_runs: 
            run_seq(i) 
        end = time.time()

        file1 = open("../output/random_networks/" + temp_trial + "/" + node_combo + '/OneSC_time.txt', 'w')
        L = 'user time is ' + str(end - start)
        file1.write(L)
        file1.close()
