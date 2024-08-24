import scanpy as sc
import pandas as pd
import numpy as np 
import onesc

import networkx as nx
import seaborn as sns
import pickle
import matplotlib.pyplot as plt

import os
import time

import math

import warnings
warnings.filterwarnings("ignore")

input_path = "../input/preprocessed_data"

subfolders = [ f.path for f in os.scandir(input_path) if f.is_dir() ]
subfolders = [x.strip("./") for x in subfolders]
subfolders = [x for x in subfolders if "ipynb" not in x]

if os.path.exists("../output/screen_generations_finer") == False:
    os.makedirs("../output/screen_generations_finer")

for temp_max_iter in np.linspace(5, 30, 6):
    for temp_num_gen in np.linspace(10, 100, 10):
        new_folder = 'numgeneration_' + str(int(temp_num_gen)) + "_maxiter_" + str(int(temp_max_iter))
        if os.path.exists(os.path.join("../output/screen_generations_finer", new_folder)) == False:
            os.makedirs(os.path.join("../output/screen_generations_finer", new_folder))
        for data_type in subfolders: 
            print(data_type)
            [training_dict, initial_clusters] = pickle.load(open('../' + data_type + "/training_package.pickle", "rb"))
            adata = sc.read_h5ad('../' + data_type + "/redefined_adata.h5ad")
            exp_train = adata.to_df()
            exp_train = exp_train.T

            start_time = time.time()
            corr_mat = onesc.calc_corr(exp_train)
            temp_grn = onesc.create_network_ensemble(training_dict, 
                                                corr_mat, 
                                                ideal_edges = round(0.4 * corr_mat.shape[1]), 
                                                num_generations = int(temp_num_gen), 
                                                max_iter = int(temp_max_iter), 
                                                num_parents_mating = 6, 
                                                sol_per_pop = 30, 
                                                reduce_auto_reg = True)


            output_path = os.path.join("../output/screen_generations_finer", new_folder, data_type.split("/")[2])
            if os.path.exists(output_path) == False:
                os.makedirs(output_path)

            time_lapse = time.time() - start_time
            pickle.dump(time_lapse, open(os.path.join(output_path, 'time_lapse.pickle'), 'wb'))
            print("--- %s seconds ---" % (time_lapse))

            temp_grn[0].to_csv(os.path.join(output_path, "OneSC_network.csv"))

