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

if os.path.exists("../output/screen_density") == False:
    os.makedirs("../output/screen_density")

for temp_density in np.linspace(0,1, 11):
    print(temp_density)
    if os.path.exists(os.path.join("../output/screen_density", str(round(temp_density, 2)))) == False:
        os.makedirs(os.path.join("../output/screen_density", str(round(temp_density, 2))))
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
                                            ideal_edges = round(temp_density * corr_mat.shape[1]), 
                                            num_generations = 300, 
                                            max_iter = 30, 
                                            num_parents_mating = 6, 
                                            sol_per_pop = 30, 
                                            reduce_auto_reg = True)


        output_path = os.path.join("../output/screen_density", str(round(temp_density, 2)), data_type.split("/")[2])
        if os.path.exists(output_path) == False:
            os.makedirs(output_path)

        time_lapse = time.time() - start_time
        pickle.dump(time_lapse, open(os.path.join(output_path, 'time_lapse.pickle'), 'wb'))
        print("--- %s seconds ---" % (time_lapse))

        temp_grn[0].to_csv(os.path.join(output_path, "OneSC_network.csv"))

