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
subfolders = [x.split("/")[2] for x in subfolders]

np.linspace(5, 30, 6)

total_max_iter = list()
total_num_gen = list()
total_avg_time = list()

for temp_max_iter in np.linspace(5, 30, 6):
    for temp_num_gen in np.linspace(10, 100, 10):
        time_list = list()
        for data_type in subfolders: 
            new_folder = 'numgeneration_' + str(int(temp_num_gen)) + "_maxiter_" + str(int(temp_max_iter)) + "/" + data_type
            output_path = os.path.join("../output/screen_generations_finer", new_folder)
            time_duration = pickle.load(open(os.path.join(output_path, "time_lapse.pickle"), "rb"))
            time_list.append(time_duration)
        avg_time = np.mean(time_list)
        total_max_iter.append(temp_max_iter)
        total_num_gen.append(temp_num_gen)
        total_avg_time.append(avg_time)

big_df = pd.DataFrame()
big_df['max_iter'] = total_max_iter
big_df['num_gen'] = total_num_gen
big_df['mean_time'] = total_avg_time

big_df.to_csv(os.path.join("../output/screen_generations_finer/generation_time_big_df.csv"))