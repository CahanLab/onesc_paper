import pandas as pd
import numpy as np 
import pySingleCellNet as pySCN
import os
import scanpy as sc 
import pickle

if "run_pySCN" not in os.listdir("../output/"):
    os.mkdir("../output/run_pySCN")

def train_pySCN(data_type):
    if data_type not in os.listdir("../output/run_pySCN"):
        os.mkdir("../output/run_pySCN/" + data_type)
    target_dir = "../output/run_pySCN/" + data_type + '/'

    train_data = sc.read_h5ad("../output/BEELINE_clusters/" + data_type + "/redefined_adata.h5ad")
    train_ncells = round(np.min(train_data.obs['cluster_id'].value_counts()) * 0.8)

    expTrain, expVal = pySCN.splitCommonAnnData(train_data, ncells=train_ncells,dLevel="cluster_id")
    [cgenesA, xpairs, tspRF] = pySCN.scn_train_simple(expTrain, 
                                               nTopGenes = 100, 
                                               nRand = 0, 
                                               nTrees = 1000,
                                               nTopGenePairs = 21,  
                                               dLevel = "cluster_id")
    pickle.dump([cgenesA, xpairs, tspRF], open(target_dir + "training_obj.pickle", "wb"))

    [cgenesA, xpairs, tspRF] = pySCN.scn_train_simple(train_data, 
                                               nTopGenes = 100, 
                                               nRand = 0, 
                                               nTrees = 1000,
                                               nTopGenePairs = 5,  
                                               dLevel = "cluster_id")
    pickle.dump([cgenesA, xpairs, tspRF], open(target_dir + "full_training_obj.pickle", "wb"))
    
    return [cgenesA, xpairs, tspRF]

data_type_list = ['dyn-TF', 'dyn-BF', 'dyn-LI', 'dyn-CY', 'dyn-BFC', 'dyn-LI_Extended']

for data_type in data_type_list:
    print(data_type)
    [cgenesA, xpairs, tspRF] = train_pySCN(data_type)

for data_type in data_type_list:
    target_dir = "../output/run_pySCN/" + data_type + '/'

    [cgenesA, xpairs, tspRF] = pickle.load(open(target_dir + "training_obj.pickle", "rb"))
    if data_type == 'dyn-LI_Extended':
        data_type = 'dyn-LI'
    query_exp = pd.read_csv("../output/compiled_simulation/" + data_type + "_exp.csv", index_col = 0)
    query_st = pd.read_csv("../output/compiled_simulation/" + data_type + "_st.csv", index_col = 0)

    query_adata = adata = sc.AnnData(query_exp.T)
    query_adata.obs = query_st

    adVal = pySCN.scn_classify(query_adata, cgenesA, xpairs, tspRF, nrand = 0)

    class_mat = adVal.to_df()
    class_mat.to_csv(target_dir + "classification_mat.csv")

    samp_tab = adVal.obs
    samp_tab.to_csv(target_dir + "samp_tab.csv")





