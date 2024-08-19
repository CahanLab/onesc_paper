import pandas as pd 
import numpy as np 
import OneCC_bool
import os
import seaborn as sns
import pickle 

# make the UMAP overlaps of the 
if 'overlap_UMAP' not in os.listdir("../output"):
    os.mkdir("../output/overlap_UMAP")

for data_type in ['dyn-BF', 'dyn-BFC', 'dyn-TF', 'dyn-LI', 'dyn-CY']:
    print(data_type)
    train_file_path = "../BEELINE-data/inputs/Synthetic/" + data_type + "/" + data_type + "-5000-1/"
    train_exp = pd.read_csv(train_file_path + "ExpressionData.csv", index_col=0)
    UMAP_obj = OneCC_bool.UMAP_embedding_train(train_exp)

    pickle.dump(UMAP_obj, open("../output/overlap_UMAP/" + data_type + "_UMAP_obj.pickle", "wb"))
    train_st = pd.read_csv(train_file_path + "PseudoTime.csv", index_col = 0)
    
    train_st = train_st.fillna(0)
    train_st['pseudoTime'] = train_st.sum(axis = 1)
    train_UMAP_coord = OneCC_bool.UMAP_embedding_apply(UMAP_obj, train_exp)
    train_UMAP_coord['sim_time'] = train_st['pseudoTime']
    train_UMAP_coord.to_csv("../output/overlap_UMAP/" + data_type + "_train_UMAP.csv")

    query_exp_path = "../output/compiled_simulation/" + data_type + "_exp.csv"
    query_exp = pd.read_csv(query_exp_path, index_col=0)

    query_st_path = "../output/compiled_simulation/" + data_type + "_st.csv"
    query_st = pd.read_csv(query_st_path, index_col=0)

    query_UMAP_coord = OneCC_bool.UMAP_embedding_apply(UMAP_obj, query_exp)
    query_UMAP_coord['sim_time'] = query_st['sim_time']

    query_UMAP_coord.to_csv("../output/overlap_UMAP/" + data_type + "_query_UMAP.csv")
    
    sns_plot = sns.scatterplot('UMAP_1', 'UMAP_2', data=query_UMAP_coord, hue='sim_time', alpha = 0.5)
    #sns_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    sns_plot.figure.savefig("../output/overlap_UMAP/" + data_type + "_OneCC_UMAP.png")
    sns_plot.clear()

    sns_plot = sns.scatterplot('UMAP_1', 'UMAP_2', data=train_UMAP_coord, hue='sim_time', alpha = 0.5)
    #sns_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    sns_plot.figure.savefig("../output/overlap_UMAP/" + data_type + "_BoolODE_UMAP.png")
    sns_plot.clear()

# run once for LI extended for BoolODE 
data_type = 'dyn-LI_Extended'

train_file_path = "../BoolODE/Debug/dyn-linear-1_extended/dyn-linear-1_extended-5000-1/"
train_exp = pd.read_csv(train_file_path + "ExpressionData.csv", index_col=0)
UMAP_obj = OneCC_bool.UMAP_embedding_train(train_exp)

pickle.dump(UMAP_obj, open("../output/overlap_UMAP/" + data_type + "_UMAP_obj.pickle", "wb"))
train_st = pd.read_csv(train_file_path + "PseudoTime.csv", index_col = 0)

train_st = train_st.fillna(0)
train_st['pseudoTime'] = train_st.sum(axis = 1)
train_UMAP_coord = OneCC_bool.UMAP_embedding_apply(UMAP_obj, train_exp)
train_UMAP_coord['sim_time'] = train_st['pseudoTime']
train_UMAP_coord.to_csv("../output/overlap_UMAP/" + data_type + "_train_UMAP.csv")

query_exp_path = "../output/compiled_simulation/" + 'dyn-LI' + "_exp.csv"
query_exp = pd.read_csv(query_exp_path, index_col=0)

query_st_path = "../output/compiled_simulation/" + 'dyn-LI' + "_st.csv"
query_st = pd.read_csv(query_st_path, index_col=0)

query_UMAP_coord = OneCC_bool.UMAP_embedding_apply(UMAP_obj, query_exp)
query_UMAP_coord['sim_time'] = query_st['sim_time']

query_UMAP_coord.to_csv("../output/overlap_UMAP/" + data_type + "_query_UMAP.csv")

query_UMAP_coord = pd.read_csv("../output/overlap_UMAP/" + data_type + "_query_UMAP.csv", index_col = 0)
sns_plot = sns.scatterplot('UMAP_1', 'UMAP_2', data=query_UMAP_coord, hue='sim_time', alpha = 0.5)
#sns_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
sns_plot.figure.savefig("../output/overlap_UMAP/" + data_type + "_OneCC_UMAP.png")
sns_plot.clear()

train_UMAP_coord = pd.read_csv("../output/overlap_UMAP/" + data_type + "_train_UMAP.csv", index_col = 0)
sns_plot = sns.scatterplot('UMAP_1', 'UMAP_2', data=train_UMAP_coord, hue='sim_time', alpha = 0.5)
#sns_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
sns_plot.figure.savefig("../output/overlap_UMAP/" + data_type + "_BoolODE_UMAP.png")
sns_plot.clear()