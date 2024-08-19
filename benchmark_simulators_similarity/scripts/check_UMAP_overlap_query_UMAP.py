import pandas as pd 
import numpy as np 
import OneCC_bool
import os
import seaborn as sns
import pickle 

# make the UMAP overlaps of the 
if 'overlap_UMAP_query_UMAP' not in os.listdir("../output"):
    os.mkdir("../output/overlap_UMAP_query_UMAP")

for data_type in ['dyn-BF', 'dyn-BFC', 'dyn-TF', 'dyn-LI', 'dyn-CY']:

    query_exp_path = "../output/compiled_simulation/" + data_type + "_exp.csv"
    query_exp = pd.read_csv(query_exp_path, index_col=0)

    query_st_path = "../output/compiled_simulation/" + data_type + "_st.csv"
    query_st = pd.read_csv(query_st_path, index_col=0)
    query_st = query_st.iloc[::20, :]

    query_exp = query_exp.loc[:, query_st.index]
    UMAP_obj = OneCC_bool.UMAP_embedding_train(query_exp)
    pickle.dump(UMAP_obj, open("../output/overlap_UMAP_query_UMAP/" + data_type + "_UMAP_obj.pickle", "wb"))

    query_UMAP_coord = OneCC_bool.UMAP_embedding_apply(UMAP_obj, query_exp)
    query_UMAP_coord['sim_time'] = query_st['sim_time']

    query_UMAP_coord.to_csv("../output/overlap_UMAP_query_UMAP/" + data_type + "_query_UMAP.csv")
    
    sns_plot = sns.scatterplot('UMAP_1', 'UMAP_2', data=query_UMAP_coord, hue='sim_time', alpha = 0.5)
    #sns_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    sns_plot.figure.savefig("../output/overlap_UMAP_query_UMAP/" + data_type + "_OneCC_UMAP.png")
    sns_plot.clear()
