import pandas as pd 
import numpy as np 
import os
from pathlib import Path

datasets = ['myeloid_progenitors_full']

if os.path.isdir("../output/max_cross_correlation") == False: 
    os.mkdir("../output/max_cross_correlation")

for dataset in datasets: 
    print(dataset)
    methods_list = os.listdir("../output/cross_correlation/" + dataset)
    methods_list = [x for x in methods_list if 'DS_Store' not in x]
    for method in methods_list:
        print(method)
        for grn_style in ['density']:
            if grn_style == 'orphan' and method == 'run_OneCC_pyEpoch':
                continue
            load_path = "../output/cross_correlation/" + dataset + "/" + method + "/" + grn_style

            cross_corr_files = os.listdir(load_path)
            cross_corr_files = [x for x in cross_corr_files if "CrossCorrelation" in x]
            sample_id = 0
            big_df = pd.DataFrame()
            for cross_corr_file in cross_corr_files:
                temp_df = pd.read_csv(load_path + '/' + cross_corr_file, index_col=0)
                temp_df['sample_id'] = str(sample_id)
                sample_id = sample_id + 1
                big_df = pd.concat([big_df, temp_df])

            top_df = pd.DataFrame()
            for temp_lineage in np.unique(big_df['lineage_type']):
                print(temp_lineage)
                lineage_df = big_df.loc[big_df['lineage_type'] == temp_lineage, :]
                best_sample_id = ''
                best_sample_cc = -99
                for temp_sample_id in np.unique(lineage_df['sample_id']):
                    average_score = np.mean(lineage_df.loc[lineage_df['sample_id'] == temp_sample_id, 'cross_corr'])
                    if best_sample_cc < average_score:
                        best_sample_cc = average_score
                        best_sample_id = temp_sample_id
                top_df = pd.concat([top_df, lineage_df.loc[lineage_df['sample_id'] == best_sample_id, :]])

            save_path = "../output/max_cross_correlation/" + dataset + "/" + method + "/" + grn_style
            Path(save_path).mkdir( parents=True, exist_ok=True )
            top_df.to_csv(save_path + "/" + "max_cross_correlation.csv")
           