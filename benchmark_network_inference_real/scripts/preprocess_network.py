import pandas as pd
import numpy as np 

def curate_network_orphan(GRN_style = 'BEELINE', data_type = 'dyn-BF'):
    if GRN_style in ['run_OneCC_full', 'run_OneCC_pyEpoch', 'run_OneCC_pyEpoch_thresh', 'run_OneCC_regulators']:
        inferred_grn = pd.read_csv(f"../Beeline_benchmark/{GRN_style}/{data_type}/OneCC_network.csv")
        inferred_grn = inferred_grn.drop(inferred_grn.columns[0],axis=1).copy()
        return inferred_grn
    elif GRN_style in ['GENIE3', 'GRISLI', 'GRNBOOST2', 'LEAP', 'PIDC', 'PPCOR', 'SCRIBE', 'SINGE']:
        inferred_grn = pd.read_csv(f"../Beeline_benchmark/outputs/real_data/{GRN_style}/{data_type}/rankedEdges.csv", sep = '\t')
        if inferred_grn.shape[0] == 0:
            return inferred_grn
        