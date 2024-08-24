from statistics import mean
from turtle import forward
import pandas as pd 
import numpy as np 
import os
import scanpy as sc
import statsmodels.api as sm
import onesc

def cluster_simulation(sim_folder_path): 
    all_files = os.listdir(sim_folder_path)
    all_files = [x for x in all_files if 'simulated_exp.csv' in x]

    steady_states_dict = dict()
    for exp_file in all_files: 
        sim_exp = pd.read_csv(sim_folder_path + exp_file, index_col=0)
        ss_profile = sim_exp.iloc[:, -1]
        ss_profile = np.array(ss_profile >= 1)
        ss_profile = ss_profile.astype(int)
        ss_profile = ss_profile.astype(str)
        ss_profile = "_".join(ss_profile)
        if ss_profile not in steady_states_dict.keys(): 
            #steady_states_dict[ss_profile] = [sim_exp]
            steady_states_dict[ss_profile] = [exp_file]
        else:
            #steady_states_dict[ss_profile].append(sim_exp)
            steady_states_dict[ss_profile].append(exp_file)
    return steady_states_dict

def curate_average(sim_folder_path, curated_dict):
    mean_curated_dict = dict()
    for ss in curated_dict.keys():
        simulated_exp_list = curated_dict[ss]
        
        mean_sim_exp = pd.DataFrame()
        for sim_exp_file in simulated_exp_list: 
            sim_exp = pd.read_csv(sim_folder_path + sim_exp_file, index_col=0)
            if mean_sim_exp.shape[0] == 0: 
                mean_sim_exp = sim_exp
            else: 
                mean_sim_exp = mean_sim_exp + sim_exp
        mean_sim_exp = mean_sim_exp / len(simulated_exp_list)
        mean_curated_dict[ss] = mean_sim_exp
    return mean_curated_dict

def ScaleData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def find_cut_off(exp_tab, ss_string):
    ss_list = ss_string.split("_")
    boolean_exp_tab = exp_tab > 1.5
    boolean_exp_tab = boolean_exp_tab.astype("int")
    
    all_index = list()
    for temp_index in boolean_exp_tab.columns: 
        if np.array_equal(np.array(boolean_exp_tab[temp_index]), np.array(ss_list).astype('int')): 
            all_index.append(temp_index)
    return int(all_index[0]) + 100

def bin_smooth(time_series, pseudoTime_bin, smooth_style = "mean", spline_ME = 0.1):
    curr_time =  np.min(time_series['PseudoTime'])
    time_list = list()
    smoothed_exp = list()
    stand_dev = list()

    while curr_time < np.max(time_series['PseudoTime']): 
        temp_time_series = time_series.loc[np.logical_and(time_series['PseudoTime'] >= curr_time, time_series['PseudoTime'] < curr_time + pseudoTime_bin), :]
        
        # in the case of 0 expression just move on to the next time frame and move on 
        if temp_time_series.shape[0] == 0:
            #smoothed_exp.append(len(smoothed_exp) - 1)
            #time_list.append(curr_time)
            curr_time = curr_time + pseudoTime_bin
            continue

        time_list.append(curr_time)
        if smooth_style == 'mean':
            smoothed_exp.append(temp_time_series['expression'].mean())
        elif smooth_style == "median":
            smoothed_exp.append(temp_time_series['expression'].median())

        curr_time = curr_time + pseudoTime_bin
        stand_dev.append(temp_time_series['expression'].std())

    # spline_w = np.divide(1, stand_dev)

    smoothed_data = pd.DataFrame()
    smoothed_data['PseudoTime'] = time_list
    smoothed_data['expression'] = smoothed_exp

    #spline_s = smoothed_data.shape[0] * (spline_ME ** 2)
    #spline_xy = UnivariateSpline(smoothed_data['PseudoTime'],smoothed_data['expression'], s = spline_s)
    #moothed_data['splined_exp'] = spline_xy(smoothed_data['PseudoTime'])
    return smoothed_data 

# input training data exp -- specific lineage
# input training data sample table -- specific lineage 
# input simulation data exp -- specific lineage 

def calc_cross_correlation(sub_sim_tab, exp_tab, samp_tab, pt_col, resolution = 0.01, top_X = 60):
    sim_st = pd.DataFrame()
    sim_st['sim_time'] = sub_sim_tab.columns
    sim_st['scaled_pt'] = ScaleData(sim_st['sim_time'].astype('int'))

    sub_st = samp_tab
    sub_exp = exp_tab
    sub_st['scaled_pt'] = ScaleData(sub_st[pt_col])
    
    temp_df = pd.DataFrame()
    temp_df['genes'] = np.array(sub_exp.index)
    temp_df['cross_corr'] = None
    temp_df['lag'] = None
    temp_df.index = sub_exp.index

    for temp_gene in sub_exp.index:
        time_series_real = pd.DataFrame()
        time_series_real['expression'] = sub_exp.loc[temp_gene, :]
        time_series_real['PseudoTime'] = sub_st['scaled_pt']

        time_series_sim = pd.DataFrame()
        time_series_sim['expression'] = sub_sim_tab.loc[temp_gene, :]
        time_series_sim['PseudoTime'] = np.array(sim_st['scaled_pt'])

        real_smoothed = bin_smooth(time_series_real, pseudoTime_bin = resolution, smooth_style = "mean", spline_ME = 0.1)
        sim_smoothed = bin_smooth(time_series_sim, resolution, smooth_style = "mean")
        
        intersect_time_point = np.intersect1d(real_smoothed['PseudoTime'], sim_smoothed['PseudoTime'])
        real_smoothed = real_smoothed.loc[real_smoothed['PseudoTime'].isin(intersect_time_point), :]
        sim_smoothed = sim_smoothed.loc[sim_smoothed['PseudoTime'].isin(intersect_time_point), :]
        
        forward_corr = sm.tsa.stattools.ccf(real_smoothed['expression'], sim_smoothed['expression'])
        forward_corr = forward_corr[0:top_X]
        
        backward_corr = sm.tsa.stattools.ccf(sim_smoothed['expression'], real_smoothed['expression'])
        backward_corr = backward_corr[0:top_X]

        if np.max(forward_corr) >= np.max(backward_corr):
            temp_df.loc[temp_gene, 'cross_corr'] = np.max(forward_corr)
            lag_index = np.where(forward_corr == np.max(forward_corr))[0][0]
            lag = real_smoothed.loc[lag_index, 'PseudoTime'] + real_smoothed.loc[0, 'PseudoTime']
            temp_df.loc[temp_gene, 'lag'] = lag
        else:
            temp_df.loc[temp_gene, 'cross_corr'] = np.max(forward_corr)
            lag_index = np.where(forward_corr == np.max(forward_corr))[0][0]
            lag = real_smoothed.loc[lag_index, 'PseudoTime'] + real_smoothed.loc[0, 'PseudoTime']
            temp_df.loc[temp_gene, 'lag'] = -lag
    return temp_df

def calc_celltype(sub_sim_tab, sim_folder_path):
    parent_folder = sim_folder_path.split("/")[0]
    lineage_files = os.listdir(parent_folder)
    lineage_files = [x for x in lineage_files if '_sample_state_profile.csv' in x]
    
    min_dist = 100
    called_ct = ''
    for lineage_file in lineage_files:
        lineage_tab = pd.read_csv(parent_folder + "/" + lineage_file, index_col=0)
        lineage_tab = lineage_tab.loc[sub_sim_tab.index, :]
        sim_ss = sub_sim_tab.iloc[:, -1]
        sim_ss = sim_ss > 1
        sim_ss = sim_ss.astype("int")
        sim_ss = np.array(sim_ss)

        real_ss = lineage_tab.iloc[:, -1]
        real_ss = np.array(real_ss)

        dist = np.linalg.norm(real_ss - sim_ss)
        if dist < min_dist:
            min_dist = dist
            called_ct = "-".join(lineage_tab.columns)
    return called_ct 
#calculate cross correlation

if 'cross_corr' not in os.listdir("../output"):
    os.mkdir("../output/cross_corr")
datasets = ['dyn-BF', 'dyn-BFC', 'dyn-TF', 'dyn-LI', 'dyn-CY']

for dataset in datasets: 
    sim_folder_path =  "../output/OneSC_sim/" + dataset + "/"
    curated_dict = cluster_simulation(sim_folder_path)

    if dataset == 'dyn-CY':
        new_list = []
        for temp_ss in curated_dict.keys():
            new_list = new_list + curated_dict[temp_ss]
        curated_dict = dict()
        curated_dict[temp_ss] = list(np.unique(new_list))
        
    curated_average = curate_average(sim_folder_path, curated_dict)
    
    exp_dat = pd.read_csv("../BEELINE-data/inputs/Synthetic/" + dataset + "/" + dataset + "-5000-1/ExpressionData.csv", index_col=0)
    samp_tab = pd.read_csv("../BEELINE-data/inputs/Synthetic/" + dataset + "/" + dataset + "-5000-1/PseudoTime.csv", index_col=0)

    compiled_test_df = pd.DataFrame()
    for ss_string in curated_average.keys():
        sim_exp = curated_average[ss_string]
        if dataset == 'dyn-CY':
            sim_exp = sim_exp.iloc[:, 0:1300]
        for temp_pt in samp_tab.columns:
            print(temp_pt)
            sub_samp_tab = samp_tab.loc[samp_tab[temp_pt].notnull(), :]
            sub_samp_tab = sub_samp_tab.sort_values(by=[temp_pt])
            sub_exp_dat = exp_dat.loc[sim_exp.index, sub_samp_tab.index]
            
            train_ss = sub_exp_dat.iloc[:, -1]
            train_ss = np.array(train_ss > 1)
            train_ss = train_ss.astype(int)
            train_ss = [str(x) for x in train_ss]
            train_ss_string = "_".join(train_ss)

            test_df = calc_cross_correlation(sim_exp, sub_exp_dat, sub_samp_tab, pt_col = temp_pt, resolution = 0.01, top_X = 40)

            test_df['sim_ss'] = ss_string
            test_df['train_ss'] = train_ss_string
            compiled_test_df = pd.concat([compiled_test_df, test_df])
    compiled_test_df.to_csv("../output/cross_corr/" + dataset + "_cross_corr.csv")

# run once for LI extended for BoolODE 
dataset = 'dyn-LI-Extended'
sim_folder_path =  "../output/OneSC_sim/" + 'dyn-LI' + "/"
curated_dict = cluster_simulation(sim_folder_path)
curated_average = curate_average(sim_folder_path, curated_dict)

exp_dat = pd.read_csv("../BoolODE/Debug/dyn-linear-1_extended/dyn-linear-1_extended-5000-1/ExpressionData.csv", index_col=0)
samp_tab = pd.read_csv("../BoolODE/Debug/dyn-linear-1_extended/dyn-linear-1_extended-5000-1/PseudoTime.csv", index_col=0)

compiled_test_df = pd.DataFrame()
for ss_string in curated_average.keys():
    sim_exp = curated_average[ss_string]
    for temp_pt in samp_tab.columns:
        print(temp_pt)
        sub_samp_tab = samp_tab.loc[samp_tab[temp_pt].notnull(), :]
        sub_samp_tab = sub_samp_tab.sort_values(by=[temp_pt])
        sub_exp_dat = exp_dat.loc[sim_exp.index, sub_samp_tab.index]
        
        train_ss = sub_exp_dat.iloc[:, -1]
        train_ss = np.array(train_ss > 1)
        train_ss = train_ss.astype(int)
        train_ss = [str(x) for x in train_ss]
        train_ss_string = "_".join(train_ss)

        test_df = calc_cross_correlation(sim_exp, sub_exp_dat, sub_samp_tab, pt_col = temp_pt, resolution = 0.01, top_X = 40)

        test_df['sim_ss'] = ss_string
        test_df['train_ss'] = train_ss_string
        compiled_test_df = pd.concat([compiled_test_df, test_df])
compiled_test_df.to_csv("../output/cross_corr/" + dataset + "_cross_corr.csv")




