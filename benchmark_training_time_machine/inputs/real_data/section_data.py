import pandas as pd 
import numpy as np 
import os

lineage_dict = dict()

t_cell_dict = dict()
t_cell_dict['lineage_0'] = ['ETP', 'DN2', 'DN3']
lineage_dict['t_cell_diff'] = t_cell_dict

EB_dict = dict()
EB_dict['lineage_0'] = ['Primed', 'PGC']
EB_dict['lineage_1'] = ['Primed', 'Ect', 'NeurEct']
EB_dict['lineage_2'] = ['Primed', 'PrimStr']
lineage_dict['day_4_EB'] = EB_dict

myeloid_dict = dict()
myeloid_dict['lineage_0'] = ['CMP', 'MEP', 'E']
myeloid_dict['lineage_1'] = ['CMP', 'GMP', 'G']
myeloid_dict['lineage_2'] = ['CMP', 'GMP', 'M']
lineage_dict['myeloid_progenitors'] = myeloid_dict

myeloid_full_dict = dict()
myeloid_full_dict['lineage_0'] = ['CMP', 'MEP', 'Erythrocytes']
myeloid_full_dict['lineage_1'] = ['CMP', 'GMP', 'Granulocytes']
myeloid_full_dict['lineage_2'] = ['CMP', 'GMP', 'Monocytes']
myeloid_full_dict['lineage_3'] = ['CMP', 'MK']
lineage_dict['myeloid_progenitors_full'] = myeloid_full_dict

folders_list = os.listdir("../../../../curate_real_training_data/output/")
folders_list = [x for x in folders_list if "." not in x]

def add_char_to_string(my_list, my_string):
    return [x + my_string for x in my_list]

for temp_folder in folders_list: 
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    exp_data = pd.read_csv("../../../../curate_real_training_data/output/" + temp_folder + "/train_exp.csv", index_col=0)
    samp_data = pd.read_csv("../../../../curate_real_training_data/output/" + temp_folder + "/samp_tab.csv", index_col=0)
    pt_data = pd.read_csv("../../../../curate_real_training_data/output/" + temp_folder + "/pseudotime.csv", index_col=0)
    pt_data = pt_data.loc[samp_data.index, :]

    if temp_folder == 'myeloid_progenitors':
        exp_data.columns = add_char_to_string([str(x) for x in list(exp_data.columns)], "_cell")
        samp_data.index = add_char_to_string([str(x) for x in list(samp_data.index)], "_cell")
        pt_data.index = add_char_to_string([str(x) for x in list(pt_data.index)], "_cell")
    elif temp_folder == 'myeloid_progenitors_full':
        exp_data.columns = add_char_to_string([str(x) for x in list(exp_data.columns)], "_cell")
        samp_data.index = add_char_to_string([str(x) for x in list(samp_data.index)], "_cell")
        pt_data.index = add_char_to_string([str(x) for x in list(pt_data.index)], "_cell")
    temp_lineage_dict = lineage_dict[temp_folder]
    if len(temp_lineage_dict) == 1: 
        new_pt_data = pd.DataFrame(data=None, index = samp_data.index, columns = ['PseudoTime'])
        new_pt_data['PseudoTime'] = pt_data['PseudoTime']
    else:
        new_col_names = ['PseudoTime' + str(x) for x in range(1, len(temp_lineage_dict) + 1)]
        new_pt_data = pd.DataFrame(data=None, index = samp_data.index, columns = new_col_names)
        for i in range(0, len(temp_lineage_dict)):
            temp_key = 'lineage_' + str(i)
            clusters = temp_lineage_dict[temp_key]
            temp_col = 'PseudoTime' + str(i + 1)

            new_pt_data.loc[samp_data['cell_types'].isin(clusters), temp_col] =  pt_data.loc[samp_data['cell_types'].isin(clusters), 'PseudoTime']
    new_pt_data.to_csv(temp_folder + "/pseudotime.csv", sep = ',')
    exp_data.to_csv(temp_folder + "/train_exp.csv", sep = ',')
    samp_data.to_csv(temp_folder + "/samp_tab.csv", sep = ',')