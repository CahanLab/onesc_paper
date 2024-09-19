import os 
import pickle
import pandas as pd 

pickle_files_list = [x for x in os.listdir("../output") if 'pickle' in x]

machine_list = list()
runtime_list = list()

for pickle_file in pickle_files_list:
    my_time = pickle.load(open('../output/' + pickle_file, 'rb'))
    runtime_list.append(my_time)
    machine_list.append(pickle_file.split("_")[0])

runtime_df = pd.DataFrame()
runtime_df['machine_type'] = machine_list
runtime_df['runtime'] = runtime_list
runtime_df.loc[4, 'machine_type'] = 'Ryzen5600x.windows'
runtime_df.loc[4, 'runtime'] = 1106.9079155921936
runtime_df['cores'] = [8, 4, 10, 16, 12]
runtime_df['OS'] = ['Ubuntu', 'Ubuntu', 'MacOS', 'Ubuntu', 'Windows']

runtime_df.to_csv("../output/compiled_runtime.csv")