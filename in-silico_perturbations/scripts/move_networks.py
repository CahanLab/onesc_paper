import pandas as pd
import os
import pickle

# get the OneCC network file 
# get the states 
output_path = '../output/files_needed'

# the whole purpose of this script is to move the OneCC inferred networks from the other place to here 
input_path = '../../benchmark_network_inference_real/output/curated_networks'
bio_systems = os.listdir(input_path)
bio_systems = [x for x in bio_systems if x != '.DS_Store']
for bio_system in bio_systems:
    OneCC_network = pd.read_csv(os.path.join(input_path, bio_system, 'run_OneSC/curated_network_density.csv'), index_col=0)
    if os.path.isdir(os.path.join(output_path, bio_system)) == False:
        os.makedirs(os.path.join(output_path, bio_system))
    OneCC_network.to_csv(os.path.join(output_path, bio_system, "curated_network.csv"))


# move the UMAP object 
input_path = '../../benchmark_network_inference_real/output/plot_UMAPs/myeloid_progenitors_full'
train_obj = pickle.load(open(os.path.join(input_path, 'train_obj.pickle'), 'rb'))
pickle.dump(train_obj, open(os.path.join(output_path + "/myeloid_progenitors_full/train_obj.pickle"), 'wb'))

# move the UMAP object 
input_path = '../../benchmark_network_inference_real/output/plot_UMAPs/myeloid_progenitors_full'
train_obj = pickle.load(open(os.path.join(input_path, 'train_obj.pickle'), 'rb'))
pickle.dump(train_obj, open(os.path.join(output_path + "/myeloid_progenitors_full/train_obj.pickle"), 'wb'))

# copy the extracted steady states 
import distutils.dir_util
from_dir = "../../benchmark_network_inference_real/output/extract_states/myeloid_progenitors_full"
to_dir = "../output/files_needed/myeloid_progenitors_full/extract_states"
distutils.dir_util.copy_tree(from_dir, to_dir)


input_path = '../../benchmark_network_inference_real/Beeline_benchmark/run_OneSC'
bio_systems = os.listdir(input_path)
bio_systems = [x for x in bio_systems if os.path.isdir(os.path.join(input_path, x)) == True]

for bio_system in bio_systems:
    state_dict = pickle.load(open(os.path.join(input_path, bio_system, "state_dict.pickle"), 'rb'))
    pickle.dump(state_dict, open(os.path.join(output_path, bio_system, "states_dict.pickle"), 'wb'))


