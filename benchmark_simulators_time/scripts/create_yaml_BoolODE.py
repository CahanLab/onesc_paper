import yaml
import os
import io

with open('../BoolODE/config-files/example-config.yaml') as f:
    skeleton_dict = yaml.safe_load(f)

skeleton_dict['global_settings']['do_post_processing'] = False
skeleton_dict['jobs'][0]['name'] = 'BoolODE_simulation'
skeleton_dict['jobs'][0]['model_definition'] = 'BoolODE_network.txt'
skeleton_dict['jobs'][0]['num_cells'] = 5

folder_list = os.listdir("../output/random_networks")

for temp_folder in folder_list: 
    network_list = os.listdir("../output/random_networks/" + temp_folder)
    print(temp_folder)
    for temp_network in network_list: 
        print(temp_network)
        active_dict = skeleton_dict.copy()
        active_dict['global_settings']['model_dir'] = "../output/random_networks/" + temp_folder + "/" + temp_network + '/BoolODE_original'
        active_dict['global_settings']['output_dir'] = "../output/random_networks/" + temp_folder + "/" + temp_network + "/BoolODE_original"
        with io.open("../output/random_networks/" + temp_folder + "/" + temp_network + '/BoolODE_config.yaml', 'w', encoding='utf8') as outfile:
            yaml.dump(active_dict, outfile, default_flow_style=False, allow_unicode=True)

# this is to load in the heaviside heaviside
skeleton_dict['global_settings']['do_post_processing'] = False
skeleton_dict['global_settings']['modeltype'] = 'heaviside'
skeleton_dict['jobs'][0]['name'] = 'BoolODE_simulation'
skeleton_dict['jobs'][0]['model_definition'] = 'BoolODE_network.txt'
skeleton_dict['jobs'][0]['num_cells'] = 5

folder_list = os.listdir("../output/random_networks")

for temp_folder in folder_list: 
    network_list = os.listdir("../output/random_networks/" + temp_folder)
    for temp_network in network_list: 
        active_dict = skeleton_dict.copy()
        active_dict['global_settings']['model_dir'] = "../output/random_networks/" + temp_folder + "/" + temp_network + '/BoolODE_heaviside'
        active_dict['global_settings']['output_dir'] = "../output/random_networks/" + temp_folder + "/" + temp_network + "/BoolODE_heaviside"
        with io.open("../output/random_networks/" + temp_folder + "/" + temp_network + '/BoolODE_config_heaviside.yaml', 'w', encoding='utf8') as outfile:
            yaml.dump(active_dict, outfile, default_flow_style=False, allow_unicode=True)