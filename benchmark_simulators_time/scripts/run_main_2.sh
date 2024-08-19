#!/bin/bash
# run BoolODE 
trial_list=(../output/random_networks/*)
for i in "${trial_list[@]}"
do
	random_networks_list=("$i"/*)
	for j in "${random_networks_list[@]}"
	do
		echo "$j"/BoolODE_config.yaml
		python ../BoolODE/boolode.py --config "$j"/BoolODE_config.yaml &>> "$j"/BoolODE_output.txt
		python ../BoolODE/boolode.py --config "$j"/BoolODE_config_heaviside.yaml &>> "$j"/BoolODE_output_heaviside.txt
	done

done
