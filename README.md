# OneSC (manuscript analysis code)
The analysis code for the OneSC manuscript

All benchmarking and simulations were done on AWS EC2 instance c5.4xlarge Ubuntu. The Python packages used are in requirements.txt. 

The analysis is seperated into folders corresponding the type of analysis. 

|Folder | Analysis Decription| Figures|
| --------- | --------------- | --------------- |
|[benchmark_network_inference_BEELINE](benchmark_network_inference_ensemble_BEELINE/)| benchmarking GRN inference against other GRN methods | Figure 2, S1-4|
|[parameters_tuning_GA](parameters_tuning_GA_ensemble)| testing various parameters of GA to see how it affects GRN inference performance| Figure S5-6|
|[benchmark_simulators_time](benchmark_simulators_time)| benchmarking simulator in terms of runtime (single core)| Figure 3A|
|benchmark_simulators_parallel_time| benchmarking simulator in terms of runtime (5 cores) | Figure 3B|
|benchmark_simulators_similarity| comparing simulator with BoolODE in terms of synthetic cell similarity| Figure 3C, S7-8| 
|curate_real_training_data| identify the dynamically expressed TFs. Process the raw scRNA-seq data|Figure S9, S11|
|benchmark_network_inference_real| apply OneSC to real mouse myeloid single-cell dataset| Figure 4B-E, S10, S12|
|in-silico_perturbations| apply OneSC to predict in-silico perturbations| Figure 5, S13|



