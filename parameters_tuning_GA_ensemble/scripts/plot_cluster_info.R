library(ggplot2)
library(stringr)
library(ggsignif)
library(cowplot)
library(dplyr)

calc_precision <- function(gs_edges, inferred_edges) { 
  num_TP = length(intersect(gs_edges, inferred_edges))
  num_FP = length(setdiff(inferred_edges, gs_edges))
  return(num_TP / (num_TP + num_FP))
}

calc_recall <- function(gs_edges, inferred_edges) { 
  num_TP = length(intersect(gs_edges, inferred_edges))
  num_FN = length(setdiff(gs_edges, inferred_edges))
  return(num_TP / (num_TP + num_FN))
}

data_types = list.dirs("../output/screen_cluster_assignment/", recursive = FALSE, full.name = FALSE)

big_df = data.frame()
for(data_type in data_types) {
  gs_grn = read.csv(file.path("../input/gs", data_type, 'refNetwork.csv'))
  parameter_sweep_list = list.dirs(file.path("../output/screen_cluster_assignment", data_type), recursive = FALSE, full.name = FALSE)
  
  for(temp_parameter in parameter_sweep_list) {
    trial_list = list.dirs(file.path("../output/screen_cluster_assignment", data_type, temp_parameter), recursive = FALSE, full.name = FALSE)
    trial_df = data.frame(row.names = trial_list)
    trial_df['precision'] = NA
    trial_df['recall'] = NA
    for(temp_trial in trial_list) {
      print(temp_trial)
      inferred_grn = read.csv(file.path("../output/screen_cluster_assignment", data_type, temp_parameter, temp_trial, paste0("OneSC_network.csv")))
      inferred_grn$X = NULL
      trial_df[temp_trial, 'precision'] = calc_precision(paste0(gs_grn$Gene1, "_" ,gs_grn$Gene2, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", gs_grn$Type))
      trial_df[temp_trial, 'recall'] = calc_recall(paste0(gs_grn$Gene1, "_" ,gs_grn$Gene2, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", gs_grn$Type))
    }
    trial_df$F1 = NA
    trial_df$F1 = (2*(trial_df$precision * trial_df$recall) / (trial_df$precision + trial_df$recall))
    trial_df[is.na(trial_df)] = 0
    mean_F1 = mean(trial_df$F1)  
    sd_F1 = sd(trial_df$F1)
    
    temp_df = data.frame(mean_F1 = c(mean_F1), 
                         sd_F1 = c(sd_F1), 
                         data_type = c(data_type))
    if(stringr::str_split(temp_parameter, "_")[[1]][2] == 'split') {
      temp_df$cluster_deviation = c(as.numeric(stringr::str_split(temp_parameter, "_")[[1]][4]))
    } else {
      temp_df$cluster_deviation = c(-as.numeric(stringr::str_split(temp_parameter, "_")[[1]][4]))
    }
    big_df = rbind(big_df, temp_df)
  }
}

for(data_type in data_types) {
  inferred_grn = read.csv(file.path("../../benchmark_network_inference_ensemble_BEELINE/output/curated_networks/", data_type, "/run_OneSC/curated_network_F1.csv"), row.names = 1)
  gs_grn = read.csv(file.path("../input/gs", data_type, 'refNetwork.csv'))
  
  precision = calc_precision(paste0(gs_grn$Gene1, "_" ,gs_grn$Gene2, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", inferred_grn$Type))
  recall = calc_recall(paste0(gs_grn$Gene1, "_" ,gs_grn$Gene2, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", inferred_grn$Type))
  
  F1 = (2*(precision * recall) / (precision + recall))
  temp_df = data.frame(mean_F1 = c(F1), 
                       sd_F1 = c(0), 
                       data_type = c(data_type), 
                       cluster_deviation = 0)
  big_df = rbind(big_df, temp_df)
}

p = ggplot(big_df, aes(x = cluster_deviation, y = mean_F1)) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_F1 - sd_F1, ymax = mean_F1 + sd_F1), width = .1) +
  facet_wrap(~data_type, ncol = 2) + 
  ylab("Mean F1") + 
  xlab("Number of Clusters Deviation") + 
  theme_half_open()

ggsave(filename = '../output/screen_cluster_assignment/cluster_deivation_performance.png', plot = p, width = 10, height = 7)
