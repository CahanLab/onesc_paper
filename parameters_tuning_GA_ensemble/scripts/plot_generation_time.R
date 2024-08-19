library(ggplot2)
library(stringr)
library(ggsignif)
library(cowplot)
library(dplyr)
library(viridis)

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

data_types = list.dirs("../input/gs", recursive = FALSE, full.name = FALSE)
parameter_sweep_list = list.dirs("../output/screen_generations/", recursive = FALSE, full.name = FALSE)

big_df = data.frame()
for(temp_parameter in parameter_sweep_list) {
    temp_param_df = data.frame()
    for(data_type in data_types) {
        gs_grn = read.csv(file.path("../input/gs", data_type, 'refNetwork.csv'))
        inferred_grn = read.csv(file.path("../output/screen_generations", temp_parameter, data_type, paste0("OneSC_network.csv")))
        inferred_grn$X = NULL
        string_split_list = stringr::str_split(temp_parameter, "_")[[1]]
        temp_df = data.frame(row.names = c(data_type))
        temp_df$data_type = c(data_type)
        temp_df$num_generation = c(string_split_list[2])
        temp_df$max_iter = c(string_split_list[4])
        temp_df$precision = NA
        temp_df$recall = NA

        temp_df[data_type, 'precision'] = calc_precision(paste0(gs_grn$Gene1, "_" ,gs_grn$Gene2, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", inferred_grn$Type))
        temp_df[data_type, 'recall'] = calc_recall(paste0(gs_grn$Gene1, "_" ,gs_grn$Gene2, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", inferred_grn$Type))
        temp_param_df = rbind(temp_param_df, temp_df)
    }
    temp_param_df$F1 = NA
    temp_param_df$F1 = (2*(temp_param_df$precision * temp_param_df$recall) / (temp_param_df$precision + temp_param_df$recall))
    
    temp_param_df = temp_param_df[temp_param_df$data_type != 'dyn-CY', ]
    additional_df = data.frame(row.names = c(temp_parameter))
    additional_df$num_generation = c(string_split_list[2])
    additional_df$max_iter = c(string_split_list[4])
    additional_df$mean_F1 = mean(temp_param_df$F1)
    additional_df$sd_F1 = sd(temp_param_df$F1)
    big_df = rbind(big_df, additional_df)
    
}

big_df = big_df[big_df$num_generation != 50, ]
big_df$max_iter = factor(x = big_df$max_iter, levels = c(5, 10, 15, 20, 25, 30))
p = ggplot(big_df, aes(num_generation, max_iter, fill= mean_F1)) + 
  geom_tile() +
  ylab("Max Number of Iterations") +
  xlab("Number of Generations") +
  scale_fill_viridis(discrete=FALSE) + 
  guides(fill=guide_legend(title="Mean F1")) +
  theme_half_open() 

ggsave(filename = '../output/screen_generations/F1_generations.png', plot = p, width = 8, height = 8)

sub_df = big_df[big_df$max_iter == 10, ]
sub_df$num_generation = as.numeric(sub_df$num_generation)
sub_df$mean_F1 = as.numeric(sub_df$mean_F1)
p = ggplot(sub_df, aes(x = num_generation, y = mean_F1)) + 
  geom_line() +
  geom_point() +
  ylab("Mean F1") +
  xlab("Number of Generations") + 
  ylim(c(0, 1)) +
  theme_half_open() + 
  ggtitle("Mean F1 across different number of generations parameter")
ggsave(filename = '../output/screen_generations/F1_generations_scatter.png', plot = p, width = 9, height = 3)

time_big_df = read.csv("../output/screen_generations/generation_time_big_df.csv", row.names = 1)
time_big_df = time_big_df[time_big_df$num_gen != 50, ]
time_big_df$max_iter = as.character(time_big_df$max_iter)
time_big_df$num_gen = as.character(time_big_df$num_gen)
time_big_df$max_iter = factor(x = time_big_df$max_iter, levels = c(5, 10, 15, 20, 25, 30))

p = ggplot(time_big_df, aes(num_gen, max_iter, fill= mean_time)) + 
  geom_tile() +
  ylab("Max Number of Iterations") +
  xlab("Number of Generations") +
  scale_fill_viridis(discrete=FALSE) + 
  guides(fill=guide_legend(title="Mean RunTime(Sec)")) +
  theme_half_open()

ggsave(filename = '../output/screen_generations/run_time_generations.png', plot = p, width = 8, height = 6)

sub_df = time_big_df[time_big_df$max_iter == 10, ]
sub_df$num_gen = as.numeric(sub_df$num_gen)
sub_df$mean_time = as.numeric(sub_df$mean_time)
p = ggplot(sub_df, aes(x = num_gen, y = mean_time)) + 
  geom_line() +
  geom_point() +
  ylab("Mean Run Time (Sec)") +
  xlab("Number of Generations") + 
  theme_half_open() + 
  ggtitle("Mean run time across different number of generations parameter")
ggsave(filename = '../output/screen_generations/run_time_generation_scatter.png', plot = p, width = 9, height = 3)
