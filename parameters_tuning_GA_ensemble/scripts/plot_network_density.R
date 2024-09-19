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

data_types = list.dirs("../input/gs", recursive = FALSE, full.name = FALSE)
parameter_sweep_list = list.dirs("../output/screen_density", recursive = FALSE, full.name = FALSE)

big_df = data.frame()
for(temp_parameter in parameter_sweep_list) {
    for(data_type in data_types) {
        gs_grn = read.csv(file.path("../input/gs", data_type, 'refNetwork.csv'))
        inferred_grn = read.csv(file.path("../output/screen_density", temp_parameter, data_type, paste0("OneSC_network.csv")))
        inferred_grn$X = NULL

        temp_df = data.frame(row.names = c(data_type))
        temp_df$data_type = c(data_type)
        temp_df$parameter = c(temp_parameter)
        temp_df$precision = NA
        temp_df$recall = NA

        temp_df[data_type, 'precision'] = calc_precision(paste0(gs_grn$Gene1, "_" ,gs_grn$Gene2, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", inferred_grn$Type))
        temp_df[data_type, 'recall'] = calc_recall(paste0(gs_grn$Gene1, "_" ,gs_grn$Gene2, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", inferred_grn$Type))

        big_df = rbind(big_df, temp_df)
    }
}

big_df$F1 = NA
big_df$F1 = (2*(big_df$precision * big_df$recall) / (big_df$precision + big_df$recall))

big_df %>% 
  group_by(parameter) %>%
  summarize(mean_F1 = mean(F1)) %>% 
  arrange(mean_F1)

p <- ggplot(big_df, aes(x=parameter, y=F1)) + 
  geom_boxplot()+
  theme_half_open() + 
  ylab("F1 Score") + 
  xlab("Ideal subnetwork density") +
  ylim(c(0, 1)) + 
  ggtitle('OneSC GRNs from multiple ideal subnetwork density') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))

ggsave(filename = '../output/screen_density/F1_comparison.png', plot = p, width = 10, height = 4)
