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

option_list = list()
option_list[['run_OneSC']] = 'OneSC'
option_list[['run_iQcell']] = 'iQcell'

datatypes = list.files("../output/curated_networks")
network_type = 'F1'
big_df = data.frame()
for(datatype in datatypes) { 
  methods = list.files(file.path("../output/curated_networks", datatype))
  temp_df = data.frame(row.names = methods)
  temp_df$methods = methods
  temp_df$datatypes = datatype
  temp_df$precision = NA
  temp_df$recall = NA
  gs_grn = read.csv(file.path("../Beeline_benchmark/inputs/example/", datatype, 'refNetwork.csv'))
  colnames(gs_grn) = c('TF', 'TG', 'Type')
  for(method in methods) { 
    inferred_grn = read.csv(file.path("../output/curated_networks", datatype, method, paste0("curated_network_", network_type, ".csv")), row.names = 1)
    temp_df[method, 'precision'] = calc_precision(paste0(gs_grn$TF, "_" ,gs_grn$TG, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", inferred_grn$Type))
    temp_df[method, 'recall'] = calc_recall(paste0(gs_grn$TF, "_" ,gs_grn$TG, "_", gs_grn$Type), paste0(inferred_grn$TF, "_", inferred_grn$TG, "_", inferred_grn$Type))
  }
  big_df = rbind(big_df, temp_df)
}

big_df$F1 = NA
big_df$F1 = (2*(big_df$precision * big_df$recall) / (big_df$precision + big_df$recall))

for(run_OneCC in names(option_list)) {
  print(run_OneCC)
  big_df[big_df$methods == run_OneCC, 'methods'] = option_list[[run_OneCC]]
}

dir.create("../output/Performance_F1")
big_df$thres_style = 'max F1'
big_df_F1 = big_df

write.csv(big_df_F1, '../output/Performance_F1/F1_F1_comparison.csv')
big_df_F1 = big_df_F1[is.na(big_df_F1$F1) == FALSE, ]
big_df_F1$methods = toupper(big_df_F1$methods)
p <- ggplot(big_df_F1, aes(x=reorder(methods,-F1, mean), y=F1)) + 
  geom_boxplot()+
  theme_half_open() + 
  ylab("F1 Score") + 
  xlab("Methods") +
  ggtitle('GRNs selected using optimal F1') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1), text = element_text(family = "Arial"))
ggsave(filename = '../output/Performance_F1/F1_F1_comparison.png', plot = p, width = 10, height = 4)

# t-test the difference between OneSC and GRNBOOST2 
t.test(big_df_F1[big_df_F1$methods == 'ONESC', 'F1'], big_df_F1[big_df_F1$methods == 'GRNBOOST2', 'F1'])
# get the average of F1, precision and recall
big_df_F1 %>% 
  group_by(methods) %>% 
  summarise_at(vars(F1), mean)

big_df_F1 %>% 
  group_by(methods) %>% 
  summarise_at(vars(precision), mean)

big_df_F1 %>% 
  group_by(methods) %>% 
  summarise_at(vars(recall), mean)

p <- ggplot(big_df_F1, aes(x=reorder(methods,-precision, mean), y=precision)) + 
  geom_boxplot()+
  theme_half_open() + 
  ylab("precision") + 
  xlab("Methods") +
  ggtitle('GRNs selected using optimal F1') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))

ggsave(filename = '../output/Performance_F1/F1_precision_comparison.png', plot = p, width = 10, height = 4)

p <- ggplot(big_df_F1, aes(x=reorder(methods,-recall), y=recall)) + 
  geom_boxplot()+
  theme_half_open() + 
  ylab("recall") + 
  xlab("Methods") +
  ggtitle('GRNs selected using optimal F1') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))

ggsave(filename = '../output/Performance_F1/F1_recall_comparison.png', plot = p, width = 10, height = 4)

big_df_F1 = big_df_F1[big_df_F1$methods != 'PCOR', ]

big_df_F1 = big_df_F1
cat_ss = big_df_F1 %>% 
  group_by(methods) %>%
  summarise(F1 = mean(F1, na.rm = TRUE)) %>% 
  arrange(F1)

big_df_F1$methods = factor(big_df_F1$methods, levels = cat_ss$methods)
write.csv(cat_ss, file = '../output/Performance_F1/category_ordering.csv')

big_df_F1$F1 = round(big_df_F1$F1, digits = 2)
p = ggplot(big_df_F1, aes(datatypes, methods, fill= F1)) + 
  geom_tile() +
  geom_text(aes(label = F1)) + 
  ylab("Methods") +
  xlab("Data Types") +
  scale_fill_viridis_c(option = 'A', alpha = 0.8) + 
  guides(fill=guide_legend(title="F1 Scores")) +
  theme_half_open() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None")
ggsave(filename = '../output/Performance_F1/F1_heatmaps.png', plot = p, width = 5, height = 6)

p = ggplot(big_df_F1, aes(datatypes, methods, fill= F1)) + 
  geom_tile() +
  geom_text(aes(label = F1)) + 
  ylab("Methods") +
  xlab("Data Types") +
  scale_fill_viridis_c(option = 'A', alpha = 0.8) + 
  guides(fill=guide_legend(title="F1 Scores")) +
  theme_half_open() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "left")
ggsave(filename = '../output/Performance_F1/F1_heatmaps_annotate.png', plot = p, width = 5, height = 6)
