library(ggplot2)
library(dplyr)
library(cowplot)

##### get the grns across different methods ##### 

GRN_methods = list.dirs("../output/curated_networks/myeloid_progenitors_full/", full.names = FALSE, recursive = FALSE)
GRN_methods = c(GRN_methods, 'IQCELL')
GRN_methods = GRN_methods[GRN_methods != 'run_OneSC']

OneSC_grn = read.csv("../output/curated_networks/myeloid_progenitors_full/run_OneSC/curated_network_density.csv", row.names = 1)
sim_df = data.frame(row.names = GRN_methods,
                    'methods' = GRN_methods, 
                    'percent_similarity' = NA)

for(tmp_method in GRN_methods) {
  if(tmp_method == 'IQCELL') {
    tmp_grn = read.csv("../Beeline_benchmark/run_iQcell/myeloid_progenitors_full/iQcell_network.csv")
    tmp_grn$X = NULL
  } else if(tmp_method == 'run_OneSC') {
    next()
  } else { 
    tmp_grn = read.csv(file.path("../output/curated_networks/myeloid_progenitors_full/", tmp_method, 'curated_network_density.csv'), row.names = 1)
  }
  percent_agree = length(intersect(paste0(OneSC_grn$TF, "_", OneSC_grn$TG, "_", OneSC_grn$Type), paste0(tmp_grn$TF, "_", tmp_grn$TG, "_", tmp_grn$Type))) / length(paste0(OneSC_grn$TF, "_", OneSC_grn$TG, "_", OneSC_grn$Type))
  sim_df[tmp_method, 'percent_similarity'] = percent_agree
}

rownames(sim_df) = toupper(rownames(sim_df))
sim_df$methods = toupper(sim_df$methods)
p = ggplot(sim_df, aes(x = reorder(methods, -percent_similarity), y = percent_similarity)) + 
  geom_bar(stat = 'identity') + 
  theme_cowplot() + 
  ylim(c(0, 1)) + 
  xlab('Methods') + 
  ylab("Percent similarity to OneSC's inferred network")

##### see if there is a correlation between percent similarity and steady states similarity ######
ss_sim = read.csv("../output/sim_ss_comparison/plot_df.csv", row.names = 1)
ss_sim_sum = ss_sim %>% 
  group_by(methods) %>% 
  summarise(average_percent = mean(percent_agreement))
ss_sim_sum = data.frame(ss_sim_sum)
rownames(ss_sim_sum) = ss_sim_sum$methods
sim_df[rownames(ss_sim_sum), 'percent_agreement'] = ss_sim_sum$average_percent
sim_df = sim_df[is.na(sim_df$methods) == FALSE, ]

p = ggplot(sim_df, aes(x = percent_similarity, y = percent_agreement)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text(aes(label = methods), vjust = -1, hjust = 0.5) +
  theme_cowplot() + 
  ylim(c(0, 1)) + 
  xlim(c(0.05, 0.51)) +
  xlab('Percent overlap with OneSC inferred network') + 
  ylab("Average top simulation terminal states similarity") + 
  annotate(
    "text",
    x = 0.3,  # X position of the text
    y = 0.25,  # Y position of the text
    label = paste("Pearson r =", round(cor(sim_df$percent_similarity, sim_df$percent_agreement), 2)),  # Text label with correlation
    size = 5,  # Font size
    color = "blue"  # Text color
  )


ggsave(filename = '../output/compare_networks_similarity/test.png', plot = p, width = 8, height = 6)
       