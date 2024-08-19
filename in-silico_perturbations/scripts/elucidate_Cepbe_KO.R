library(ggplot2)
library(cowplot)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

my_colors = RColorBrewer::brewer.pal(name = 'Set2', n = 8)
names(my_colors) = sort(c("CMP", "GMP", 'MEP', 'Monocytes', 'Erythrocytes', 'Granulocytes', 'MK', 'other'))

plot_df = data.frame()
# gather the wildtype similarity 
wt_UMAP = read.csv("../../benchmark_network_inference_real/output/plot_UMAPs/myeloid_progenitors_full/UMAP_coord.csv", row.names = 1)
wt_UMAP = wt_UMAP[wt_UMAP$steady_states == 'Granulocytes', ]
wt_UMAP$GMP_status = NA

sim_df = read.csv("../../benchmark_network_inference_real/output/plot_UMAPs/myeloid_progenitors_full/correlation_mat.csv", row.names = 1, header = FALSE)
colnames(sim_df) = sim_df[1, ]
sim_df = sim_df[2:nrow(sim_df), ]
sim_df = sim_df[, rownames(wt_UMAP)]

wt_UMAP$closest_ct = rownames(sim_df)[apply(sim_df, MARGIN = 2, FUN = which.max)]
wt_UMAP[wt_UMAP$closest_ct  == 'GMP', 'GMP_status'] = 1
wt_UMAP[wt_UMAP$closest_ct  != 'GMP', 'GMP_status'] = 0

summary_df = wt_UMAP %>% 
  group_by(., sim_time) %>% 
  summarize(avg = mean(GMP_status))

colnames(summary_df) = c('sim_time', 'proportion_GMP', 'condition')
summary_df$condition = 'wild type'
plot_df = rbind(plot_df, summary_df)

# gather the Cepba knockout 
ko_UMAP = read.csv("../output/plot_UMAPs/knockout/Cebpe/UMAP_coord.csv", row.names = 1)
ko_UMAP = ko_UMAP[ko_UMAP$sample %in% wt_UMAP$sample, ]

ko_sim_df = read.csv("../output/plot_UMAPs/knockout/Cebpe/correlation_mat.csv", row.names = 1, header = FALSE)
colnames(ko_sim_df) = ko_sim_df[1, ]
ko_sim_df = ko_sim_df[2:nrow(ko_sim_df), ]
ko_sim_df = ko_sim_df[, rownames(ko_UMAP)]

ko_UMAP$closest_ct = rownames(ko_sim_df)[apply(ko_sim_df, MARGIN = 2, FUN = which.max)]
ko_UMAP[ko_UMAP$closest_ct  == 'GMP', 'GMP_status'] = 1
ko_UMAP[ko_UMAP$closest_ct  != 'GMP', 'GMP_status'] = 0

summary_df = ko_UMAP %>% 
  group_by(., sim_time) %>% 
  summarize(avg = mean(GMP_status))

colnames(summary_df) = c('sim_time', 'proportion_GMP', 'condition')
summary_df$condition = 'Cebpe KO'
plot_df = rbind(plot_df, summary_df)


ggplot(plot_df, aes(x = sim_time, y = proportion_GMP, color = condition)) + 
  geom_line() + 
  theme_half_open() + 
  scale_color_brewer(palette = "Set1") +
  ggtitle("Cells that become granulocytes") + 
  xlim(c(0, 1500))

