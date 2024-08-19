library(ggplot2)
library(cowplot)
library(ggrepel)

data_types = list.dirs("../output/plot_UMAPs", recursive = FALSE, full.names = FALSE)

data_type = 'myeloid_progenitors_full'
UMAP_coord = read.csv(file.path('../output/plot_UMAPs/', data_type, "UMAP_coord.csv"), row.names = 1)
correlation_df = read.csv(file.path("../output/plot_UMAPs/", data_type, "dist_mat.csv"), row.names = 1)
UMAP_coord$cell_types = rownames(correlation_df)[apply(correlation_df, MARGIN = 2, FUN = which.min)]
UMAP_coord[UMAP_coord$steady_states == 'other', 'cell_types'] = 'other'

for(time_step in unique(UMAP_coord$sim_time)) {

  sub_UMAP_coord = UMAP_coord[UMAP_coord$sim_time <= time_step, ]
  p = ggplot(sub_UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color=sim_time)) + 
    geom_point() +
    theme_half_open() + 
    labs(color = 'Simulation Steps') + 
    scale_colour_viridis_c() + 
    xlim(c(-5, 20)) + 
    ylim(c(-6, 18))
  ggsave(filename = file.path(file.path('../output/plot_UMAPs/', data_type, 'animation_folder', paste0(time_step, "_UMAP.png"))), plot = p, width = 5, height = 4)
  
}

p1 = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color=cell_types)) + 
  geom_point() +
  scale_color_brewer(palette = "Set2") +
  theme_half_open() + 
  labs(color = "Cell Types") + 
  ylab('UMAP_2') + 
  xlab('UMAP_1') +
  theme(text = element_text(size = 14), 
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        legend.position = 'right')   

p2 = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color=sim_time)) + 
  geom_point() +
  theme_half_open() + 
  labs(color = 'Simulation Steps') + 
  scale_colour_viridis_c() + 
  ylab('UMAP_2') + 
  xlab('UMAP_1') +
  theme(text = element_text(size = 14), 
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        legend.position = 'right')   


ggsave(filename = file.path(file.path('../output/plot_UMAPs/', data_type, "UMAP_ct.png")), plot = p1, width = 5, height = 4)
ggsave(filename = file.path(file.path('../output/plot_UMAPs/', data_type, "UMAP_sim_time.png")), plot = p2, width = 5, height = 4)

##### load in the expression profiles ######
exp_df = read.csv("../output/plot_UMAPs/myeloid_progenitors_full/big_sim_df.csv", header = FALSE, row.names = 1)
exp_df = t(exp_df)
exp_df = as.data.frame(exp_df)
rownames(exp_df) = exp_df$V1
exp_df$V1 = NULL
for(temp_gene in colnames(exp_df)) {
  plotting_UMAP_coord = UMAP_coord[, c('UMAP_1', 'UMAP_2', temp_gene)]
  colnames(plotting_UMAP_coord) = c('UMAP_1', 'UMAP_2', 'expression')
  plotting_UMAP_coord = plotting_UMAP_coord[order(plotting_UMAP_coord$expression), ]
  p_gene = ggplot(plotting_UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color=expression)) + 
    geom_point() +
    theme_half_open() + 
    labs(color = temp_gene) + 
    scale_colour_viridis_c(option = 'plasma')
  ggsave(filename = file.path(file.path('../output/plot_UMAPs/', data_type, paste0(temp_gene, "_UMAP.png"))), plot = p_gene, width = 4, height = 4)
}
