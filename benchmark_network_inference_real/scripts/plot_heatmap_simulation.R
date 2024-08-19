library(ggplot2)
library(cowplot)
library(ggrepel)
library(data.table)
data_type = 'myeloid_progenitors_full'

##### this is to just plot the correlation similarity #####
UMAP_coord = read.csv(file.path('../output/plot_UMAPs/', data_type, "UMAP_coord.csv"), row.names = 1)
correlation_df = fread(file.path("../output/plot_UMAPs/", data_type, "correlation_mat.csv"))
rownames(correlation_df) = correlation_df$V1
correlation_df$V1 = NULL
UMAP_coord$cell_types = rownames(correlation_df)[apply(correlation_df, MARGIN = 2, FUN = which.max)]

#UMAP_coord[UMAP_coord$steady_states == 'other', 'cell_types'] = 'other'

#UMAP_coord = UMAP_coord[UMAP_coord$steady_states != 'other', ]

p = ggplot(UMAP_coord, aes(sim_time, sample, fill= cell_types)) + 
  geom_tile() + 
  scale_fill_brewer(palette = "Set2") + 
  theme_half_open() + 
  facet_grid(
    rows = vars(steady_states),
    scales = "free_y",
    space = "free_y",
    switch = "x"
  ) + 
  theme(
    panel.spacing = unit(x = 0.2, units = "lines"),
    strip.background = element_blank(), 
    strip.text.y = element_text(angle = 0),
    #axis.text.x=element_blank(), #remove x axis labels
    #axis.ticks.x=element_blank(), #remove x axis ticks
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  ) + 
  ylab("Simulation Runs") + 
  xlab("Simulation Steps")
ggsave(filename = file.path(file.path('../output/plot_UMAPs/', data_type, "corr_sim_heatmap.png")), plot = p, width = 8, height = 12)

##### this is to just plot the distance #####
UMAP_coord = read.csv(file.path('../output/plot_UMAPs/', data_type, "UMAP_coord.csv"), row.names = 1)
correlation_df = fread(file.path("../output/plot_UMAPs/", data_type, "dist_mat.csv"))
rownames(correlation_df) = correlation_df$V1
correlation_df$V1 = NULL
UMAP_coord$cell_types = rownames(correlation_df)[apply(correlation_df, MARGIN = 2, FUN = which.min)]

#UMAP_coord[UMAP_coord$steady_states == 'other', 'cell_types'] = 'other'

#UMAP_coord = UMAP_coord[UMAP_coord$steady_states != 'other', ]

p = ggplot(UMAP_coord, aes(sim_time, sample, fill= cell_types)) + 
  geom_tile() + 
  scale_fill_brewer(palette = "Set2") + 
  theme_half_open() + 
  facet_grid(
    rows = vars(steady_states),
    scales = "free_y",
    space = "free_y",
    switch = "x"
  ) + 
  theme(
    panel.spacing = unit(x = 0.2, units = "lines"),
    strip.background = element_blank(), 
    strip.text.y = element_text(angle = 0),
    #strip.text.y = element_blank(),
    #axis.text.x=element_blank(), #remove x axis labels
    #axis.ticks.x=element_blank(), #remove x axis ticks
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  ) +
  ylab("Simulation Runs") + 
  xlab("Simulation Steps")
ggsave(filename = file.path(file.path('../output/plot_UMAPs/', data_type, "dist_sim_heatmap.png")), plot = p, width = 8, height = 12)

