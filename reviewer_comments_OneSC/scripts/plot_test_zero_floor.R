library(ggplot2)
library(cowplot)
library(ggrepel)

data_type = 'myeloid_progenitors_full'
UMAP_coord = read.csv(file.path('../output/test_zero_floor_OneSC/myeloid_progenitors_full/wildtype_UMAPs/UMAP_coord.csv'), row.names = 1)
correlation_df = read.csv(file.path("../output/test_zero_floor_OneSC/myeloid_progenitors_full/wildtype_UMAPs/dist_mat.csv"), row.names = 1)
UMAP_coord$cell_types = rownames(correlation_df)[apply(correlation_df, MARGIN = 2, FUN = which.min)]
UMAP_coord[UMAP_coord$steady_states == 'other', 'cell_types'] = 'other'


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


ggsave(filename = file.path(file.path("../output/test_zero_floor_OneSC/myeloid_progenitors_full/wildtype_UMAPs/UMAP_ct.png")), plot = p1, width = 5, height = 4)
ggsave(filename = file.path(file.path("../output/test_zero_floor_OneSC/myeloid_progenitors_full/wildtype_UMAPs/UMAP_sim_time.png")), plot = p2, width = 5, height = 4)
