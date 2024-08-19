library(ggplot2)
library(cowplot)
library(ggrepel)

data_types = list.dirs("../output/plot_UMAPs", recursive = FALSE, full.names = FALSE)

for(data_type in data_types) {
  UMAP_coord = read.csv(file.path('../output/plot_UMAPs/', data_type, "UMAP_coord.csv"), row.names = 1)
  correlation_df = read.csv(file.path("../output/plot_UMAPs/", data_type, "correlation_mat.csv"), row.names = 1)
  UMAP_coord$cell_types = rownames(correlation_df)[apply(correlation_df, MARGIN = 2, FUN = which.max)]
  UMAP_coord[UMAP_coord$steady_states == 'other', 'cell_types'] = 'other'

  p1 = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color=cell_types)) + 
    geom_point() +
    scale_color_brewer(palette = "Set2") +
    theme_half_open() + 
    labs(color = "Cell Types")
  
  p2 = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color=sim_time)) + 
    geom_point() +
    theme_half_open() + 
    labs(color = 'Simulation Steps') + 
    scale_colour_viridis_c()

  ggsave(filename = file.path(file.path('../output/plot_UMAPs/', data_type, "UMAP_ct.png")), plot = p1, width = 7, height = 5)
  ggsave(filename = file.path(file.path('../output/plot_UMAPs/', data_type, "UMAP_sim_time.png")), plot = p2, width = 7, height = 5)
}
