library(ggplot2)
library(cowplot)
library(ggrepel)
library(RColorBrewer)

data_types_list = list.dirs("../output/plot_UMAPs", recursive = FALSE, full.names = FALSE)
my_colors = RColorBrewer::brewer.pal(name = 'Set2', n = 8)
names(my_colors) = sort(c("CMP", "GMP", 'MEP', 'Monocytes', 'Erythrocytes', 'Granulocytes', 'MK', 'other'))

for(data_type in data_types_list) {
  genes_list = list.dirs(file.path("../output/plot_UMAPs/", data_type), recursive = FALSE, full.name = FALSE)
  for(temp_gene in genes_list) {
    UMAP_coord = read.csv(file.path("../output/plot_UMAPs", data_type, temp_gene, "UMAP_coord.csv"), row.names = 1)
    correlation_df = read.csv(file.path("../output/plot_UMAPs/", data_type, temp_gene, "correlation_mat.csv"), row.names = 1)
    UMAP_coord$cell_types = rownames(correlation_df)[apply(correlation_df, MARGIN = 2, FUN = which.max)]
    UMAP_coord[UMAP_coord$steady_states == 'other', 'cell_types'] = 'other'
    
    p1 = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color=cell_types)) + 
      geom_point() +
      scale_color_manual(values = my_colors) +
      theme_half_open() + 
      labs(color = "Cell Types") + 
      ggtitle(paste0(temp_gene, " ", data_type))
    
    p2 = ggplot(UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color=sim_time)) + 
      geom_point() +
      theme_half_open() + 
      labs(color = 'Simulation Steps') + 
      scale_colour_viridis_c() + 
      ggtitle(paste0(temp_gene, " ", data_type))
    
    ggsave(filename = file.path(file.path("../output/plot_UMAPs", data_type, temp_gene, "UMAP_ct.png")), plot = p1, width = 5, height = 3)
    ggsave(filename = file.path(file.path("../output/plot_UMAPs", data_type, temp_gene, "UMAP_sim_time.png")), plot = p2, width = 5, height = 3)
    
  }
}


##### load in the expression profiles ######
'''
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
'''