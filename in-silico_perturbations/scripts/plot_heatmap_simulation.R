library(ggplot2)
library(cowplot)
library(ggrepel)
library(data.table)
data_type = 'myeloid_progenitors_full'
data_types_list = list.dirs("../output/plot_UMAPs", recursive = FALSE, full.names = FALSE)
my_colors = RColorBrewer::brewer.pal(name = 'Set2', n = 8)
names(my_colors) = sort(c("CMP", "GMP", 'MEP', 'Monocytes', 'Erythrocytes', 'Granulocytes', 'MK', 'other'))

for(data_type in data_types_list) {
  genes_list = list.dirs(file.path("../output/plot_UMAPs/", data_type), recursive = FALSE, full.name = FALSE)
  for(temp_gene in genes_list) {
    UMAP_coord = read.csv(file.path('../output/plot_UMAPs/', data_type, temp_gene, "UMAP_coord.csv"), row.names = 1)
    correlation_df = fread(file.path("../output/plot_UMAPs/", data_type, temp_gene, "correlation_mat.csv"))
    rownames(correlation_df) = correlation_df$V1
    correlation_df$V1 = NULL
    UMAP_coord$cell_types = rownames(correlation_df)[apply(correlation_df, MARGIN = 2, FUN = which.max)]
    
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
    ggsave(filename = file.path(file.path("../output/plot_UMAPs", data_type, temp_gene, "corr_sim_heatmap.png")), plot = p, width = 8, height = 12)
    
  }
}

