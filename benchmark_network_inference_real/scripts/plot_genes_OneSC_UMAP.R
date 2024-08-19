library(stringr)
library(ggplot2)
library(reshape2)
data_types = list.dirs("../output/run_correlation/", full.names = FALSE, recursive = FALSE)

# plot the train UMAP 
for(data_type in data_types) {
  query_cluster_df = read.csv(paste0("../output/run_correlation/", data_type, "/sample_tab.csv"), row.names = 1)

  p = ggplot(query_cluster_df, aes(x=UMAP_1, y=UMAP_2, color = cluster_id)) + 
    geom_point() + 
    theme_bw()
  
  ggsave(filename = file.path('../output/run_correlation', data_type, "query_UMAP.png"), width = 8, height = 6)
}

for(data_type in data_types) { 
  print(data_type)
  file_list = list.files(paste0("../output/run_correlation/", data_type))
  mean_files = file_list[grep("mean_corr", file_list)]
  for(mean_file in mean_files) {
    
    ss_name = stringr::str_remove_all(mean_file, '_mean_corr.csv')

    avg_corr = read.csv(file.path(paste0("../output/run_correlation/", data_type), mean_file), row.names = 1)
    avg_corr$clusters = rownames(avg_corr)
    
    max_corr = read.csv(file.path(paste0("../output/run_correlation/", data_type), paste0(ss_name, '_max_corr.csv')), row.names = 1)
    max_corr$clusters = rownames(max_corr)
    
    min_corr = read.csv(file.path(paste0("../output/run_correlation/", data_type), paste0(ss_name, '_min_corr.csv')), row.names = 1)
    min_corr$clusters = rownames(min_corr)
    
    avg_long <- melt(avg_corr, id.vars = c("clusters"))
    avg_long$variable = str_remove_all(avg_long$variable, "X")
    
    plot_df = avg_long
    colnames(plot_df) = c("clusters", "sim_time", "mean_correlation")
    
    plot_df$sim_time = as.numeric(plot_df$sim_time)
    #plot_df = plot_df[plot_df$sim_time < 2000, ]
    p<-ggplot(plot_df, aes(x=sim_time, y=mean_correlation, group=clusters)) +
      geom_line(aes(color=clusters)) + 
      ylim(c(0, 1)) + theme_bw() + 
      ylab("average correlation similarity") + 
      xlab("simulation time")
    ggsave(filename = file.path(paste0("../output/run_correlation/", data_type), paste0(ss_name, '_averagePlot.png')), plot = p, width = 8, height = 6)
  }
}
