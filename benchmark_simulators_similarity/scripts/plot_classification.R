library(stringr)
library(ggplot2)

data_types = list.dirs("../output/run_pySCN/", full.names = FALSE, recursive = FALSE)

# plot the train UMAP 
for(data_type in data_types) {
  train_UMAP_coord = read.csv(paste0("../output/overlap_UMAP/", data_type, "_train_UMAP.csv"), row.names = 1)
  train_cluster_df = read.csv(paste0("../output/BEELINE_clusters/", data_type,"/compiled_sampTab_pt.csv"), row.names = 1)
  train_UMAP_coord$cluster_id = train_cluster_df$cluster_id
  p = ggplot(train_UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = cluster_id)) + 
    geom_point() + 
    theme_bw()
  ggsave(filename = file.path('../output/run_pySCN', data_type, "train_UMAP.png"), width = 8, height = 6)
}

for(data_type in data_types) {
  query_UMAP_coord = read.csv(paste0("../output/overlap_UMAP/", data_type, "_query_UMAP.csv"), row.names = 1)
  query_cluster_df = read.csv(paste0("../output/run_pySCN/", data_type, "/samp_tab.csv"), row.names = 1)
  query_cluster_df = query_cluster_df[rownames(query_UMAP_coord), ]
  
  query_UMAP_coord$cluster_id = query_cluster_df$SCN_class
  p = ggplot(query_UMAP_coord, aes(x=UMAP_1, y=UMAP_2, color = cluster_id)) + 
    geom_point() + 
    theme_bw()
  ggsave(filename = file.path('../output/run_pySCN', data_type, "query_UMAP.png"), width = 8, height = 6)
}
