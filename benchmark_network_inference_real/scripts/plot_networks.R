library(igraph)

grn_df = read.csv("../output/curated_networks/myeloid_progenitors_full/run_OneSC/curated_network_orphan.csv", row.names = 1)

Group = grn_df$Type

# Build network
net <- graph_from_data_frame(d=grn_df, directed=TRUE)
G = net
plot(G, vertex.shape= "rectangle", edge.arrow.size=.3, 
     edge.color=ifelse(E(G)$Type =="+", "#DE4D7E", "#4AA0E7"), 
     vertex.color="#f7f7b7", 
     vertex.size2=10, 
     vertex.size = 20,
     vertex.frame.color="#f7f7b7", 
     vertex.label.color="black", 
     vertex.label.cex=1.4, vertex.label.dist=0, edge.curved=0.2, 
     edge.width = 2, 
     layout=layout_nicely)
legend("topleft",bty = "n",
      legend=unique(Group),
      fill=c("#4AA0E7", "#DE4D7E"), border=NA)

