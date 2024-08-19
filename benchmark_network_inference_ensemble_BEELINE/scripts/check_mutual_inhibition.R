library(ggplot2)
library(stringr)

option_list = list()
option_list[['run_OneCC_full']] = 'OneCC_all_genes'
option_list[['run_OneCC_regulators']] = 'OneCC_preNet_act_window'
option_list[['run_OneCC_pyEpoch']] = 'OneCC_preNet_pyEpoch_top_reg'
option_list[['run_OneCC_pyEpoch_thresh']] = 'OneCC_preNet_pyEpoch_thresh'

datatypes = list.files("../output/curated_networks")
network_type = 'F1'

big_df = data.frame()

calc_correct_mutual_inhibition <- function(gs_grn, inferred_grn) {
  gs_grn$edge = paste0(gs_grn$TF, "_", gs_grn$TG)
  gs_grn$edge_rev = paste0(gs_grn$TG, "_", gs_grn$TF)
  
  sub_grn = gs_grn[grepl('-', gs_grn$Type), ] # only get inhibition 
  sub_grn = sub_grn[sub_grn$edge %in% sub_grn$edge_rev, ]
  
  if(nrow(sub_grn) == 0) { 
    return(c(9000, 9000))  
  }
  
  inferred_grn_sub = inferred_grn[inferred_grn$Type == '-', ]
  if(nrow(inferred_grn_sub) == 0) {
    return(c(0, 0))
  }
  
  inferred_grn_sub$edge = paste0(inferred_grn_sub$TF, "_", inferred_grn_sub$TG)
  inferred_grn_sub$edge_rev = paste0(inferred_grn_sub$TG, "_", inferred_grn_sub$TF)
  inferred_grn_sub = inferred_grn_sub[inferred_grn_sub$edge %in% inferred_grn_sub$edge_rev, ]
  

  correct_percent = length(intersect(inferred_grn_sub$edge, sub_grn$edge)) / length(sub_grn$edge)
  wrong_percent = length(setdiff(inferred_grn_sub$edge, sub_grn$edge)) / length(sub_grn$edge)
  
  return(c(correct_percent, wrong_percent))
}
for(datatype in datatypes) { 
  methods = list.files(file.path("../output/curated_networks", datatype))
  temp_df = data.frame(row.names = methods)
  temp_df$methods = methods
  temp_df$datatypes = datatype
  temp_df$percent_correct = NA
  temp_df$percent_error = NA
  gs_grn = read.csv(file.path("../Beeline_benchmark/inputs/example/", datatype, 'refNetwork.csv'))
  colnames(gs_grn) = c('TF', 'TG', 'Type')
  
  for(method in methods) { 
    inferred_grn = read.csv(file.path("../output/curated_networks", datatype, method, paste0("curated_network_", network_type, ".csv")), row.names = 1)
    correct_outputs = calc_correct_mutual_inhibition(gs_grn, inferred_grn)
    temp_df[method, 'percent_correct'] = correct_outputs[1]
    temp_df[method, 'percent_error'] = correct_outputs[2]
  }
  big_df = rbind(big_df, temp_df)
}

dir.create("../output/mutual_inhibition")
big_df$thres_style = 'max F1'
big_df = big_df[big_df$percent_correct < 999, ]

for(run_OneCC in names(option_list)) {
  print(run_OneCC)
  big_df[big_df$methods == run_OneCC, 'methods'] = option_list[[run_OneCC]]
}

p <- ggplot(big_df, aes(x=reorder(methods,-percent_correct), y=percent_correct)) + 
  geom_boxplot()+
  theme_bw() + 
  ylab("mutual inhibition - percent correct") + 
  xlab("methods") +
  ggtitle('GRNs selected using optimal F1') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))

ggsave(filename = '../output/mutual_inhibition/mutual_inhibition_percent_correct.png', plot = p, width = 10, height = 5)

p <- ggplot(big_df, aes(x=reorder(methods,percent_error), y=percent_error)) + 
  geom_boxplot()+
  theme_bw() + 
  ylab("mutual inhibition - percent error") + 
  xlab("methods") +
  ggtitle('GRNs selected using optimal F1') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))

ggsave(filename = '../output/mutual_inhibition/mutual_inhibition_percent_error.png', plot = p, width = 10, height = 5)


