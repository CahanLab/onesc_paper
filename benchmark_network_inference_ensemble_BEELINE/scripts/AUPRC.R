library(stringr)
library(tidyr)
library(tibble)
library(dplyr)
library(DescTools)
library(ggplot2)
library(cowplot)

calc_precision <- function(gs_edges, inferred_edges) { 
  num_TP = length(intersect(gs_edges, inferred_edges))
  num_FP = length(setdiff(inferred_edges, gs_edges))
  return(num_TP / (num_TP + num_FP))
}

calc_recall <- function(gs_edges, inferred_edges) { 
  num_TP = length(intersect(gs_edges, inferred_edges))
  num_FN = length(setdiff(gs_edges, inferred_edges))
  return(num_TP / (num_TP + num_FN))
}

dir.create("../output/AUPRC/")

curate_network <- function(GRN_style = 'BEELINE', data_type = 'dyn-BF') { 
  gs_grn = read.csv(file.path("../Beeline_benchmark/inputs/example/", data_type, 'refNetwork.csv'))
  colnames(gs_grn) = c('TF', 'TG', 'Type')
  if(GRN_style %in% c('run_OneSC')) {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/", GRN_style, data_type, 'OneSC_network.csv'))
    inferred_grn$X = NULL
    inferred_grn$EdgeWeight = 1
    return(inferred_grn)
  } else if (GRN_style == 'run_iQcell') {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/", GRN_style, data_type, 'iQcell_network.csv'))
    inferred_grn$X = NULL
    inferred_grn$EdgeWeight = 1
    return(inferred_grn)
  } else if(GRN_style %in% c('GENIE3', 'GRISLI', 'GRNBOOST2', 'LEAP', 'PIDC', 'SCRIBE', 'SINGE')) { 
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/example/", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    
    PPCOR_data = read.csv(file.path("../Beeline_benchmark/outputs/example", data_type, 'PPCOR/rankedEdges.csv'), sep = '\t')
    rownames(PPCOR_data) = paste0(PPCOR_data$Gene1, "_", PPCOR_data$Gene2)
    inferred_grn = inferred_grn[!duplicated(inferred_grn), ]
    rownames(inferred_grn) = paste0(inferred_grn$Gene1, "_", inferred_grn$Gene2)
    
    PPCOR_data = PPCOR_data[rownames(inferred_grn), ]
    inferred_grn[PPCOR_data$EdgeWeight < 0, 'EdgeWeight'] = -1 * inferred_grn[PPCOR_data$EdgeWeight < 0, 'EdgeWeight']

    inferred_grn$Type = NA
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    return(inferred_grn)
  } else if(GRN_style %in% c('PPCOR', 'PCOR', 'pyEpoch')) { 
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    
    inferred_grn$Type = NA
    
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    return(inferred_grn)
  } else if(GRN_style == 'GRNVBEM') {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    
    raw_scores = list.files(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style), full.names = FALSE, recursive = FALSE) %>%
      grep("outFile", .) %>%
      list.files(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style), full.names = FALSE, recursive = FALSE)[.] %>%
      file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style, .) %>%
      lapply(., read.csv, sep='\t') %>% 
      bind_rows
    
    for(temp_index in rownames(inferred_grn)) {
      sub_raw = raw_scores[raw_scores$Parent == inferred_grn[temp_index, 'Gene1'] & raw_scores$Child == inferred_grn[temp_index, 'Gene2'], ]
      temp_sign = sign(sub_raw$Weight[which.min(sub_raw$Probability - inferred_grn[temp_index, 'EdgeWeight'])])
      inferred_grn[temp_index, 'EdgeWeight'] = inferred_grn[temp_index, 'EdgeWeight'] * temp_sign
    }
    inferred_grn$Type = NA
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    return(inferred_grn)
  } else if(GRN_style == 'SCODE') {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    
    # get the raw folder ,
    raw_scores = data.frame()
    raw_folder_list = list.dirs(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style), full.names = FALSE, recursive = FALSE)
    for(raw_folder in raw_folder_list) {
      raw_mat = read.csv(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style, raw_folder, 'meanA.txt'), sep = '\t', header = FALSE)
      exp_tab = read.csv(file.path("../Beeline_benchmark/inputs/example/", data_type, 'ExpressionData.csv'), row.names = 1)
      
      rownames(raw_mat) = rownames(exp_tab)
      colnames(raw_mat) = rownames(exp_tab)
      long_tab = raw_mat %>%
        rownames_to_column("Gene1") %>%
        gather(Gene2, EdgeWeight, -Gene1)
      
      raw_scores = rbind(raw_scores, long_tab)
    }
    
    for(temp_index in rownames(inferred_grn)) {
      sub_raw = raw_scores[raw_scores$Gene1 == inferred_grn[temp_index, 'Gene1'] & raw_scores$Gene2 == inferred_grn[temp_index, 'Gene2'], ]
      temp_sign = sign(sub_raw$EdgeWeight[which.min(sub_raw$EdgeWeight - inferred_grn[temp_index, 'EdgeWeight'])])
      inferred_grn[temp_index, 'EdgeWeight'] = inferred_grn[temp_index, 'EdgeWeight'] * temp_sign
    }

    inferred_grn$Type = NA
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    return(inferred_grn)
  } else if(GRN_style == 'SINCERITIES') {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    
    raw_scores = list.files(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style), full.names = FALSE, recursive = FALSE) %>%
      grep("outFile", .) %>%
      list.files(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style), full.names = FALSE, recursive = FALSE)[.] %>%
      file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style, .) %>%
      lapply(., read.csv, sep=',') %>% 
      bind_rows
    
    for(temp_index in rownames(inferred_grn)) {
      sub_raw = raw_scores[raw_scores$SourceGENES == inferred_grn[temp_index, 'Gene1'] & raw_scores$TargetGENES == inferred_grn[temp_index, 'Gene2'], ]
      temp_sign_bool = sub_raw$Edges[which.min(sub_raw$Interaction - inferred_grn[temp_index, 'EdgeWeight'])] == 'activation'
      if(temp_sign_bool == TRUE) {
        temp_sign = 1
      } else {
        temp_sign = -1
      }
      inferred_grn[temp_index, 'EdgeWeight'] = inferred_grn[temp_index, 'EdgeWeight'] * temp_sign
    }
    inferred_grn$Type = NA
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    return(inferred_grn)
  }
}

datatypes = list.files("../Beeline_benchmark/outputs/example/")

calc_AUPRC <- function(gs_grn, inferred_network, grn_style) {
  if(grn_style %in% c("run_OneSC", "run_iQcell")) {
    precision = calc_precision(paste0(gs_grn$TF, "_", gs_grn$TG, "_", gs_grn$Type), 
                               paste0(inferred_network$TF, "_", inferred_network$TG, "_", inferred_network$Type))
    recall = calc_recall(paste0(gs_grn$TF, "_", gs_grn$TG, "_", gs_grn$Type), 
                         paste0(inferred_network$TF, "_", inferred_network$TG, "_", inferred_network$Type))
    return(precision * recall)
  } else {
    threshold_list = seq(from = 0, to = max(abs(inferred_network$EdgeWeight)), by = max(abs(inferred_network$EdgeWeight)) / 100)
    
    PR_df = data.frame(threshold = threshold_list)
    PR_df$precision = NA
    PR_df$recall = NA
    
    for(threshold in threshold_list) {
      sub_inferred_grn = inferred_network[abs(inferred_network$EdgeWeight) >= threshold, ]
      PR_df[PR_df$threshold == threshold, 'precision'] = calc_precision(paste0(gs_grn$TF, "_", gs_grn$TG, "_", gs_grn$Type), 
                                                                        paste0(sub_inferred_grn$TF, "_", sub_inferred_grn$TG, "_", sub_inferred_grn$Type))
      PR_df[PR_df$threshold == threshold, 'recall'] = calc_recall(paste0(gs_grn$TF, "_", gs_grn$TG, "_", gs_grn$Type), 
                                                                  paste0(sub_inferred_grn$TF, "_", sub_inferred_grn$TG, "_", sub_inferred_grn$Type))
    }
    
    return(DescTools::AUC(PR_df$recall, PR_df$precision))
  }

}

# run this because last time it was interrupted for whatever reasons 
big_df = data.frame()
for(data_type in datatypes) { 
  grn_styles = list.files(file.path("../Beeline_benchmark/outputs/example/", data_type))
  grn_styles = c(grn_styles, 'run_OneSC', 'run_iQcell')
  dir.create(file.path("../output/AUPRC", data_type))
  
  gs_grn = read.csv(file.path("../Beeline_benchmark/inputs/example/", data_type, 'refNetwork.csv'))
  colnames(gs_grn) = c('TF', 'TG', 'Type')
  
  for(grn_style in grn_styles) { 
    folder_path = file.path("../output/AUPRC", data_type, grn_style)
    dir.create(folder_path, recursive = TRUE)
    
    inferred_network = curate_network(grn_style, data_type)
    withr::with_dir(folder_path, {
      write.csv(inferred_network, file = 'curated_network_orphan.csv')
    })
    AUPRC = calc_AUPRC(gs_grn, inferred_network, grn_style)
    temp_df = data.frame('data_type' = c(data_type), 
                         'grn_style' = c(grn_style), 
                         'AUPRC' = c(AUPRC))
    big_df = rbind(big_df, temp_df)
    
  }
  print(data_type)
}

big_df[big_df$grn_style == 'run_OneSC', 'grn_style'] = 'OneSC'
big_df[big_df$grn_style == 'run_iQcell', 'grn_style'] = 'iQcell'

# remove the cycle 
big_df = big_df[big_df$data_type != "dyn-CY", ]
p <- ggplot(big_df, aes(x=reorder(grn_style,-AUPRC, mean), y=AUPRC)) + 
  geom_boxplot()+
  theme_half_open() + 
  ylab("AUPRC") + 
  xlab("Methods") +
  ggtitle('GRN Methods') + 
  ylim(c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))
ggsave(filename = '../output/AUPRC/AUPRC.png', plot = p, width = 10, height = 4)





