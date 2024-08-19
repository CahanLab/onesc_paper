library(stringr)
library(tidyr)
library(tibble)
library(dplyr)
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

find_best_thresh <- function(gs_grn, inferred_grn) { 
  gs_grn$gold_standard = paste0(gs_grn$TF, "_", gs_grn$TG, "_", gs_grn$Type)
  inferred_grn$temp_type = NA
  inferred_grn[inferred_grn$EdgeWeight >= 0, 'temp_type'] = "+"
  inferred_grn[inferred_grn$EdgeWeight < 0, 'temp_type'] = "-"
  
  inferred_grn$Edges = paste0(inferred_grn$Gene1, "_", inferred_grn$Gene2, "_", inferred_grn$temp_type)
  inferred_grn$EdgeWeight = abs(inferred_grn$EdgeWeight)
  possible_threshes = unique(inferred_grn$EdgeWeight)
  gs_edges = gs_grn$gold_standard
  
  best_thresh = possible_threshes[1]
  best_F1 = 0 
  
  for(thresh in possible_threshes) { 
    inferred_edges = inferred_grn[inferred_grn$EdgeWeight >= thresh, 'Edges'] 
    print(length(inferred_edges))
    prec = calc_precision(gs_edges, inferred_edges)
    recall = calc_recall(gs_edges, inferred_edges)
    if(prec+recall == 0) { 
      new_F1 = 0 
    }
    else { 
      new_F1 = (2 * prec * recall) / (prec + recall)
    }
    if(new_F1 > best_F1) { 
      best_F1 = new_F1
      best_thresh = thresh
    }
  }
  return(best_thresh)
}

curate_network_F1 <- function(GRN_style = 'BEELINE', data_type = 'dyn-BF') { 
  gs_grn = read.csv(file.path("../Beeline_benchmark/inputs/example/", data_type, 'refNetwork.csv'))
  colnames(gs_grn) = c('TF', 'TG', 'Type')
  if(GRN_style %in% c('run_OneSC')) {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/", GRN_style, data_type, 'OneSC_network.csv'))
    inferred_grn$X = NULL
    return(inferred_grn)
  } else if (GRN_style == 'run_iQcell') {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/", GRN_style, data_type, 'iQcell_network.csv'))
    inferred_grn$X = NULL
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
    
    thresh = find_best_thresh(gs_grn, inferred_grn)
    inferred_grn = inferred_grn[abs(inferred_grn$EdgeWeight) >= thresh, ]
    
    inferred_grn$Type = NA
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    inferred_grn = inferred_grn[, c('TF', 'TG', 'Type')]
    return(inferred_grn)
  } else if(GRN_style %in% c('PPCOR', 'PCOR', 'pyEpoch')) { 
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    
    thresh = find_best_thresh(gs_grn, inferred_grn)
    inferred_grn = inferred_grn[abs(inferred_grn$EdgeWeight) >= thresh, ]
    
    inferred_grn$Type = NA
    
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    inferred_grn = inferred_grn[, c('TF', 'TG', 'Type')]
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
    thresh = find_best_thresh(gs_grn, inferred_grn)
    inferred_grn = inferred_grn[abs(inferred_grn$EdgeWeight) >= thresh, ]
    inferred_grn$Type = NA
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    inferred_grn = inferred_grn[, c('TF', 'TG', 'Type')]
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
    thresh = find_best_thresh(gs_grn, inferred_grn)
    inferred_grn = inferred_grn[abs(inferred_grn$EdgeWeight) >= thresh, ]
    inferred_grn$Type = NA
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    inferred_grn = inferred_grn[, c('TF', 'TG', 'Type')]
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
    thresh = find_best_thresh(gs_grn, inferred_grn)
    inferred_grn = inferred_grn[abs(inferred_grn$EdgeWeight) >= thresh, ]
    inferred_grn$Type = NA
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    inferred_grn = inferred_grn[, c('TF', 'TG', 'Type')]
    return(inferred_grn)
  }
}

curate_network_orphan <- function(GRN_style = 'BEELINE', data_type = 'dyn-BF') { 
  if(GRN_style == 'BEELINE') { 
    gs_grn = read.csv(file.path("../Beeline_benchmark/inputs/example/", data_type, 'refNetwork.csv'))
    colnames(gs_grn) = c('TF', 'TG', 'Type')
    return(gs_grn)  
  } else if(GRN_style %in% c('run_OneSC')) {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/", GRN_style, data_type, 'OneSC_network.csv'))
    inferred_grn$X = NULL
    return(inferred_grn)
  } else if (GRN_style == 'run_iQcell') {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/", GRN_style, data_type, 'iQcell_network.csv'))
    inferred_grn$X = NULL
    return(inferred_grn)
  } else if(GRN_style %in% c('GENIE3', 'GRISLI', 'GRNBOOST2', 'LEAP', 'PIDC', 'PPCOR', 'SCRIBE', 'SINGE', 'PCOR')) { 
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/example/", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    latest_index = which(inferred_grn$EdgeWeight == inferred_grn$EdgeWeight[1])[length(which(inferred_grn$EdgeWeight == inferred_grn$EdgeWeight[1]))]
    for(gene in unique(inferred_grn$Gene2)) { 
      if(which(inferred_grn$Gene2 == gene)[1] > latest_index) { 
        latest_index = which(inferred_grn$Gene2 == gene)[1]
      }
    }
    inferred_grn = inferred_grn[1:latest_index, ]
    inferred_grn$Type = NA
    PPCOR_data = read.csv(file.path("../Beeline_benchmark/outputs/example", data_type, 'PPCOR/rankedEdges.csv'), sep = '\t')
    rownames(PPCOR_data) = paste0(PPCOR_data$Gene1, "_", PPCOR_data$Gene2)
    inferred_grn = inferred_grn[!duplicated(inferred_grn), ]
    rownames(inferred_grn) = paste0(inferred_grn$Gene1, "_", inferred_grn$Gene2)
    
    PPCOR_data = PPCOR_data[rownames(inferred_grn), ]
    inferred_grn[PPCOR_data$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[PPCOR_data$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    inferred_grn = inferred_grn[, c('TF', 'TG', 'Type')]
    return(inferred_grn)
  } else if(GRN_style %in% c('GRNVBEM', 'SCODE', 'SINCERITIES', 'pyEpoch')) { 
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/example", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    latest_index = which(inferred_grn$EdgeWeight == inferred_grn$EdgeWeight[1])[length(which(inferred_grn$EdgeWeight == inferred_grn$EdgeWeight[1]))]
    for(gene in unique(inferred_grn$Gene2)) { 
      if(which(inferred_grn$Gene2 == gene)[1] > latest_index) { 
        latest_index = which(inferred_grn$Gene2 == gene)[1]
      }
    }
    inferred_grn = inferred_grn[1:latest_index, ]
    inferred_grn$Type = NA
    
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    inferred_grn = inferred_grn[, c('TF', 'TG', 'Type')]
    return(inferred_grn)
  }
}

datatypes = list.files("../Beeline_benchmark/outputs/example/")

# run this because last time it was interrupted for whatever reasons 
dir.create(file.path("../output/curated_networks"))
for(data_type in datatypes) { 
  grn_styles = list.files(file.path("../Beeline_benchmark/outputs/example/", data_type))
  grn_styles = c(grn_styles, 'run_OneSC', 'run_iQcell')
  dir.create(file.path("../output/curated_networks", data_type))
  
  for(grn_style in grn_styles) { 
    folder_path = file.path("../output/curated_networks", data_type, grn_style)
    dir.create(folder_path)
    
    curated_network = curate_network_orphan(grn_style, data_type)
    withr::with_dir(folder_path, {
      write.csv(curated_network, file = 'curated_network_orphan.csv')
    })
    
    curated_network = curate_network_F1(grn_style, data_type)
    
    withr::with_dir(folder_path, {
      write.csv(curated_network, file = 'curated_network_F1.csv')
    })
  }
  print(data_type)
}






