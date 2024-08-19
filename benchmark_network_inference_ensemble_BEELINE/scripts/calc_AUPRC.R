library(ggplot2)
library(stringr)
library(DescTools)

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

find_AUPRC <- function(gs_grn, inferred_grn) { 
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