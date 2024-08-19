library(stringr)

curate_network_orphan <- function(GRN_style = 'BEELINE', data_type = 'dyn-BF') { 
  if(GRN_style %in% c('run_OneSC')) {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/", GRN_style, data_type, 'OneSC_network.csv'))
    inferred_grn$X = NULL
    return(inferred_grn)
  } else if(GRN_style %in% c('GENIE3', 'GRISLI', 'GRNBOOST2', 'LEAP', 'PIDC', 'PPCOR', 'SCRIBE', 'SINGE')) { 
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/real_data/", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    if(nrow(inferred_grn) == 0) {
      return(inferred_grn)
    }
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
    PPCOR_data = read.csv(file.path("../Beeline_benchmark/outputs/real_data/", data_type, 'PPCOR/rankedEdges.csv'), sep = '\t')
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
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/real_data/", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
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

curate_network_density <- function(GRN_style = 'BEELINE', data_type = 'dyn-BF', OneCC_style = 'run_OneSC') { 
  OneCC_net = read.csv(file.path("../Beeline_benchmark/", OneCC_style, data_type, 'OneSC_network.csv'))
  OneCC_net$X = NULL
  edge_num = nrow(OneCC_net)
  
  if(GRN_style %in% c('run_OneSC', 'run_OneCC_pyEpoch', 'run_OneCC_pyEpoch_thresh', 'run_OneCC_regulators')) {
    inferred_grn = read.csv(file.path("../Beeline_benchmark/", GRN_style, data_type, 'OneSC_network.csv'))
    inferred_grn$X = NULL
    return(inferred_grn)
  } else if(GRN_style %in% c('GENIE3', 'GRISLI', 'GRNBOOST2', 'LEAP', 'PIDC', 'PPCOR', 'SCRIBE', 'SINGE')) { 
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/real_data/", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    if(nrow(inferred_grn) == 0) {
      return(inferred_grn)
    }
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    
    if(edge_num < nrow(inferred_grn)){
      inferred_grn = inferred_grn[1:edge_num, ]
    }
    inferred_grn$Type = NA
    PPCOR_data = read.csv(file.path("../Beeline_benchmark/outputs/real_data/", data_type, 'PPCOR/rankedEdges.csv'), sep = '\t')
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
    inferred_grn = read.csv(file.path("../Beeline_benchmark/outputs/real_data/", data_type, GRN_style, 'rankedEdges.csv'), sep = '\t')
    inferred_grn = inferred_grn[order(abs(inferred_grn$EdgeWeight), decreasing = TRUE), ]
    rownames(inferred_grn) = seq(1, nrow(inferred_grn))
    if(edge_num < nrow(inferred_grn)){
      inferred_grn = inferred_grn[1:edge_num, ]
    }
    inferred_grn$Type = NA
    
    inferred_grn[inferred_grn$EdgeWeight >= 0, 'Type'] = "+"
    inferred_grn[inferred_grn$EdgeWeight < 0, 'Type'] = "-"
    
    colnames(inferred_grn)[grep("Gene1", colnames(inferred_grn))] = 'TF'
    colnames(inferred_grn)[grep("Gene2", colnames(inferred_grn))] = 'TG'
    inferred_grn = inferred_grn[, c('TF', 'TG', 'Type')]
    return(inferred_grn)
  }
}

datatypes = list.files("../Beeline_benchmark/outputs/real_data/")

datatypes = c('myeloid_progenitors_full')
dir.create(file.path("../output/curated_networks"))

for(data_type in datatypes) { 
  grn_styles = list.files(file.path("../Beeline_benchmark/outputs/real_data/", data_type))
  grn_styles = c(grn_styles, 'run_OneSC')
  dir.create(file.path("../output/curated_networks", data_type))
  
  for(grn_style in grn_styles) { 

    folder_path = file.path("../output/curated_networks", data_type, grn_style)
    curated_network = curate_network_orphan(grn_style, data_type)
    if(nrow(curated_network) == 0) {
      next()
    }
    dir.create(folder_path)
    
    withr::with_dir(folder_path, {
      write.csv(curated_network, file = 'curated_network_orphan.csv')
    })
  }
  print(data_type)
}

for(data_type in datatypes) { 
  grn_styles = list.files(file.path("../Beeline_benchmark/outputs/real_data/", data_type))
  grn_styles = c(grn_styles, 'run_OneSC')
  dir.create(file.path("../output/curated_networks", data_type))
  
  for(grn_style in grn_styles) { 
    
    folder_path = file.path("../output/curated_networks", data_type, grn_style)
    curated_network = curate_network_density(grn_style, data_type)
    if(nrow(curated_network) == 0) {
      next()
    }
    dir.create(folder_path)
    
    withr::with_dir(folder_path, {
      write.csv(curated_network, file = 'curated_network_density.csv')
    })
  }
  print(data_type)
}





