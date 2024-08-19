library(stringr)
library(BoolNet)

write_boolean <- function(sub_network) { 
  target_string = unique(sub_network$TG)
  function_string = ''
  pos_genes = as.vector(sub_network[sub_network$Type == "+", 'TF'])
  if(length(pos_genes) == 0) { 
    pos_genes = ''
  } else{ 
    pos_genes = toString(pos_genes)
    pos_genes = str_replace_all(pos_genes, ", ", " | ")
    pos_genes = paste0("(", pos_genes, ")")
  }
  
  neg_genes = as.vector(sub_network[sub_network$Type == "-", 'TF'])
  if(length(neg_genes) == 0) { 
    neg_genes = ''
  } else{ 
    neg_genes = toString(neg_genes)
    neg_genes = str_replace_all(neg_genes, ", ", " | ")
    neg_genes = paste0("(", neg_genes, ")")
  }
  
  if(pos_genes == '') { 
    function_string = paste0(" ! ", neg_genes)
  } else if (neg_genes == '') { 
    function_string = paste0(pos_genes)
  } else { 
    function_string = paste0(pos_genes, " &", " ! ", neg_genes)
  }
  function_string = paste0(target_string, ", ", function_string)
  return(function_string)
}

write_boolean_network <- function(curated_network, save_path) { 
  
  network_lines = c("targets, functions")
  for(TG in unique(curated_network$TG)) { 
    sub_network = curated_network[curated_network$TG == TG, ]
    network_lines = c(network_lines, write_boolean(sub_network))
  }
  
  fileConn<-file(save_path)
  writeLines(network_lines, fileConn)
  close(fileConn)
}

get_stochastic_ss <- function(net, mean_avg, num_runs = 10000){ 
  ss_df = matrix(data = NA, nrow = length(net$genes), ncol = num_runs)
  rownames(ss_df) = net$genes
  for(i in seq(1, num_runs)) { 
    initial_states = c()
    for(gene in net$genes) { 
      initial_states = c(initial_states, mean_avg[gene])
    }
    initial_states = as.numeric(initial_states >= 1)
    attr = getAttractors(net, type = 'asynchronous', startStates = list(initial_states))
    attr =attractorsToLaTeX(attr)
    attr = attr$`1`
    if(is.null(attr) == TRUE) {
      next
    }
    ss_df[, i] = attr[, 1]
  }
  return(ss_df)
}

dir.create(file.path("../output/golden_Boolean_sim"))

data_type_run = list.dirs("../output/curated_networks", recursive = FALSE, full = FALSE)

states_list = list()

# Dyn-LI
temp_states_list = list()
states = c(1, 0, 0, 0, 0, 0, 0)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6", "g7")
temp_states_list[['init']] = states
states = c(0, 0, 0, 0, 0, 0, 1)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6", "g7")
temp_states_list[['steady_0']] = states
states_list[['dyn-LI']] = temp_states_list

# Dyn-LL
temp_states_list = list()
states = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
names(states) = c("g1", "g10", "g11", "g12", "g13", "g14", "g15", "g16", "g17", "g18", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9")
temp_states_list[['init']] = states
states = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
names(states) = c("g1", "g10", "g11", "g12", "g13", "g14", "g15", "g16", "g17", "g18", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9")
temp_states_list[['steady_0']] = states
states_list[['dyn-LL']] = temp_states_list

# Dyn-BF
temp_states_list = list()
states = c(1, 0, 0, 0, 0, 0, 1)
names(states) = c("g1", "g2", "g3", "g4", "g6", "g7", "g8")
temp_states_list[['init']] = states
states = c(0, 0, 0, 1, 0, 1, 0)
names(states) = c("g1", "g2", "g3", "g4", "g6", "g7", "g8")
temp_states_list[['steady_0']] = states
states = c(0, 0, 0, 0, 1, 0, 1)
names(states) = c("g1", "g2", "g3", "g4", "g6", "g7", "g8")
temp_states_list[['steady_1']] = states
states_list[['dyn-BF']] = temp_states_list

# Dyn-BFC
temp_states_list = list()
states = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", 'g9', 'g10')
temp_states_list[['init']] = states
states = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", 'g9', 'g10')
temp_states_list[['steady_0']] = states
states_list[['dyn-BFC']] = temp_states_list

# Dyn-TF
temp_states_list = list()
states = c(1, 0, 0, 0, 0, 0, 0, 0)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8")
temp_states_list[['init']] = states
states = c(0, 0, 0, 0, 0, 1, 1, 1)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8")
temp_states_list[['steady_0']] = states
states = c(0, 0, 0, 0, 1, 0, 0, 1)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8")
temp_states_list[['steady_1']] = states
states = c(0, 0, 0, 1, 0, 0, 0, 1)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8")
temp_states_list[['steady_2']] = states
states_list[['dyn-TF']] = temp_states_list

# Dyn-CY
temp_states_list = list()
states = c(1, 1, 1, 0, 0, 0)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6")
temp_states_list[['init']] = states
states = c(-1, -1, -1, -1, -1, -1)
names(states) = c("g1", "g2", "g3", "g4", "g5", "g6")
temp_states_list[['steady_0']] = states
states_list[['dyn-CY']] = temp_states_list

# HSC 
temp_states_list = list()
states = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
names(states) = c('Cebpa','Gata2', 'Pu1', 'Gata1', 'Fog1', 'Eklf', 'Fli1', 'Scl', 'cJun', 'EgrNab', 'Gfi1')
temp_states_list[['init']] = states
states_list[['HSC']] = temp_states_list

# GSD 
temp_states_list = list()
states = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
names(states) = c('UGR', 'CBX2', 'GATA4', 'WT1mKTS', 'WT1pKTS', 'NR5A1', 'NR0B1', 'SRY', 'SOX9', 'FGF9', 'PGD2', 'DMRT1', 'DHH', 'DKK1', 'AMH', 'WNT4', 'RSPO1', 'FOXL2', 'CTNNB1')
temp_states_list[['init']] = states
states_list[['GSD']] = temp_states_list

# VSC 
temp_states_list = list()
states = c(0, 0, 0, 0, 0, 0, 0, 0)
names(states) = c('Nkx61', 'Nkx62', 'Nkx22', 'Pax6', 'Dbx1', 'Dbx2', 'Olig2', 'Irx3')
temp_states_list[['init']] = states
states_list[['VSC']] = temp_states_list

# mCAD 
temp_states_list = list()
states = c(1, 0, 0, 0, 0)
names(states) = c('Sp8', 'Fgf8', 'Emx2', 'Pax6', 'Coup')
temp_states_list[['init']] = states
states_list[['mCAD']] = temp_states_list

#data_type_run = data_type_run[data_type_run != 'dyn-CY']
for(data_type in data_type_run) {
  print(data_type)
  method = 'gold_standard'
  curated_network = read.csv(file.path("../Beeline_benchmark/inputs/example/", data_type, 'refNetwork.csv'))
  colnames(curated_network) = c('TF', 'TG', 'Type')
  out_path = file.path("../output/golden_Boolean_sim", data_type)
  dir.create(out_path, recursive = TRUE)
  
  write_boolean_network(curated_network, file.path(out_path, 'boolean_net_F1.txt'))
  net <- loadNetwork(file.path(out_path, 'boolean_net_F1.txt'))
  ss_df = get_stochastic_ss(net, states_list[[data_type]][['init']], num_runs = 50000) 
  
  if(sum(!duplicated(t(ss_df))) == 1) { 
    single_vector = ss_df[, !duplicated(t(ss_df))]
    ss_df = matrix(data = NA, ncol = 1, nrow = length(single_vector))
    rownames(ss_df) = names(single_vector)
    ss_df[, 1] = single_vector
  } else { 
    ss_df = ss_df[, !duplicated(t(ss_df))]
  }
  write.csv(ss_df, file = file.path(out_path, 'unique_steady_states_profiles.csv'))
}

###### check the statistics of gold standard networks 
for(data_type in data_type_run) {
  print(data_type)
  method = 'gold_standard'
  curated_network = read.csv(file.path("../Beeline_benchmark/inputs/example/", data_type, 'refNetwork.csv'))
  colnames(curated_network) = c('TF', 'TG', 'Type')
  print(paste0('number of edges ', nrow(curated_network)))
  print(paste0('number of nodes ', length(unique(curated_network$TF))))
}
