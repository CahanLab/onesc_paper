library(ggplot2)
library(stringr)

data_types = list.files("../output/steady_states/")

TARGET_dir = '../output/steady_states_plot'
dir.create(TARGET_dir)

compile_gold_ss <- function(data_type = 'HSC') {
  lineage_files = list.files(file.path("../output/extract_states/", data_type))
  lineage_files = lineage_files[grepl("states.csv", lineage_files)]
  ss_list = list()
  for(lineage_file in lineage_files) {
    lineage_df = read.csv(file.path("../output/extract_states/", data_type, lineage_file), row.names = 1)
    ss_list[[colnames(lineage_df)[ncol(lineage_df)]]] = lineage_df[, colnames(lineage_df)[ncol(lineage_df)]]
  }
  
  gold_ss = data.frame(Reduce(rbind, ss_list))
  if(data_type != 't_cell_diff') {
    gold_ss = t(gold_ss)
  }
  colnames(gold_ss) = names(ss_list)
  rownames(gold_ss) = rownames(lineage_df)
  return(gold_ss)
}

# number of correct steady states 
big_num_df = data.frame(data_type = NULL, 
                        method = NULL, 
                        style = NULL, 
                        num_ss = NULL,
                        gold_num_ss = NULL)
for(data_type in data_types) {
  print(data_type)
  methods_list = list.files(file.path("../output/steady_states", data_type))
  for(method in methods_list) {
    style_list = list.files(file.path("../output/steady_states", data_type, method))
    for(style in style_list) {
      compiled_ss = read.csv(file.path("../output/steady_states", data_type, method, style, 'unique_steady_states.csv'), row.names = 1)
      gold_ss = compile_gold_ss(data_type)
      temp_num_df = data.frame(data_type = data_type, 
                              method = method, 
                              style = style, 
                              num_ss = ncol(compiled_ss),
                              gold_num_ss = ncol(gold_ss))
      big_num_df = rbind(big_num_df, temp_num_df)
    }
  }
}
big_num_df$ss_diff = big_num_df$num_ss - big_num_df$gold_num_ss
big_num_df[grepl('OneCC', big_num_df$method), 'method'] = 'OneCC_window'

for(style in style_list) {
  if(style == 'orphan') {
    temp_big_df = big_num_df[big_num_df$style == style | big_num_df$method == 'OneCC_pyEpoch', ]
    
  } else {
    temp_big_df = big_num_df[big_num_df$style == style, ]
  }
  p = ggplot(data=big_num_df, aes(x=reorder(method, ss_diff), y=ss_diff, fill=data_type)) +
    geom_bar(stat="identity", position=position_dodge()) + 
    theme_bw() + 
    ylab("Num of steady states diff") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(filename = file.path(TARGET_dir, style, 'num_ss_diff.png'), plot = p, width = 8, height= 6)

}

write.csv(big_num_df, file = file.path(TARGET_dir, 'number_ss_compilation.csv'))
#######################################
# number of correct steady states 
big_ss_df = data.frame(data_type = NULL, 
                        method = NULL, 
                        style = NULL, 
                        ss_name = NULL, 
                        percent_error = NULL)
for(data_type in data_types) {
  print(data_type)
  methods_list = list.files(file.path("../output/steady_states", data_type))
  for(method in methods_list) {
    style_list = list.files(file.path("../output/steady_states", data_type, method))
    for(style in style_list) {
      compiled_ss = read.csv(file.path("../output/steady_states", data_type, method, style, 'unique_steady_states.csv'), row.names = 1)
      colnames(compiled_ss) = seq(1, ncol(compiled_ss))
      gold_ss = compile_gold_ss(data_type)
      
      for(ss in colnames(gold_ss)) {
        percent_error = min(apply(abs(compiled_ss - gold_ss[, ss]), FUN = sum, MARGIN = 2)) / nrow(gold_ss)
        
        temp_ss_df = data.frame(data_type = data_type, 
                                 method = method, 
                                 style = style, 
                                 ss_name = ss, 
                                 percent_error = percent_error)
        big_ss_df = rbind(big_ss_df, temp_ss_df)
        
      }
    }
  }
}
big_ss_df[grepl('OneCC', big_ss_df$method), 'method'] = 'OneCC_window'
for(data_type in unique(big_ss_df$data_type)) {
  print(data_type)
  temp_ss_df = big_ss_df[big_ss_df$data_type == data_type, ]
  temp_ss_df$percent_correctness = 1 - temp_ss_df$percent_error
  p <- ggplot(temp_ss_df, aes(x=reorder(method, -percent_correctness), y=percent_correctness, fill=ss_name)) +
       geom_bar(stat="identity", position=position_dodge())+theme_minimal() + 
       scale_fill_brewer(palette="Dark2") +     
       ylab("top percent correctness") + 
       xlab("methods") +
       theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(filename = file.path(TARGET_dir, paste0(data_type, "_correctness.png")), plot = p, width = 8, height = 6)
}

write.csv(big_ss_df, file = file.path(TARGET_dir, 'correctness_ss_compilation.csv'))

  
