library(stringr)
library(ggplot2)
library(cowplot)

data_types = list.files("../output/max_cross_correlation")

TARGET_dir = '../output/max_cross_correlation_plot'
dir.create(TARGET_dir)

big_df = data.frame()
for(data_type in data_types) {
  print(data_type)
  methods_list = list.files(file.path("../output/max_cross_correlation", data_type))
  for(method in methods_list) {
    style_list = list.files(file.path("../output/max_cross_correlation", data_type, method))
    for(style in style_list) {
      compiled_cc = read.csv(file.path("../output/max_cross_correlation", data_type, method, style, 'max_cross_correlation.csv'))
      compiled_cc$X = NULL
      compiled_cc$data_type = data_type
      compiled_cc$method = method 
      compiled_cc$style = style
      big_df = rbind(big_df, compiled_cc)
    }
  }
  big_df[big_df$method == 'run_OneCC_full', 'method'] = 'OneCC'
  p = ggplot(big_df, aes(x = reorder(method, -cross_corr), y = cross_corr, fill = lineage_type)) + 
    geom_boxplot() + 
    theme_half_open() + 
    xlab("methods") + 
    ylab("gene expression cross-correlation") + 
    guides(fill=guide_legend(title="Trajectories")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(filename = file.path(TARGET_dir, paste0(data_type, ".png")), plot = p, width = 10, height = 4.5)
}


##### look at min dtw #####
big_df = data.frame()
for(data_type in data_types) {
  print(data_type)
  methods_list = list.files(file.path("../output/min_dtw", data_type))
  for(method in methods_list) {
    style_list = list.files(file.path("../output/min_dtw", data_type, method))
    for(style in style_list) {
      compiled_cc = read.csv(file.path("../output/min_dtw", data_type, method, style, 'max_dtw.csv'))
      compiled_cc$X = NULL
      compiled_cc$data_type = data_type
      compiled_cc$method = method 
      compiled_cc$style = style
      big_df = rbind(big_df, compiled_cc)
    }
  }
  big_df[big_df$method == 'run_OneCC_full', 'method'] = 'OneCC'
  p = ggplot(big_df, aes(x = reorder(method, dtw), y = dtw, fill = lineage_type)) + 
    geom_boxplot() + 
    theme_half_open() + 
    xlab("methods") + 
    ylab("gene expression DTW") + 
    guides(fill=guide_legend(title="Trajectories")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(filename = file.path(TARGET_dir, paste0(data_type, "_dtw.png")), plot = p, width = 10, height = 4.5)
}


for(style in style_list) {
  if(style == 'orphan') {
    temp_big_df = big_df[big_df$style == style | big_df$method == 'OneCC_pyEpoch', ]
    
  } else {
    temp_big_df = big_df[big_df$style == style, ]
  }
  temp_big_df = temp_big_df[temp_big_df$cross_corr != -999.00000000, ]
  p = ggplot(temp_big_df, aes(x=reorder(method, cross_corr), y=cross_corr, fill=data_type)) +
    geom_boxplot() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    ylab("cross correlation") + 
    xlab("methods")
  ggsave(filename = file.path(TARGET_dir, style, 'cross_correlation.png'), plot = p, width = 8, height= 6)
  for(data_type in unique(temp_big_df$data_type)) {
    sub_temp_big_df = temp_big_df[temp_big_df$data_type == data_type, ]
    p = ggplot(sub_temp_big_df, aes(x=reorder(method, cross_corr), y=cross_corr, fill=lineage_type)) +
      geom_boxplot() + theme_bw() + 
      scale_fill_brewer(palette="Dark2") +     
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      ylab("cross correlation") + 
      xlab("methods")
    ggsave(filename = file.path(TARGET_dir, style, paste0(data_type, '_cross_correlation.png')), plot = p, width = 8, height= 6)
    
  }
}

