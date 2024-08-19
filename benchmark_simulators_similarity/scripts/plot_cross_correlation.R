library(ggplot2)
library(stringr)
library(cowplot)

file_list = list.files("../output/cross_corr/")
file_list = file_list[grepl(".png", file_list) == FALSE]

plot_df = data.frame()
for(temp_file in file_list) {
  data_type = stringr::str_split(temp_file, "_")[[1]][1]
  cross_cor_df = read.csv(file.path("../output/cross_corr/", temp_file))
  if(temp_file %in% c("dyn-CY_cross_corr.csv", "dyn-LI_cross_corr.csv" ) == FALSE) {
    cross_cor_df = cross_cor_df[cross_cor_df$sim_ss == cross_cor_df$train_ss, ]
  }
  cross_cor_df$data_type = data_type
  plot_df = rbind(plot_df, cross_cor_df)
}

plot_df[plot_df$cross_corr > 1, 'cross_corr'] = 1
plot_df = plot_df[plot_df$data_type != 'dyn-LI', ]
p<-ggplot(plot_df, aes(x=data_type, y=cross_corr)) +
  geom_boxplot(lwd=1) + 
  ylab("Cross Correlation of \n Gene Expression Dynamics") +
  xlab("Single-cell Data Types") +
  ylim(c(0, 1)) +
  theme_half_open() + 
  theme(text = element_text(size=20), 
        axis.text = element_text(size=20))
ggsave(filename = '../output/cross_corr/cross_correlation.png', plot = p, width = 10, height = 4.5)

p<-ggplot(plot_df, aes(x=data_type, y=lag)) +
  geom_boxplot(lwd=1) + 
  ylab("Lag") +
  xlab("Single-cell Data Structures") +
  ylim(c(0, 1)) +
  theme_bw() + 
  theme(text = element_text(size=20), 
        axis.text = element_text(size=20))
ggsave(filename = '../output/cross_corr/lag.png', plot = p, width = 10, height = 4.5)
