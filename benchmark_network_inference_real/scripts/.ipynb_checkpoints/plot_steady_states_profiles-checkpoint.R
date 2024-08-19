library(ggplot2)
library(cowplot)
out_path = '../output/extract_states'

convert_to_long <- function(wide_df) {
  long_df = data.frame()
  for(state in colnames(wide_df)) {
    temp_long_df = data.frame("state" = state, 
                              "bool_exp" = wide_df[, state], 
                              "genes" = rownames(wide_df))
    long_df = rbind(long_df, temp_long_df)
  }
  return(long_df)  
}

color_scheme = c("#f5edf0", "#d64f68")
names(color_scheme) = c("off", 'on')

##### plot out for the CMP maturation #####
temp_path = list.dirs(out_path, recursive = FALSE)[1]
plot_df = data.frame()
for(temp_file in list.files(file.path(temp_path))) {
  wide_df = read.csv(file.path(temp_path, temp_file), row.names = 1)
  temp_long = convert_to_long(wide_df)
  plot_df = rbind(plot_df, temp_long)
}

plot_df[plot_df['bool_exp'] == 0, 'bool_exp'] = 'off'
plot_df[plot_df['bool_exp'] == 1, 'bool_exp'] = 'on'
plot_df$state = factor(x = plot_df$state, levels = c('Monocytes', 'Granulocytes', 'Erythrocytes', 'MK', 'GMP', 'MEP', 'CMP'))
p = ggplot(plot_df, aes(x = genes, y = state, fill = bool_exp)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_manual(values = color_scheme) +
  guides(fill=guide_legend(title="Boolean Status")) +
  coord_fixed() + 
  theme_half_open() + 
  ylab("Cell States") + 
  xlab("Transcription Factors") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = file.path(out_path, "CMP_states.png"), plot = p, width = 8, height = 10)
