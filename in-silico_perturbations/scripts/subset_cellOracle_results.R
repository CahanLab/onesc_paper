library(ggplot2)
library(cowplot)
library(ggrepel)
library(tidyr)


cellOracle = read.csv("../literature_validations/ CellOracle/hemapoesis.csv", row.names = 1)

##### load the plots #####
##### set up the working directory #####
out_path = '../output/ss_proportion_plot'

##### plot out all the perturbation ######
all_files = list.files("../output/ss_proportion/")
for(temp_file in all_files) {
  proportion_df = read.csv(file.path("../output/ss_proportion", temp_file))
  
  long_prop_df = proportion_df %>%
    pivot_longer(-X)
  colnames(long_prop_df) = c("TF", "states", "counts")
  if(grepl('day_4_EB', temp_file)) {
    plot_title = 'EB'  
  } else if(grepl('myeloid', temp_file)) {
    plot_title = 'Mouse Hematopoietic'  
  } else if(grepl("t_cell", temp_file)) {
    plot_title = 'T-cells'  
  }
  
  if(grepl('knockout', temp_file)) {
    perturb_type = 'knockout'
  } else if(grepl('overexpression', temp_file)) {
    perturb_type = 'overexpression'
  }
  
  long_prop_df$cat = NA
  long_prop_df[long_prop_df$TF == 'wild_type', 'cat'] = 'WT'
  if(perturb_type == 'knockout') {
    long_prop_df[long_prop_df$TF != 'wild_type', 'cat'] = 'KO'
  } else if(perturb_type == 'overexpression') {
    long_prop_df[long_prop_df$TF != 'wild_type', 'cat'] = 'OE'
  }
  
  p = ggplot(long_prop_df, aes(fill=states, y=counts, x=TF)) + 
    geom_bar(position="stack", stat="identity") + 
    theme_half_open() + 
    scale_fill_brewer(palette = "Set2") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    facet_grid(
      rows = vars(cat),
      scales = "free_y",
      space = "free_y",
      switch = "x"
    ) + 
    theme(
      panel.spacing = unit(x = 0.2, units = "lines"),
      strip.background = element_blank(), 
      strip.text.y = element_blank()
    ) +
    ylab(paste0('# / 100 simulations')) +
    xlab(paste0(perturb_type, " of TFs")) +
    coord_flip() + 
    guides(fill=guide_legend(title="Final States"))
  
}

##### this is the data #####
cat_temp = unique(long_prop_df$TF)[unique(long_prop_df$TF) != 'wild_type']
sort(cat_temp, decreasing = TRUE)

cellOracle = cellOracle[sort(cat_temp, decreasing = TRUE), ]
write.csv(cellOracle, file = '../output/selected_CellOracle/celloracle_data.csv')
