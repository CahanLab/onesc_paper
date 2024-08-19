library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

##### set up the working directory #####
out_path = '../output/ss_proportion_plot'
dir.create(out_path, recursive = TRUE)

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
      strip.text.y = element_blank(), 
      text = element_text(size = 20), 
      axis.text = element_text(size = 18)
    ) +
    ylab(paste0('# / 200 simulations')) +
    xlab(paste0(perturb_type, " of TFs")) +
    coord_flip() + 
    guides(fill=guide_legend(title="Final States"))
 
  ggsave(filename = file.path(out_path, paste0(perturb_type, "_", plot_title, "_plots.png")), plot = p, width = 6, height = 8) 
}

##### select the plots that are important ######
all_files = list.files("../output/ss_proportion/")
compiled_long_df = data.frame()
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
  
  compiled_long_df = rbind(compiled_long_df, long_prop_df)
}

compiled_long_df = compiled_long_df[!duplicated(compiled_long_df), ]
compiled_long_df$category = NA
compiled_long_df$category = paste0(compiled_long_df$TF, "_", compiled_long_df$cat)
compiled_long_df[compiled_long_df$TF == 'wild_type', 'category'] = 'wild_type'

interesting_cond = c("Cebpa_KO", 
                     "wild_type", 
                     "Cebpe_KO", 
                     "Gata1_KO", 
                     "Gfi1b_KO", 
                     "Irf8_KO")
sub_compiled_df = compiled_long_df[compiled_long_df$category %in% interesting_cond, ]

p = ggplot(sub_compiled_df, aes(fill=states, y=counts, x=category)) + 
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
    strip.text.y = element_blank(), 
    legend.position = 'bottom', 
    text = element_text(size = 30), 
    axis.text = element_text(size = 30), 
    plot.margin = margin(r = 500)
  ) +
  ylab(paste0('# / 200 simulations')) +
  xlab(paste0("Perturbations")) +
  coord_flip() + 
  guides(fill=guide_legend(title="Final States"))
ggsave(filename = file.path(out_path, paste0("selected_markers_plots.png")), plot = p, width = 15, height = 15) 

