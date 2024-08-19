library(ggplot2)
library(cowplot)
folder_list = list.dirs("../output/random_networks", recursive = FALSE, full.names = FALSE)

grab_BoolODE_time <- function(file_path) {
  if(file.info(file_path)$size == 0) {
    return(NA)
  }
  
  myData = read.delim(file_path, header = FALSE)
  if(sum(grepl("Simulations took ", myData$V1)) == 0) {
    return(NA)
  } else { 
    time_string = myData$V1[grepl("Simulations took ", myData$V1)]
    time_string = stringr::str_remove(time_string, "Simulations took ")
    time_string = stringr::str_remove(time_string, " s")
    
    time_num = as.numeric(time_string)
    return(time_num)
  }
}

grab_OneCC_time <- function(file_path) {
  myData = read.delim(file_path, header = FALSE)
  time_string = myData[1, ]
  time_string = stringr::str_remove(time_string, "user time is ")
  time_num = as.numeric(time_string)
  return(time_num)
}

big_plot_df = data.frame()
for(temp_trial in folder_list) {
  print(temp_trial)
  edge_config_list = list.dirs(file.path("../output/random_networks", temp_trial), recursive = FALSE, full.names = FALSE)

  OneCC_time_df = data.frame(matrix(nrow = length(edge_config_list), ncol = 5))
  colnames(OneCC_time_df) = c("num_nodes", "num_edges", "time", 'style', 'trial')

  BoolODE_time_df = data.frame(matrix(nrow = length(edge_config_list), ncol = 5))
  colnames(BoolODE_time_df) = c("num_nodes", "num_edges", "time", 'style', 'trial')

  BoolODE2_time_df = data.frame(matrix(nrow = length(edge_config_list), ncol = 5))
  colnames(BoolODE2_time_df) = c("num_nodes", "num_edges", "time", 'style', 'trial')

  i = 1
  for(temp_edge_config in edge_config_list) {
    num_nodes = as.numeric(stringr::str_split(temp_edge_config, "_")[[1]][2])
    num_edges = as.numeric(stringr::str_split(temp_edge_config, "_")[[1]][4])

    focus_dir = file.path("../output/random_networks", temp_trial, temp_edge_config)
    if("BoolODE_output_heaviside.txt" %in% list.files(focus_dir)) {
      BoolODE_heavy_time = grab_BoolODE_time(file.path(focus_dir, "BoolODE_output_heaviside.txt"))
    } else {
      BoolODE_heavy_time = NA
    }

    if(length(BoolODE_heavy_time) > 1) {
      BoolODE_heavy_time = BoolODE_heavy_time[2]
    }
    BoolODE2_time_df[i, 'time'] = BoolODE_heavy_time
    BoolODE2_time_df[i, 'style'] = 'BoolODE Heavyside'
    BoolODE2_time_df[i, 'num_nodes'] = num_nodes
    BoolODE2_time_df[i, 'num_edges'] = num_edges
    BoolODE2_time_df[i, 'trial'] = temp_trial

    if("BoolODE_output.txt" %in% list.files(focus_dir)) {
      BoolODE_time = grab_BoolODE_time(file.path(focus_dir, "BoolODE_output.txt"))
    } else {
      BoolODE_time = NA
    }
    if(length(BoolODE_time) > 1) {
      BoolODE_time = BoolODE_time[2]
    }
    BoolODE_time_df[i, 'time'] = BoolODE_time
    BoolODE_time_df[i, 'style'] = 'BoolODE'
    BoolODE_time_df[i, 'num_nodes'] = num_nodes
    BoolODE_time_df[i, 'num_edges'] = num_edges
    BoolODE_time_df[i, 'trial'] = temp_trial

    if("OneSC_time.txt" %in% list.files(focus_dir)) {
      OneCC_time = grab_OneCC_time(file.path(focus_dir, "OneSC_time.txt"))
    } else { 
      OneCC_time = NA
    }
    
    OneCC_time_df[i, 'time'] = OneCC_time
    OneCC_time_df[i, 'style'] = 'OneSC'
    OneCC_time_df[i, 'num_nodes'] = num_nodes
    OneCC_time_df[i, 'num_edges'] = num_edges
    OneCC_time_df[i, 'trial'] = temp_trial
    i = i + 1
  }
  big_plot_df = rbind(big_plot_df, BoolODE2_time_df, BoolODE_time_df, OneCC_time_df)
}

# to add in the auto-regulation 
#OneCC_time_df$num_edges = OneCC_time_df$num_nodes + OneCC_time_df$num_edges
#BoolODE_time_df$num_edges = BoolODE_time_df$num_nodes + BoolODE_time_df$num_edges

curate_plot_df <- function(big_plot_df, num_nodes = 5) {
  sub_big_plot = big_plot_df[big_plot_df$num_nodes == num_nodes, ]
  return_df = data.frame()
  for(temp_style in unique(sub_big_plot$style)) {
    temp_sub_plot = sub_big_plot[sub_big_plot$style == temp_style, ]

    temp_return_df = data.frame('row.names' = unique(temp_sub_plot$num_edges), 
                                num_edges = unique(temp_sub_plot$num_edges), 
                                med_time = NA, 
                                sd_time = NA, 
                                style = temp_style)
    for(temp_edge in unique(temp_sub_plot$num_edges)) {
      temp_sub_plot2 = temp_sub_plot[temp_sub_plot$num_edges == temp_edge, ]
      temp_sub_plot2 = temp_sub_plot2[!is.na(temp_sub_plot2$time), ]
      temp_return_df[as.character(temp_edge), 'med_time'] = median(temp_sub_plot2$time)
      temp_return_df[as.character(temp_edge), 'sd_time'] = sd(temp_sub_plot2$time)
    }
    return_df = rbind(return_df, temp_return_df)
  }
  return(return_df)
}

for(num_nodes in unique(OneCC_time_df$num_nodes)) {
  plot_df = curate_plot_df(big_plot_df, num_nodes)
  plot_df = plot_df[!is.na(plot_df$med_time), ]
  p<-ggplot(plot_df, aes(x=num_edges, y=med_time, group=style)) +
    geom_line(aes(color=style)) +
    geom_point(aes(color=style)) + 
    geom_ribbon(aes(y = med_time, ymin = med_time - sd_time, ymax = med_time + sd_time, fill = style), alpha = .2) +
    ylim(c(0, 400)) +
    ylab("median time (seconds) - 10 trials") + 
    xlab("number of edges") +
    scale_color_brewer(palette="Dark2") +
    scale_fill_brewer(palette="Dark2") +
    guides(color=guide_legend(title="Methods"), fill='none') +
    ggtitle(paste0(num_nodes, " nodes in the network")) +
    theme_half_open() + 
    theme(text = element_text(size=20)) 
  ggsave(paste0("../output/random_networks/", num_nodes, "_performance.png"), plot = p, width = 8, height = 6)
}

total_plot_df = data.frame()
for(num_nodes in unique(OneCC_time_df$num_nodes)) {
  plot_df = curate_plot_df(big_plot_df, num_nodes)
  plot_df = plot_df[!is.na(plot_df$med_time), ]
  plot_df$network_size = paste0(num_nodes, " nodes networks")
  total_plot_df = rbind(total_plot_df, plot_df)
}

total_plot_df[total_plot_df$style == 'BoolODE Heavyside', 'style'] = 'BoolODE\nSoft-Heaviside' 
total_plot_df[total_plot_df$network_size == '5 nodes networks', 'network_size'] = 'i) 5 nodes networks'
total_plot_df[total_plot_df$network_size == '10 nodes networks', 'network_size'] = 'ii) 10 nodes networks'
total_plot_df[total_plot_df$network_size == '15 nodes networks', 'network_size'] = 'iii) 15 nodes networks'
total_plot_df[total_plot_df$network_size == '20 nodes networks', 'network_size'] = 'iv) 20 nodes networks'

total_plot_df$network_size = factor(total_plot_df$network_size, levels = c('i) 5 nodes networks', 'ii) 10 nodes networks', 'iii) 15 nodes networks', 'iv) 20 nodes networks'))
  
p<-ggplot(total_plot_df, aes(x=num_edges, y=med_time, group=style)) +
  geom_line(aes(color=style), linewidth = 2) +
  geom_point(aes(color=style), size = 3) + 
  geom_ribbon(aes(y = med_time, ymin = med_time - sd_time, ymax = med_time + sd_time, fill = style), alpha = .2) +
  ylim(c(0, 400)) +
  ylab("Median Time (Seconds) - 10 Trials") + 
  xlab("Number of Edges") +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  guides(color=guide_legend(title="Methods"), fill='none') +
  ggtitle('Runtime Comparisons - Single Core') +
  theme_half_open() + 
  facet_wrap(~ network_size, ncol=2, scale = 'free_x') +
  theme(text = element_text(size=30), strip.background = element_blank()) 
ggsave(paste0("../output/random_networks/", "combined_performance.png"), plot = p, width = 11, height = 8)
