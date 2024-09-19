library(stringr)
library(ggplot2)
library(cowplot)
library(ggbreak)
library(dplyr)
library(cowplot)

data_types = list.dirs("../Beeline_benchmark/outputs/example/", full.names = FALSE, recursive = FALSE)

##### load all the other GRN inference methods #####
big_time_df = data.frame()
for(data_type in data_types) { 
    methods_list = list.dirs(file.path("../Beeline_benchmark/outputs/example", data_type), full.names = FALSE, recursive =  FALSE)
    temp_time_df = data.frame(row.names = methods_list) 
    temp_time_df$methods = methods_list
    temp_time_df$data_type = data_type
    temp_time_df$user_time = NA
    for(method in methods_list) { 
      time_files = list.files(file.path("../Beeline_benchmark/outputs/example", data_type, method))
      time_files = time_files[grep("time", time_files)]
      
      use_time = 0
      for(time_file in time_files) {
        res <- readLines(file.path("../Beeline_benchmark/outputs/example", data_type, method, time_file))
        for(file_line in res) { 
          file_line = toupper(file_line)
          if(grepl("USER TIME", file_line)) { 
            string_list = stringr::str_split(file_line, " ")[[1]]
          }
        }
        use_time = use_time + as.numeric(string_list[length(string_list)])
      }
      temp_time_df[method, 'user_time'] = use_time
    }
    big_time_df = rbind(big_time_df, temp_time_df)
}

##### load in the runtime for OneSC #####
OneSC_style = list.dirs("../Beeline_benchmark", full.names = FALSE, recursive = FALSE)
OneSC_style = OneSC_style[grep('OneSC', OneSC_style)]

option_list = list()
option_list[['run_OneSC']] = 'OneSC'

for(temp_style in OneSC_style) {
  data_types = list.dirs(paste0("../Beeline_benchmark/", temp_style, "/"), full.names = FALSE, recursive = FALSE)
  data_types = data_types[data_types != ".ipynb_checkpoints"]
  temp_time_df = data.frame(row.names = data_types) 
  temp_time_df$data_type = data_types
  temp_time_df$methods = option_list[[temp_style]]
  temp_time_df$user_time = NA
  
  for(data_type in data_types) { 
    res = readLines(file.path("../Beeline_benchmark", temp_style, data_type, 'time.txt'))
    for(file_line in res) { 
      file_line = toupper(file_line)
      if(grepl("USER TIME", file_line)) { 
        string_list = stringr::str_split(file_line, " ")[[1]]
        temp_time_df[data_type, 'user_time'] = string_list[length(string_list)]
      }
    }
  }
  
  big_time_df = rbind(big_time_df, temp_time_df)
  big_time_df$user_time = as.numeric(big_time_df$user_time)
}

##### load in the runtime for iQcell #####
iQcell_style = list.dirs("../Beeline_benchmark", full.names = FALSE, recursive = FALSE)
iQcell_style = iQcell_style[grep('iQcell', iQcell_style)]

option_list = list()
option_list[['run_iQcell']] = 'iQcell'

for(temp_style in iQcell_style) {
  data_types = list.dirs(paste0("../Beeline_benchmark/", temp_style, "/"), full.names = FALSE, recursive = FALSE)
  data_types = data_types[data_types != ".ipynb_checkpoints"]
  temp_time_df = data.frame(row.names = data_types) 
  temp_time_df$data_type = data_types
  temp_time_df$methods = option_list[[temp_style]]
  temp_time_df$user_time = NA
  
  for(data_type in data_types) { 
    res = readLines(file.path("../Beeline_benchmark", temp_style, data_type, 'time.txt'))
    for(file_line in res) { 
      file_line = toupper(file_line)
      if(grepl("USER TIME", file_line)) { 
        string_list = stringr::str_split(file_line, " ")[[1]]
        temp_time_df[data_type, 'user_time'] = string_list[length(string_list)]
      }
    }
  }
  
  big_time_df = rbind(big_time_df, temp_time_df)
  big_time_df$user_time = as.numeric(big_time_df$user_time)
}

big_time_df$OneSC_True = NA
big_time_df[big_time_df$methods == 'OneSC', 'OneSC_True'] = TRUE
big_time_df[big_time_df$methods != 'OneSC', 'OneSC_True'] = FALSE

big_time_df %>% 
  group_by(methods) %>% 
  summarize(avg_value = mean(user_time)) %>% 
  arrange(avg_value)

big_time_df$methods = toupper(big_time_df$methods)
p <- ggplot(big_time_df, aes(x=data_type, y=user_time, fill = methods)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  #scale_fill_manual(values = RColorBrewer::brewer.pal(12, 'Set3')) + 
  geom_text(aes(label = ifelse(OneSC_True, "*", "")), 
          position = position_dodge(width = .9), vjust = 1, size = 20 / .pt) +
  ylab("User Time (Sec)") + 
  xlab("Data Types") + 
  theme_half_open() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    facet_grid(
    cols = vars(data_type),
    scales = "free_x",
    space = "free_x",
    switch = "x"
  ) + 
  theme(
    panel.spacing = unit(x = 0.2, units = "lines"),
    strip.background = element_blank(), 
    axis.text.x=element_blank(), #remove x axis labels
    axis.ticks.x=element_blank() #remove x axis ticks
  ) + 
  scale_y_break(c(2000, 4000), scales = 1) 

dir.create("../output/time_comparison")
ggsave("../output/time_comparison/run_time_comparison.png", plot = p, width = 10, height = 5)

##### plot out number of nodes vs runtime ######
OneSC_df = big_time_df[big_time_df$methods == 'ONESC', ]
num_nodes = c()
for(tmp_datatype in OneSC_df$data_type) {
  network = read.csv(paste0("../Beeline_benchmark/inputs/example/", tmp_datatype, '/refNetwork.csv'))
  all_genes = unique(c(network$Gene1, network$Gene2))
  num_nodes = c(num_nodes, length(all_genes))
}
OneSC_df$num_nodes = num_nodes

p = ggplot(data = OneSC_df, aes(x = num_nodes, y = user_time)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = 'black') + 
  theme_cowplot() + 
  ylab("Runtime (seconds)") + 
  xlab("Number of network genes") + 
  annotate(
    "text",
    x = 12,  # X position of the text
    y = 1000,  # Y position of the text
    label = paste("Pearson r =", round(cor(OneSC_df$num_nodes, OneSC_df$user_time), 2)),  # Text label with correlation
    size = 5,  # Font size
    color = "blue"  # Text color
  )+
  ggtitle("OneSC runtime of BEELINE data on 16 cores machine")

ggsave("../output/time_comparison/run_time_machine.png", plot = p, width = 10, height = 5)

  