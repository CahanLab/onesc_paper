library(stringr)
library(BoolNet)
library(ggplot2)
library(viridis)
library(cowplot)
library(dplyr)

data_type_run = list.dirs("../output/Boolean_sim", recursive = FALSE, full = FALSE)

check_num <- function(ss_df, temp_states_list, data_type) {
    if(data_type == 'dyn-CY') {
      expected_ss = 0
    } else { 
      expected_ss = sum(grepl("steady", names(temp_states_list)))
    }
    sim_ss = ncol(ss_df)
    if(length(ss_df[is.na(ss_df) == FALSE, ]) == 0) {
      return(0 - expected_ss)
    } else {
      return(sim_ss - expected_ss)
    }
}

check_avg_agreement <- function(ss_df, temp_states_list) {
    ss_names = names(temp_states_list)[grepl("steady", names(temp_states_list))]
    total_percent = 0
    if(ncol(ss_df) == 0) {
      return(NA)
    }
    for(ss_name in ss_names) {
        gold_ss = temp_states_list[[ss_name]]
        gold_ss = gold_ss[rownames(ss_df)]
        diff_ss_df = abs(ss_df - gold_ss)
        percent_agreement = 1 - min(apply(diff_ss_df, MARGIN = 2, FUN = sum)) / nrow(ss_df)
        total_percent = total_percent + percent_agreement
    } 
    return(max(0, total_percent / length(ss_names)))
}

plot_df = data.frame()

for(data_type in data_type_run) {
    print(data_type)
    method_list = list.dirs(file.path("../output/Boolean_sim", data_type), recursive = FALSE, full = FALSE)
    steady_states = read.csv(file.path("../output/golden_Boolean_sim", data_type, "unique_steady_states_profiles.csv"), row.names = 1)
    temp_states_list = list()
    if (data_type == 'dyn-CY') {
      states = c(-1, -1, -1, -1, -1, -1)
      names(states) = c("g1", "g2", "g3", "g4", "g5", "g6")
      temp_states_list[['steady_0']] = states
    } else if(ncol(steady_states) == 1) { 
      state_vector = as.vector(steady_states[, 1])
      names(state_vector) = rownames(steady_states)
      temp_states_list[[paste0('steady_', '1')]] = state_vector
    } else { 
      steady_states = steady_states[, c(lapply(steady_states, MARGIN = 2, FUN = sum) > 0)]
      for(i in seq(1, ncol(steady_states))) {
        state_vector = as.vector(steady_states[, i])
        names(state_vector) = rownames(steady_states)
        temp_states_list[[paste0('steady_', i)]] = state_vector
      }
    }

    for(method in method_list) {
        ss_df = read.csv(file.path("../output/Boolean_sim", data_type, method, "unique_steady_states_profiles.csv"), row.names = 1)
        ss_df = ss_df[, colSums(is.na(ss_df))<nrow(ss_df), drop = FALSE]
        temp_plot = data.frame("data_type" = data_type, 
                               "method" = method, 
                               "steady_states_num" = check_num(ss_df, temp_states_list, data_type), 
                               "average_steady_states_agreement" = check_avg_agreement(ss_df, temp_states_list))
        plot_df = rbind(plot_df, temp_plot)
    }   
}

plot_df[plot_df$method == "run_OneSC", "method"] = "OneSC"

plot_df$average_steady_states_agreement = round(plot_df$average_steady_states_agreement, 2)

plot_df$method = toupper(plot_df$method)

cat_ss = plot_df %>% 
  group_by(method) %>%
  summarise(MeanOfBoolean = mean(average_steady_states_agreement, na.rm = TRUE)) %>% 
  arrange(MeanOfBoolean)

cat_ss = plot_df %>% 
  group_by(method) %>% 
  summarize(mean_of_ss = mean(abs(steady_states_num))) %>% 
  arrange(mean_of_ss)
cat_ss = read.csv(file = '../output/Performance_F1/category_ordering.csv', row.names = 1)

plot_df = plot_df[plot_df$method != 'PCOR', ]

plot_df$method = factor(plot_df$method, levels = cat_ss$method)

# not significant 
t.test(plot_df[plot_df$method == 'ONESC', 'average_steady_states_agreement'], plot_df[plot_df$method == 'GENIE3', 'average_steady_states_agreement'])
t.test(plot_df[plot_df$method == 'ONESC', 'average_steady_states_agreement'], plot_df[plot_df$method == 'SCODE', 'average_steady_states_agreement'])

p = ggplot(plot_df, aes(data_type, method, fill= average_steady_states_agreement)) + 
  geom_tile() +
  geom_text(aes(label = average_steady_states_agreement)) + 
  ylab("Methods") +
  xlab("Data Types") +
  scale_fill_viridis(discrete=FALSE, alpha = 0.8) + 
  guides(fill=guide_legend(title="")) +
  theme_half_open() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None")
        
ggsave(filename = file.path("../output/Boolean_sim/agreement.png"), plot = p, height = 6, width = 5)

p = ggplot(plot_df, aes(data_type, method, fill= steady_states_num)) + 
  geom_tile() +
  geom_text(aes(label = steady_states_num)) + 
  ylab("Methods") +
  xlab("Data Types") + 
  guides(fill=guide_legend(title="Number of \nSteady States \nDifferences")) +
  scale_fill_gradient2(low = "steelblue",
                      mid = 'white',
                      high = "salmon",
                      midpoint = 0,
                      guide = "colorbar") + 
  theme_half_open() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "None")

ggsave(filename = file.path("../output/Boolean_sim/num_ss.png"), plot = p, height = 6, width = 5)

# make the plots with annotations 
p = ggplot(plot_df, aes(data_type, method, fill= average_steady_states_agreement)) + 
  geom_tile() +
  geom_text(aes(label = average_steady_states_agreement)) + 
  ylab("Methods") +
  xlab("Data Types") +
  scale_fill_viridis(discrete=FALSE, alpha = 0.8) + 
  guides(fill=guide_legend(title="")) +
  theme_half_open() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "left")

ggsave(filename = file.path("../output/Boolean_sim/agreement_annotate.png"), plot = p, height = 6, width = 4)

p = ggplot(plot_df, aes(data_type, method, fill= steady_states_num)) + 
  geom_tile() +
  geom_text(aes(label = steady_states_num)) + 
  ylab("Methods") +
  xlab("Data Types") + 
  guides(fill=guide_legend(title="Number of \nSteady States \nDifferences")) +
  scale_fill_gradient2(low = "steelblue",
                       mid = 'white',
                       high = "salmon",
                       midpoint = 0,
                       guide = "colorbar") + 
  theme_half_open() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "left")

ggsave(filename = file.path("../output/Boolean_sim/num_ss_annotate.png"), plot = p, height = 6, width = 4)
