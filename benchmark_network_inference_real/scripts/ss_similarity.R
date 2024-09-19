library(ggplot2)
library(stringr)
library(dplyr)
library(viridis)
library(cowplot)

output_path = "../output/sim_ss_comparison"
dir.create(output_path, recursive = TRUE)
  
datatype = 'myeloid_progenitors_full'
ss_path = "../output/extract_states/myeloid_progenitors_full/"

curate_ss <- function(ss_path) {
  state_df = data.frame()
  for(temp_file in list.files(ss_path)) {
    temp_state_df = read.csv(file.path(ss_path, temp_file), row.names = 1)
    if(ncol(state_df) == 0) {
      state_df = temp_state_df
    } else {
      state_df = cbind(state_df, temp_state_df)
    }
  }
  state_df = state_df[!duplicated(as.list(state_df))]
  return(state_df)
}

calc_max_sim <- function(reference_state, sim_states) {
  all_percent_agree = c()
  for(unique_ss in colnames(sim_states)) {
    num_diff = sum(abs(sim_states[, unique_ss] - reference_state))
    percent_agree = 1 - (num_diff / length(reference_state))
    all_percent_agree = c(all_percent_agree, percent_agree)
  }
  return(max(all_percent_agree))
}

state_df = curate_ss(ss_path)

plot_df = data.frame()
for(temp_method in list.files(file.path("../output/steady_states/", datatype))) {
  print(temp_method)
  sim_states = read.csv(file.path("../output/steady_states/", datatype, temp_method, "density/unique_steady_states.csv"), row.names = 1)
  for(temp_state in colnames(state_df)) {
    print(temp_state)
    reference_state = state_df[, temp_state]
    temp_plot_df = data.frame('methods' = c(temp_method), 
                              'states' = c(temp_state), 
                              'percent_agreement' = calc_max_sim(reference_state, sim_states))
    plot_df = rbind(plot_df, temp_plot_df)
  }
}

##### get the iqcell steady states ######
pt_paths = list.dirs("../output/iQcell_simulations/myeloid_progenitors_full/", recursive = FALSE, full.names = FALSE)
multiple_read_df <- function(x) {
  return(read.csv(file.path("../output/iQcell_simulations/myeloid_progenitors_full/", x, "in_out/ABNfiles/9.2_profilesOutput.csv"), stringsAsFactors=FALSE))  
}

all_iqcell_states = lapply(pt_paths, multiple_read_df) %>%
  Reduce(function(x, y) rbind(x, y), .) %>%
  t()
all_iqcell_states = all_iqcell_states[rownames(state_df), ]
all_iqcell_states = as.data.frame(all_iqcell_states)
all_iqcell_states =  all_iqcell_states %>% mutate_all(as.numeric)
#all_iqcell_states = as.numeric(all_iqcell_states)
for(temp_state in colnames(state_df)) {
  print(temp_state)
  reference_state = state_df[, temp_state]
  temp_plot_df = data.frame('methods' = c('iQcell'), 
                            'states' = c(temp_state), 
                            'percent_agreement' = calc_max_sim(reference_state, all_iqcell_states))
  plot_df = rbind(plot_df, temp_plot_df)
}

plot_df = plot_df[plot_df$states %in% c('Erythrocytes', 'Granulocytes', 'Monocytes', 'MK'), ]
plot_df$percent_agreement = round(plot_df$percent_agreement, 2)
plot_df[plot_df$methods == 'run_OneSC', 'methods'] = 'OneSC'
plot_df$methods = toupper(plot_df$methods)

cat_ss = plot_df %>% 
  dplyr::group_by(methods) %>% 
  summarise_at(vars(percent_agreement), mean) %>% 
  arrange(percent_agreement)
plot_df$methods = factor(plot_df$methods, levels = cat_ss$methods)
write.csv(plot_df, file.path(output_path, 'plot_df.csv'))

p = ggplot(plot_df, aes(x = methods, y = states, fill= percent_agreement)) + 
  geom_tile() +
  geom_text(aes(label = percent_agreement)) + 
  ylab("Terminal States") +
  xlab("Methods") +
  scale_fill_viridis(discrete=FALSE) + 
  guides(fill=guide_legend(title="Top \nTerminal States \nSimilarity")) +
  theme_half_open() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  coord_flip()
ggsave(filename = file.path(output_path, paste0(datatype, ".png")), plot = p, width = 5, height = 8)

p = ggplot(plot_df, aes(x = methods, y = states, fill= percent_agreement)) + 
  geom_tile() +
  geom_text(aes(label = percent_agreement)) + 
  ylab("Terminal States") +
  xlab("Methods") +
  scale_fill_viridis(discrete=FALSE) + 
  guides(fill=guide_legend(title="Top \nTerminal \nStates \nSimilarity")) +
  theme_half_open() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  coord_flip()
ggsave(filename = file.path(output_path, paste0(datatype, "_copy.png")), plot = p, width = 5, height = 8)
