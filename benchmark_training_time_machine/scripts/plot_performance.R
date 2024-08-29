library(ggplot2)
library(dplyr)
library(cowplot)

performance_df = read.csv("../output/compiled_runtime.csv", row.names = 1)

performance_df$label = paste0(performance_df$machine_type, "\n", performance_df$cores, " cores ", performance_df$OS)

p = ggplot(performance_df, aes(x = reorder(label, cores), y = runtime)) +
  geom_bar(stat = "identity") +
  ggtitle("Runtime for inferring myeloid networks") +
  xlab("Machine Type") +
  ylab("Runtime (seconds)") +
  theme_cowplot()

ggsave('../output/runtime_plot.png', plot = p, width = 9, height = 5)
