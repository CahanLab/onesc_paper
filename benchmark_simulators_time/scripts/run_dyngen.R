library(dyngen)
library(tidyverse)

module_info <- tribble(
  ~module_id, ~basal, ~burn, ~independence,
  "A1", 1, TRUE, 1,
  "A2", 0, TRUE, 1,
  "A3", 1, TRUE, 1,
  "B1", 0, FALSE, 1,
  "B2", 1, TRUE, 1,
  "C1", 0, FALSE, 1,
  "C2", 0, FALSE, 1,
  "C3", 0, FALSE, 1,
  "D1", 0, FALSE, 1,
  "D2", 0, FALSE, 1,
  "D3", 1, TRUE, 1,
  "D4", 0, FALSE, 1,
  "D5", 0, FALSE, 1
)

module_network <- tribble(
  ~from, ~to, ~effect, ~strength, ~hill,
  "A1", "A2", 1L, 10, 2,
  "A2", "A3", -1L, 10, 2,
  "A2", "B1", 1L, 1, 2,
  "B1", "B2", -1L, 10, 2,
  "B1", "C1", 1L, 1, 2,
  "B1", "D1", 1L, 1, 2,
  "C1", "C1", 1L, 10, 2,
  "C1", "D1", -1L, 100, 2,
  "C1", "C2", 1L, 1, 2,
  "C2", "C3", 1L, 1, 2,
  "C2", "A2", -1L, 10, 2,
  "C2", "B1", -1L, 10, 2,
  "C3", "A1", -1L, 10, 2,
  "C3", "C1", -1L, 10, 2,
  "C3", "D1", -1L, 10, 2,
  "D1", "D1", 1L, 10, 2,
  "D1", "C1", -1L, 100, 2,
  "D1", "D2", 1L, 1, 2,
  "D1", "D3", -1L, 10, 2,
  "D2", "D4", 1L, 1, 2,
  "D4", "D5", 1L, 1, 2,
  "D3", "D5", -1L, 10, 2
)

expression_patterns <- tribble(
  ~from, ~to, ~module_progression, ~start, ~burn, ~time,
  "sBurn", "sA", "+A1,+A2,+A3,+B2,+D3", TRUE, TRUE, 60,
  "sA", "sB", "+B1", FALSE, FALSE, 60,
  "sB", "sC", "+C1,+C2|-A2,-B1,+C3|-C1,-D1,-D2", FALSE, FALSE, 80,
  "sB", "sD", "+D1,+D2,+D4,+D5", FALSE, FALSE, 120,
  "sC", "sA", "+A1,+A2", FALSE, FALSE, 60
)

expression_patterns <- tribble(
  ~from, ~to, ~module_progression, ~start, ~burn, ~time,
  "sBurn", "sA", "+A1,+A2,+A3,+B2,+D3", TRUE, TRUE, 60,
  "sA", "sB", "+B1", FALSE, FALSE, 60,
  "sB", "sC", "+C1,+C2|-A2,-B1,+C3|-C1,-D1,-D2", FALSE, FALSE, 80,
  "sB", "sD", "+D1,+D2,+D4,+D5", FALSE, FALSE, 120,
  "sC", "sA", "+A1,+A2", FALSE, FALSE, 60
)

backbone <- backbone(
  module_info = module_info,
  module_network = module_network,
  expression_patterns = expression_patterns
)

config <- initialise_model(
  backbone = backbone,
  num_tfs = nrow(backbone$module_info),
  num_targets = 500,
  num_hks = 500,
  simulation_params = simulation_default(
    experiment_params = simulation_type_wild_type(num_simulations = 100),
    total_time = 600
  ),
  verbose = FALSE
)

plot_backbone_modulenet(config)
