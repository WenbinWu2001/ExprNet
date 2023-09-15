## demo

rm(list = ls())
devtools::install(dependencies = "never")

library(ExprNet)


network <- read_graph(here::here("demo/network_info", "network"), format = "edgelist")
data_type1 <- readr::read_csv(here::here("demo/data", paste0("LGG", ".csv")))
data_type2 <- readr::read_csv(here::here("demo/data", paste0("GBM", ".csv")))
edge_pair_selected <- c("1-8", "1-15", "2-16", "3-16", "5-10", "5-16", "8-11", "8-13", "8-14", "8-15", "13-15")

analysis_ExprNet(data_type1, data_type2, network,
                 edge_pair_selected = edge_pair_selected,
                 type1_name = "LGG", type2_name = "GBM",
                 subnet_label = "Demo_GO0006306_DNA_methylation(LGG-GBM)",
                 save_edge_res = TRUE, save_plot = TRUE, save_dir = here::here("demo"),
                 vertex.label.cex = 1, vertex.size = 10, edge.width = 7)


# ## Stepwise Demo
# save_dir <- here::here("demo")
# network <- read_graph(here::here("demo/network_info", "network"), format = "edgelist")
# data_type1 <- readr::read_csv(here::here("demo/data", paste0("LGG", ".csv")))
# data_type2 <- readr::read_csv(here::here("demo/data", paste0("GBM", ".csv")))
# vertex_idx_selected <- c(1, 2, 3, 5, 8, 10, 11, 13, 14, 15, 16)
# edge_pair_selected <- c("1-8", "1-15", "2-16", "3-16", "5-10", "5-16", "8-11", "8-13", "8-14", "8-15", "13-15")
#
# ## use of compute_edge_t_stat
# res <- compute_edge_t_stat(data_type1, data_type2, network, type1_name = "LGG", type2_name = "GBM", save_edge_res = TRUE, save_dir = save_dir)
# edge_t_stat <- res$edge_t_stat
# edge_dist_mat <- res$edge_dist_mat
#
# ## use of compute_AT
# compute_AT(edge_t_stat, edge_dist_mat, network, vertex_idx_selected = vertex_idx_selected)  # select with vertex indices
# compute_AT(edge_t_stat, edge_dist_mat, network, edge_pair_selected = edge_pair_selected, save_plot = TRUE, save_dir = here::here("demo"),
#            subnet_label = "Demo_GO0006306_DNA_methylation(LGG-GBM)", vertex.label.cex = 1, vertex.size = 10, edge.width = 7)  # select with edge pairs
