#' Computing and testing AT1 and AT2, and plotting sub-network
#'
#' @inheritParams compute_edge_t_stat
#' @inheritParams compute_AT
#'
#' @return The same return as [compute_AT()]. A list of the following: \cr
#' subnet_label: label of the subnet; num_edges: number of edges in the sub-network; \cr
#' vertex_idx_selected, edge_pair_selected: the sub-network selection arguments provided by you; \cr
#' t_stat: t statistics of edges in the sub-network, t_stat_perc: percentiles of these t statistics among all edges in the original network;\cr
#' AT1: Value of AT1; pval_AT1: p-value of AT1; \cr
#' AT2: Value of AT2; pval_AT2: p-value of AT2)
#' @export
#'
#' @examples
#' network <- read_graph(here::here("demo/network_info", "network"), format = "edgelist")
#' data_type1 <- readr::read_csv(here::here("demo/data", paste0("LGG", ".csv")))
#' data_type2 <- readr::read_csv(here::here("demo/data", paste0("GBM", ".csv")))
#' edge_pair_selected <- c("1-5", "1-10", "5-10", "8-10", "2-11", "3-11", "4-11")
#'
#' # compute AT's
#' AT_res <- analysis_ExprNet(data_type1 = data_type1, data_type2 = data_type2, network = network,
#'                            edge_pair_selected = edge_pair_selected, type1_name = "LGG", type2_name = "GBM",
#'                            subnet_label = "Demo_GO0006306_DNA_methylation(LGG-GBM)",
#'                            save_edge_res = TRUE, save_plot = TRUE, save_dir = here::here("demo"),
#'                            vertex.label.cex = 1, vertex.size = 10, edge.width = 7)

analysis_ExprNet <- function(data_type1, data_type2, network,
                             vertex_idx_selected = NULL, edge_pair_selected = NULL,
                             alpha_t_test = 0.05, AT2_perm_test = TRUE, num_perm = 500, num_cores = parallel::detectCores(),
                             type1_name = "Type1", type2_name = "Type2",
                             subnet_label = "subnet", plot_subnet = TRUE,
                             save_edge_res = FALSE, save_plot = FALSE, save_dir = here::here(), ...) {

  edge_res <- compute_edge_t_stat(data_type1, data_type2, network, type1_name, type2_name, alpha_t_test,
                                  save_edge_res, save_dir)  # compute edge length and t statistics
  AT_res <- compute_AT(edge_res$edge_t_stat, edge_res$edge_dist_mat, network,
                       vertex_idx_selected, edge_pair_selected,
                       AT2_perm_test, num_perm, num_cores,
                       type1_name, type2_name,
                       subnet_label, plot_subnet, save_plot, save_dir, ...)  # compute AT1 & AT2

  return (AT_res)
}


