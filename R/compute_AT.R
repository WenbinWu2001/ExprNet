#' Given t test results, computing and testing AT1 and AT2, and plotting sub-network
#'
#' @description Given the results from [compute_edge_t_stat()], compute AT1 and AT2.
#' Calculate the p-value of AT1 using the Irwin-Hall (or normal) null distribution and the p-value of AT2 using permutation test.
#' Plot and save the sub-network.
#' A red edge indicates that the t statistic of this edge is larger than 0, i.e. mean of its length is larger in phenotype 1 than in phenotype 2. A blue edge indicates the opposite.
#'
#' @inheritParams compute_edge_t_stat
#' @param edge_t_stat The result from [compute_edge_t_stat()].
#' @param edge_dist_mat The result from [compute_edge_t_stat()].
#' @param vertex_idx_selected A numeric vector of vertex indices that you want to include in the sub-network. All edges on the network between these vertices will be included.
#' @param edge_pair_selected A character vector of edges that you want to include in the sub-network, each element in the form *"vertex1-vertex2"* (e.g. "1-3").
#' @param AT2_perm_test Whether to conduct permutation test for AT2. If *FALSE*, then p-value of AT2 will be returned as NA.
#' @param num_perm Number of permutations in the permutation test.
#' @param num_cores Number of cores to register for parallel computing in permutation test. Should be no larger than the number of available cores on the computer.
#' @param subnet_label A label for the sub-network, used for naming the files and setting the title in the plot.
#' @param plot_subnet Whether to plot the sub-network.
#' @param save_plot Whether to save the sub-network.
#' @param ... Parameters for plots. For example, you may specify the size for the vertices and width for the edges.
#'
#' @return A list of the following: \cr
#' subnet_label: label of the subnet; num_edges: number of edges in the sub-network; \cr
#' vertex_idx_selected, edge_pair_selected: the sub-network selection arguments provided by you; \cr
#' t_stat: t statistics of edges in the sub-network, t_stat_perc: percentiles of these t statistics among all edges in the original network; \cr
#' AT1: Value of AT1; pval_AT1: p-value of AT1; \cr
#' AT2: Value of AT2; pval_AT2: p-value of AT2)
#'
#' @export
#'
#' @importFrom here here
#' @importFrom foreach foreach "%dopar%"
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom grDevices dev.off pdf
#' @importFrom stats ecdf pnorm sd t.test
#' @import igraph
#'
#'
#' @examples
#' library(ExprNet)
#' network <- network_demo
#' data_type1 <- data_LGG_demo
#' data_type2 <- data_GBM_demo
#' edge_pair_selected <- c("1-8", "1-15", "2-16", "3-16", "5-10",
#'                         "5-16", "8-11", "8-13", "8-14", "8-15", "13-15")
#' save_dir <- here::here("demo")
#'
#' # compute the t-statistics and percentiles
#' res <- compute_edge_t_stat(data_type1, data_type2, network,
#'                            type1_name = "LGG", type2_name = "GBM", save_dir = save_dir)
#' edge_t_stat <- res$edge_t_stat
#' edge_dist_mat <- res$edge_dist_mat
#'
#' # compute the AT's
#' compute_AT(edge_t_stat = edge_t_stat, edge_dist_mat = edge_dist_mat, network = network,
#'            edge_pair_selected = edge_pair_selected,
#'            AT2_perm_test = TRUE, num_perm = 250, num_cores = 2,
#'            subnet_label = "Demo_GO0006306_DNA_methylation(LGG-GBM)",
#'            vertex.label.cex = 1, vertex.size = 10, edge.width = 7)

compute_AT <- function(edge_t_stat, edge_dist_mat, network,
                       vertex_idx_selected = NULL, edge_pair_selected = NULL,
                       AT2_perm_test = TRUE, num_perm = 500, num_cores = parallel::detectCores(),
                       type1_name = "Type1", type2_name = "Type2",
                       subnet_label = "subnet", plot_subnet = TRUE, save_plot = FALSE, save_dir = here::here(), ...) {

  ## This function computes AT1, AT2 and its p-value and optionally plots the sub-network and saves it. ##
  network <- as.undirected(network, mode = "collapse")

  if (!is.null(vertex_idx_selected)) {
    # select edges b/w ONLY selected genes
    row_selector <- ((edge_t_stat$V1 %in% vertex_idx_selected) & (edge_t_stat$V2 %in% vertex_idx_selected))
  }
  if (!is.null(edge_pair_selected)) {
    # select specified edges only
    if (!is.null(vertex_idx_selected)) warning("The sub-network is specified by edge_pair_selected (overriding vertex_idx_selected)")
    avail_edges <- apply(cbind(edge_t_stat$V1, edge_t_stat$V2), 1, paste, collapse = "-")
    row_selector <- (avail_edges %in% edge_pair_selected)
  }
  if (is.null(vertex_idx_selected) & is.null(edge_pair_selected)) {
    stop("Please specify either vertex_idx_selected or edge_pair_selected.")
  }

  if (sum(row_selector)==0) {  # no selected edge in the network
    warning("No selected edge in the network. Please check you input.")
    return (NULL)
  }

  edge_t_stat_selected <- edge_t_stat[row_selector, c("V1", "V2", "t_stat", "t_stat_perc")]  # select edges b/w ONLY selected genes

  num_edges <- nrow(edge_t_stat_selected)
  cat("\nThe subnetwork consists of ", num_edges, " edges.\n")


  ### Compute AT1 ###
  cat("\n---Computing AT1---\n")
  AT1 <- ifelse(num_edges >= 1, 1 - mean(edge_t_stat_selected$t_stat_perc), NA)
  # if (is.na(AT1) | (num_edges < 4)) return (NULL)  # not sufficient edges
  if (is.na(AT1)) {
    warnings("AT1 is NA. Please check selection of edges.\n")
    return (NULL)
  }

  # Compute p-value for AT1
  if (num_edges <= 10) {
    cdf <- pirwin.hall(num_edges*(1-AT1), num_edges)  # if num_edges is small, use Irwin-Hall(num_edges)
  } else {
    cdf <- pnorm(num_edges*(1-AT1), mean = num_edges/2, sd = sqrt(num_edges/12))  # if num_edges is large, approximated by N(num_edges/2, num_edges/12)
  }
  pval_AT1 <- 2*min(cdf, 1-cdf)
  cat("AT1 = ", round(AT1, 2), ", pval_AT1 = ", round(pval_AT1, 4), "\n")

  ### Compute AT2 ###
  cat("\n---Computing AT2---\n")
  AT2 <- ifelse(num_edges >= 1, mean(abs(edge_t_stat_selected$t_stat_perc - 0.5)), NA)
  if (is.na(AT2)) {
    warnings("AT2 is NA. Please check selection of edges.\n")
    return (NULL)
  }

  # Compute p-value for AT2 (Permutation test)
  pval_AT2 <- NA
  if (AT2_perm_test == TRUE) {  # permutation test for AT2
    cat("The permutation test may take some time, especially for high dimension. Please stay tuned.\n")

    num_cores <- ifelse(num_cores <= parallel::detectCores(), num_cores, parallel::detectCores())
    doParallel::registerDoParallel(num_cores)  # parallel computing
    cat("Parallel Computing: ", num_cores, " cores registered.\n")

    '%dopar%' <- foreach::'%dopar%'
    perm_idx <- NULL  # accommodate to notes in R CHECK
    AT2_perm <- foreach::foreach(perm_idx = 1:num_perm, .combine = 'c') %dopar% {
      AT2_one_perm(perm_idx, edge_t_stat_selected, edge_dist_mat)
    }
    doParallel::stopImplicitCluster()

    pval_AT2 <- min(mean(AT2_perm <= AT2, na.rm = TRUE), mean(AT2_perm >= AT2, na.rm = TRUE))
  }
  cat("AT2 computed.\nAT2 = ", round(AT2, 2), ", pval_AT2 = ", round(pval_AT2, 4), "\n")


  # Plot the subnetwork
  if (plot_subnet == TRUE) {
    # select genes and edges
    pos_edges_idx <- edge_t_stat_selected[edge_t_stat_selected$t_stat > 0, c("V1", "V2")]  # vertex indices for edges with POSITIVE t-values
    neg_edges_idx <- edge_t_stat_selected[edge_t_stat_selected$t_stat < 0, c("V1", "V2")]  # vertex indices for edges with NEGATIVE t-values

    # transform to igraph edge objects
    P_pos <- P_neg <- c()
    if (nrow(pos_edges_idx) >= 1) {  # if no pos edges selected, then skip it (otherwise as.character(df null row) will cause problems)
      for (i in 1:nrow(pos_edges_idx)) P_pos <- c(P_pos, as.character(pos_edges_idx[i,]))  # read the document of E to understand this
    }
    if (nrow(neg_edges_idx) >= 1) {
      for (i in 1:nrow(neg_edges_idx)) P_neg <- c(P_neg, as.character(neg_edges_idx[i,]))
    }

    selected_edges <- E(network, c(P_pos, P_neg))  # all edges involved in the pathway
    pos_edges <- E(network, P_pos)
    neg_edges <- E(network, P_neg)

    # Plot edges in different colors (edges with positive t-values: red; edges with negative t-values: blue)
    E(network)$color <- "grey"
    E(network)$color[E(network) %in% pos_edges] <- "red"
    E(network)$color[E(network) %in% neg_edges] <- "blue"

    subgraph <- subgraph.edges(network, selected_edges, delete.vertices = TRUE)  # the subgraph related to the pathway

    plot(subgraph, main = paste0(subnet_label),
         sub = paste0("AT1 = ", round(AT1, 2), ", pval_AT1 = ", round(pval_AT1, 4), "\n",
                      "AT2 = ", round(AT2, 2), ", pval_AT2 = ", round(pval_AT2, 4), "\n",
                      " (t statistic +: red, -: blue)"),
         ...,                                                        # user-given plot parameters, can override the below default parameters
         vertex.label.cex = 0.5, vertex.size = 3, edge.width = 0.5)  # default plot parameters
  }

  # Save the plot
  if (save_plot == TRUE) {
    cat("Please find results in the subfolder result_ExprNet")
    if (dir.exists(save_dir) == FALSE) { # save_dir is invalid, prompt warning and save to default directory
      save_dir <- here::here()  # by default
      warning(paste0("Invalid save directory.\nResults are saved under this directory: ", save_dir))
    }

    path <- paste0(save_dir, "/result_ExprNet/")
    dir.create(path) # create a subfolder result_ExprNet to save results
    pdf(paste0(path, "AT1_", subnet_label, "-", type1_name, "_", type2_name, ".pdf"))

    plot(subgraph, main = paste0(subnet_label),
         sub = paste0("AT1 = ", round(AT1, 2), ", pval_AT1 = ", round(pval_AT1, 4), "\n",
                      "AT2 = ", round(AT2, 2), ", pval_AT2 = ", round(pval_AT2, 4), "\n",
                      " (t statistic +: red, -: blue)"),
         ...,                                                        # user-given plot parameters, can override the below default parameters
         vertex.label.cex = 0.5, vertex.size = 3, edge.width = 0.5)  # default plot parameters

    dev.off()
    cat("---Plot saved---\n")
  }

  return (list(subnet_label = subnet_label, num_edges = num_edges,
               vertex_idx_selected = vertex_idx_selected, edge_pair_selected = edge_pair_selected,
               t_stat = edge_t_stat_selected$t_stat, t_stat_perc = edge_t_stat_selected$t_stat_perc,
               AT1 = AT1, pval_AT1 = pval_AT1,
               AT2 = AT2, pval_AT2 = pval_AT2))
}
