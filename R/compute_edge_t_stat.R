#' Computing edge length, t statistics and t percentiles
#'
#' @description Given data for two phenotypes and a network, compute the length of each edge and conduct two-sample t tests on all edges.
#'
#' @param data_type1 data.frame for phenotype 1. Each row corresponds to a vertex feature.
#' Each column corresponds to a sample. The first column are the numeric vertex indices.
#' @param data_type2 data.frame for phenotype 2, in the same format as *data_type1*.
#' @param network An \link[igraph]{igraph-package} graph object. It will be converted to an undirected graph by default.
#' @param type1_name The name for phenotype 1, used for naming files of the results.
#' @param type2_name The name for phenotype 2.
#' @param alpha_t_test Significance level of the t test for each edge. It does not affect the computation of AT1 and AT2.
#' @param save_edge_res Logical. Whether the results should be saved. Set to TRUE if you plan to conduct analysis on the t statistics of the edges.
#' @param save_dir Directory to save the results. Results are saved in a subfolder *result_ExprNet* will be created under the directory.
#'
#' @return A list consisting of *edge_t_stat* and *edge_dist_mat*.
#' *edge_t_stat* is a data.frame storing results of the t test.
#' Each row gives the following properties of an edge: Vertex1(V1), Vertex2(V2), (mean, sd) of each type, difference in means, t statistic, p-value, percentile of the t statistic.
#' *edge_dist_mat* is a matrix containing all the sample distances of each edge.
#' Each column is the computed edge length of an edge *"vertex1-vertex2"* for all samples (Type1 and Type2). The first column is the sample label.
#' @export
#'
#' @examples
#' library(ExprNet)
#' network <- network_demo
#' data_type1 <- data_LGG_demo
#' data_type2 <- data_GBM_demo
#'
#' # compute the t-statistics and percentiles
#' res <- compute_edge_t_stat(data_type1, data_type2, network, type1_name = "LGG", type2_name = "GBM")
#' res

compute_edge_t_stat <- function(data_type1, data_type2, network, type1_name = "Type1", type2_name = "Type2",
                           alpha_t_test = 0.05, save_edge_res = FALSE, save_dir = here::here()) {
  ## This function computes edge length (distance between any pair of connected vertices) and t-stat and optionally saves the results. ##

  colnames(data_type1)[1] <- "vertex_idx"  # the first column should be vertex index (numeric)
  colnames(data_type2)[1] <- "vertex_idx"
  network <- as.undirected(network, mode = "collapse")
  edge_list <- get.edgelist(network)   # all edges on the network

  num_type1 <- dim(data_type1[,-1])[2]  # number of type 1 samples (The first column is vertex index)
  num_type2 <- dim(data_type2[,-1])[2]  # number of type 2 samples (The first column is vertex index)
  message(paste0("\nData read successfully.\n",
             dim(data_type1[,-1])[1], " features\n",
             num_type1, " samples for phenotype 1\n",
             num_type2, " samples for phenotype 2."))
  sample_label <- c(rep(1, num_type1), rep(2, num_type2))  # the true sample label (1 for type1, 2 for type2), corresponding to rows in edge_dist_mat

  message(paste0("\nGraph imported successfully.\n",
             "There are ", length(V(network)), " vertices and ", nrow(edge_list), " edges in the graph."))

  ## store results as a (#edge x 10) dataframe:
  ## each row corresponds to results of an edge:
  ## Vertex1(V1), Vertex2(V2), (mean, sd) of each type, difference in means, t statistic, p-value, percentile of the t statistic
  ## Remark: distance between nodes Vi and Vj is calculated as the expression on Vi - expression on Vj (the direction actually doesn't matter)
  edge_t_stat <- as.data.frame(matrix(NA, nrow = nrow(edge_list), ncol = 10))
  colnames(edge_t_stat) <- c("V1", "V2",
                               paste0("mean_", type1_name), paste0("mean_", type2_name),
                               paste0("sd_", type1_name), paste0("sd_", type2_name),
                               "diff_in_means", "t_stat", "p_value", "t_stat_perc")


  message("\n---Computing edge length and t statistics---")
  ## Store the results in a (#samples in type1 + #samples in type2) x #edges matrix.
  ## Each row gives distances (edge length) for a sample. Each col corresponds to the edge "vertex1-vertex2" of which the distance is computed.
  ## The first col is sample label.
  edge_dist_mat <- matrix(NA, nrow = num_type1 + num_type2, ncol = dim(edge_list)[1])
  colnames(edge_dist_mat) <- apply(edge_list, 1, paste, collapse = "-")

  for (i in (1:nrow(edge_list))) {
    vertex1 <- edge_list[i, 1]  # numeric (vertex index)
    vertex2 <- edge_list[i, 2]
    vertex_pair <- paste0(vertex1, "-", vertex2)

    samples_type1 <-
      as.numeric(data_type1[data_type1$vertex_idx == vertex1, -1]) - as.numeric(data_type1[data_type1$vertex_idx == vertex2, -1])
    samples_type2 <-
      as.numeric(data_type2[data_type2$vertex_idx == vertex1, -1]) - as.numeric(data_type2[data_type2$vertex_idx == vertex2, -1])

    edge_dist_mat[sample_label == 1, vertex_pair] <- samples_type1
    edge_dist_mat[sample_label == 2, vertex_pair] <- samples_type2

    test_res <- t.test(x = samples_type1, y = samples_type2,
                       alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE,
                       conf.level = 1-alpha_t_test)  # two-sided two sample t test

    edge_t_stat[i, ] <- c(vertex1, vertex2,
                            mean(samples_type1), mean(samples_type2), sd(samples_type1), sd(samples_type2),
                            mean(samples_type1) - mean(samples_type2), test_res$statistic, test_res$p.value,
                            "NA")  # save results

    if ((i %in% as.integer(nrow(edge_list)/5 * 1:5))) message(paste0(i, " / ", nrow(edge_list), " Edges Computed"))
  }

  edge_t_stat$t_stat_perc <- ecdf(edge_t_stat$t_stat)(edge_t_stat$t_stat)
  edge_dist_mat <- cbind(sample_label, edge_dist_mat)  # Add a leftmost col as sample labels

  message(paste0("\nAmong ", nrow(edge_list), " edge distances, \n",
               sum(edge_t_stat$p_value < alpha_t_test), " of them have significant differences at ", alpha_t_test, " level.\n",
               sum(edge_t_stat$diff_in_means > 0), " of them > 0.\n",
               sum(edge_t_stat$diff_in_means < 0), " of them < 0.\n"))

  # save results
  if (save_edge_res == TRUE) {
    message("\n---Saving results---")
    message("Please find results in the subfolder result_ExprNet\n")
    if (dir.exists(save_dir) == FALSE) { # save_dir is invalid, prompt warning and save to default directory
      save_dir <- here::here()  # by default
      warning(paste0("Invalid save directory.\nResults are saved under this directory: ", save_dir))
    }

    path <- file.path(save_dir, "result_ExprNet")
    dir.create(path)  # create a subfolder result_ExprNet to save results
    save(edge_dist_mat, file = file.path(path, paste0(type1_name, "_", type2_name, "_edge_distances.Rdata")))
    save(edge_t_stat, file = file.path(path, paste0(type1_name, "_", type2_name, "_edge_t_stat.Rdata")))
  }

  return (list(edge_t_stat = edge_t_stat, edge_dist_mat = edge_dist_mat))
}
