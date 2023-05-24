#' Distribution function of the Irwin-Hall distribution
#' @description Computes the distribution function of Irwin-Hall with parameter n.
#' @param x A number between 0 and n.
#' @param n The number of i.i.d. U(0,1) to sum.
#'
#' @return The CDF of the Irwin-Hall(n) at x.
#' @export
#' @references
#' Irwin-Hall Distribution. (n.d.). Randomservices.org. Retrieved May 8, 2023, from https://www.randomservices.org/random/special/IrwinHall.html
#' @examples
#' pirwin.hall(1, 2)
pirwin.hall <- function(x, n) {
  k = 0:n
  return ( 1/2 + 1/(2*factorial(n)) * sum( (-1)^k * choose(n, k) * sign(x-k) * (x-k)^n ) )
}


col_t_test_fun <- function(col_data, label) {
  # conducts a two sample t-test on a data vector [col_data] with assigned label (2 factors). Return the t-statistic.
  t.test(formula = col_data ~ as.factor(label),
         alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE)$statistic
}


AT2_one_perm <- function(perm_idx, edge_t_stat_selected, edge_dist_mat) {
  # Compute AT2 with one random permutation of the label.
  perm_label <- sample(edge_dist_mat[,1])  # recall that the first column is sample label

  t_stat_perm <- apply(edge_dist_mat[,-1], 2, col_t_test_fun, label = perm_label)  # t-statistic for ALL edges with permuted label
  t_perc_perm <- ecdf(t_stat_perm)(t_stat_perm)  # percentile of each t-statistic among ALL t-statistic

  names(t_stat_perm) <- names(t_perc_perm) <- colnames(edge_dist_mat[,-1])
  t_perc_perm_selected <- t_perc_perm[apply(cbind(edge_t_stat_selected$V1, edge_t_stat_selected$V2), 1, paste, collapse = "-")]

  AT2 <- mean(abs(t_perc_perm_selected - 0.5))  # Compute AT2 associated with the permuted label

  return (AT2)
}


