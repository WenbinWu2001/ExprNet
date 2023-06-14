#' LGG Demo Data
#'
#' An LGG normalized gene expression matrix with selected genes.
#'
#' @format ## `data_LGG_demo`
#' A data frame with 17 rows (genes) and 6 columns (the first column for vertex indices + 5 samples):
#' \describe{
#' The demo dataset gives the expression levels of genes which are involved in GO:0006306 (DNA methylation) and GO:0051781 (positive regulation of cell division) and have at least one neighbor on \link[demo network]{network_demo}.
#' Each row corresponds to a gene. Each column corresponds to a sample. The first column are the numeric vertex indices. \cr
#' The dataset is the extracted LGG part of an aggregated dataset of LGG, GBM, LUSC and LUAD datasets on TCGA after quantile normalization.
#' }
#' @source <https://www.who.int/teams/global-tuberculosis-programme/data>
"data_LGG_demo"
