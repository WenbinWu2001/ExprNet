#' LGG Demo Data
#'
#' An LGG normalized gene expression matrix with selected genes.
#'
#' @format ## `data_LGG_demo`
#' A data frame with 17 rows (genes) and 6 columns (the first column for vertex indices + 5 samples).
#' @details
#' The demo dataset gives the expression levels of genes which are involved in GO:0006306 (DNA methylation) and GO:0051781 (positive regulation of cell division) and have at least one neighbor on [network_demo].
#' Each row corresponds to a gene. Each column corresponds to a sample. The first column are the numeric vertex indices.\cr\cr
#' The dataset is the extracted LGG part of an aggregated dataset of LGG, GBM, LUSC and LUAD datasets on TCGA after quantile normalization.\cr
#' LGG and GBM are two groups of brain tumors and are usually studied comparatively (as in this demo study).
#' @source https://xenabrowser.net/datapages/?dataset=TCGA-LGG.htseq_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
"data_LGG_demo"


#' GBM Demo Data
#'
#' A GBM normalized gene expression matrix with selected genes.
#'
#' @format ## `data_GBM_demo`
#' A data frame with 17 rows (genes) and 6 columns (the first column for vertex indices + 5 samples).
#' @details
#' The demo dataset gives the expression levels of genes which are involved in GO:0006306 (DNA methylation) and GO:0051781 (positive regulation of cell division) and have at least one neighbor on [network_demo].
#' Each row corresponds to a gene. Each column corresponds to a sample. The first column are the numeric vertex indices.\cr\cr
#' The dataset is the extracted GBM part of an aggregated dataset of LGG, GBM, LUSC and LUAD datasets on TCGA after quantile normalization.\cr
#' LGG and GBM are two groups of brain tumors and are usually studied comparatively (as in this demo study).
#' @source https://xenabrowser.net/datapages/?dataset=TCGA-GBM.htseq_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
"data_GBM_demo"


#' LUSC Demo Data
#'
#' An LUSC normalized gene expression matrix with selected genes.
#'
#' @format ## `data_LUSC_demo`
#' A data frame with 17 rows (genes) and 6 columns (the first column for vertex indices + 5 samples).
#' @details
#' The demo dataset gives the expression levels of genes which are involved in GO:0006306 (DNA methylation) and GO:0051781 (positive regulation of cell division) and have at least one neighbor on [network_demo].
#' Each row corresponds to a gene. Each column corresponds to a sample. The first column are the numeric vertex indices.\cr\cr
#' The dataset is the extracted LUSC part of an aggregated dataset of LGG, GBM, LUSC and LUAD datasets on TCGA after quantile normalization.\cr
#' LUSC and LUAD are two groups of lung tumors and are usually studied comparatively (as in this demo study).
#' @source https://xenabrowser.net/datapages/?dataset=TCGA-LUSC.htseq_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
"data_LUSC_demo"


#' LUAD Demo Data
#'
#' An LUAD normalized gene expression matrix with selected genes.
#'
#' @format ## `data_LUAD_demo`
#' A data frame with 17 rows (genes) and 6 columns (the first column for vertex indices + 5 samples).
#' @details
#' The demo dataset gives the expression levels of genes which are involved in GO:0006306 (DNA methylation) and GO:0051781 (positive regulation of cell division) and have at least one neighbor on [network_demo].
#' Each row corresponds to a gene. Each column corresponds to a sample. The first column are the numeric vertex indices.\cr\cr
#' The dataset is the extracted LUAD part of an aggregated dataset of LGG, GBM, LUSC and LUAD datasets on TCGA after quantile normalization.\cr
#' LUSC and LUAD are two groups of lung tumors and are usually studied comparatively (as in this demo study).
#' @source https://xenabrowser.net/datapages/?dataset=TCGA-LUAD.htseq_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
"data_LUAD_demo"


#' Demo Network
#'
#' A gene regulatory network extracted from a HINT (High-quality INTeractomes) binary compilation of protein-protein interactions.
#'
#' @format ## `network_demo`
#' An \link[igraph]{igraph-package} graph object.
#'
#' @source http://hint.yulab.org/download/HomoSapiens/binary/hq/
"network_demo"
