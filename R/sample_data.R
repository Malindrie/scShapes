#' @title Sample data for analysis
#'
#' @description Toy example data list of scRNA-seq counts,
#' information on covariates, and library sizes for randomly
#' generated 20 genes to illustrate how to use the functions
#' of the package \code{scShapes}
#'
#' @name scData
#'
#' @docType data
#'
#' @usage data(scData)
#'
#' @format A list of three lists labeled 'counts', 'covariates',
#' 'lib_size'. 'counts' an RNA-seq counts matrix of 20 genes
#' and 500 cells; 'covariates' a dataframe of covariates corresponding
#' to the RNA-seq counts where rows are cells of the counts matrix;
#' 'lib_size' a numeric vector of library sizes corresponding to the
#' columns of the RNA-seq counts matrix.
#'
#' @return An RData object
#'
#' @examples
#'
#' # load toy  example data
#'
#' data(scData)
NULL
