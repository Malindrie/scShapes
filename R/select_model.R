#' select_model
#'
#' This function is used to select the distribution of best fit for scRNA-seq count data
#'
#' @param lrt.value A list of genes with the p-values from performing the GOF tests
#' from \code{gof_model}
#'
#' @export
#'
#' @importFrom stats p.adjust
#'
#' @return A list of selected model distributions for genes scShapes selects.
#'
#' @examples
#'
#' data(scData)
#'
#' # apply the select_model function to the best fit model from the results of
#' # the gof_model function
#'
#' scData_models <- fit_models(counts=scData$counts, cexpr=scData$covariates,
#'                             lib.size=scData$lib_size, BPPARAM=bpparam())
#' scData_bicvals <- model_bic(scData_models)
#' scData_least.bic <- lbic_model(scData_bicvals, scData$counts)
#' scData_gof <- gof_model(scData_least.bic, cexpr=scData$covariates, lib.size=scData$lib_size,
#' BPPARAM=bpparam())
#'
#' scData_fit <- select_model(scData_gof)


select_model <- function(lrt.value){

  #p-value adjustment using BH correction
  pval_lrt <- lapply(lrt.value, function(x) p.adjust(x, method="BH"))

  #Check for significance of p-values
  P_genes <- unlist(lapply(unlist(pval_lrt[["P_lrt"]]), function (x) x[x > 0.05]))
  NB_genes <- unlist(lapply(pval_lrt[["NB_lrt"]], function (x) x[x > 0.05]))
  ZIP_genes <- unlist(lapply(pval_lrt[["ZIP_lrt"]], function (x) x[x < 0.05]))
  ZINB_genes <- unlist(lapply(pval_lrt[["ZINB_lrt"]], function (x) x[x < 0.05]))


  #grep the names of genes following each distribution

  gene_dist <- list()

  gene_dist[["P_genes"]] <- names(P_genes)
  gene_dist[["NB_genes"]] <- names(NB_genes)
  gene_dist[["ZIP_genes"]] <- names(ZIP_genes)
  gene_dist[["ZINB_genes"]] <- names(ZINB_genes)


  return(gene_dist)
}
