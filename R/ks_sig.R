#' ks_sig
#'
#' This function is used to select genes significant from the
#' \code{ks_test}.
#'
#' @param ks.pval.unadj The output from \code{ks_test} which is a list
#' of p-values from the KS test with gene names.
#'
#' @export
#'
#' @importFrom stats p.adjust
#' @importFrom utils stack
#'
#' @return List object containing the significant gene indices from the KS test,
#' their adjusted p-values
#'
#' @examples
#'
#' data(scData)
#'
#' # apply the ks_test function to subset genes belonging to the
#' # family of ZINB distributions.
#'
#' scData_KS <- ks_test(counts=scData$counts, cexpr=scData$covariates, lib.size=scData$lib.size)
#'
#' # apply the ks_sig function to select genes significant from
#' # the Kolmogorov Smirnov test.
#'
#' scData_KS_sig <- ks_sig(scData_KS)


ks_sig <- function(ks.pval.unadj){


  #Remove genes that failed the KS test
  ks.pval.unadj <- stack(ks.pval.unadj)
  ks.pval.unadj$values <- as.numeric(ks.pval.unadj$values)
  ks.pval.unadj <- ks.pval.unadj[!is.na(ks.pval.unadj$values), ]

  ks.pval.unadj_vec <- as.numeric(ks.pval.unadj$values)
  names(ks.pval.unadj_vec) <- ks.pval.unadj$ind

  ks.pval <- p.adjust(ks.pval.unadj_vec, method="BH")

  #Select genes that pass the KS test
  sig_genes_ks <- which(ks.pval > 0.01)

  return(list(genes=sig_genes_ks, p=ks.pval, p.unadj=ks.pval.unadj))

}


