#' ks_test
#'
#' this function is used to perform Kolmogorv-Smirnov test on the
#' filtered sparse counts matrix from \code{filter_counts} to select genes
#' belonging to the family of ZINB distributions
#'
#' @param counts A non-negative integer matrix of scRNA-seq filtered read counts.
#' The rows of the matrix are genes and columns are samples/cells.
#'
#' @param cexpr A dataframe that contains the covariate values.
#' The rows of the dataframe are the corresponding samples/cells from the counts
#' matrix from \code{filter_counts}.
#' The cells of the dataframe are the covariates to be included in the GLM.
#'
#' @param formula A regression formula to fit the covariates in the ZINB GLM.
#'
#' @param workers Number of workers to be used in parallel computation
#' using \code{future.apply}, with argument \code{multisession}.
#'
#' @param seed Seed number to be used in parallel computation
#' using \code{future.apply}, with argument \code{multisession}
#'
#' @export
#'
#' @import dplyr
#' @importFrom Matrix Matrix
#' @importFrom future.apply future_apply
#' @importFrom pscl zeroinfl
#' @importFrom VGAM rzinegbin
#' @importFrom dgof ks.test
#'
#' @return @return List object containing the significant gene indices from the KS test,
#' their adjusted p-values


ks_test <- function(counts, cexpr, formula=NULL, workers=NULL, seed=NULL){

  if(is.null(workers)) {
    workers <- min(4, parallel::availableCores())
  }
  if(is.null(seed)) {
    seed <- 0xBEEF
  }

  #Formulate a simple additive model using all the covariates in 'cexpr'
  covariates <- names(cexpr)
  if(is.null(formula_string)) {
    message(sprintf("Formulating the default additive model..."))
    formula <- 'x ~ 1 |'
    if(!identical(covariates, character(0))) {
      for (covar in covariates) {
        formula <- paste(formula_string, sprintf(' + %s', covar))
      }
    }
  }

  #set-up multisession
  plan(multisession, workers = workers)

  #Calculate library size for each cell
  lib.size <- apply(counts,2, function(x) sum(x))

  #KS test with simulated p-values
  KS_ZINB <- function(x, lib.size){

    library(pscl)
    m1 <- try(zeroinfl(formula, offset=log(lib.size), dist = "negbin"), silent = TRUE)

    if(!(class(m1) == "try-error")){
      pi_ML = predict(m1, type = "zero")
      theta_ML = m1$theta
      mean_ML = predict(m1, type = "count")
      var_ML = mean_ML + (mean_ML ^ 2 / theta_ML)
      ccc = rbind(pi_ML, theta_ML, mean_ML)
    }
    else {
      ccc = "NA"
    }


    if(!(class(ccc) == "character")){
      library(VGAM)
      pp <- try(rzinegbin(n = length(x), size = ccc[2,], munb = ccc[3,], pstr0 = ccc[1,]), silent = TRUE)
    }
    else {
      pp <- "NA"
    }

    if(!(class(pp) == "character")){

      library(dgof)
      D <- try(ks.test(x, ecdf(pp), simulate.p.value = TRUE)$p.value, silent = TRUE)
    }
    else {
      D <- "NA"
    }

    if(!(class(D) == "character")){

      p_value <- D
    }
    else {
      p_value <- "NA"
    }

    return(p_value)
  }

  ks.pval.unadj <- future_apply(counts, MARGIN = 1L, FUN = KS_ZINB, lib.size <- lib.size, future.seed = seed)

  #Remove genes that filed the KS test
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
