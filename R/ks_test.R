#' ks_test
#'
#' This function is used to perform Kolmogorv-Smirnov test on the
#' filtered sparse counts matrix from \code{filter_counts} to select genes
#' belonging to the family of ZINB distributions
#'
#' @param counts A non-negative integer matrix of scRNA-seq filtered read counts
#' from \code{filter_counts}.
#' The rows of the matrix are genes and columns are samples/cells.
#'
#' @param cexpr A dataframe that contains the covariate values.
#' The rows of the dataframe are the corresponding samples/cells from the counts
#' matrix from \code{filter_counts}.
#' The cells of the dataframe are the covariates to be included in the GLM.
#'
#' @param lib.size A numeric vector that contains the total number of counts
#' per cell from the counts matrix from \code{filter_counts}.
#'
#' @param formula A regression formula to fit the covariates in the ZINB GLM.
#'
#' @param BPPARAM configuration parameter related to the method of parallel execution.
#' For further information on how to set-up parallel execution refer to
#' \code{BiocParallel} vignette.
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom stats as.formula
#' @importFrom stats ecdf
#' @importFrom Matrix Matrix
#' @importFrom pscl zeroinfl
#' @importFrom VGAM predict
#' @importFrom VGAM rzinegbin
#' @importFrom dgof ks.test
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#'
#' @return List object containing the p-values from the
#' KS test.
#'
#' @examples
#'
#' #' # load toy  example data
#'
#' data(scData)
#'
#' # apply the ks_test function to subset genes belonging to the
#' # family of ZINB distributions.
#'
#' scData_KS <- ks_test(counts=scData$counts, cexpr=scData$covariates, lib.size=scData$lib_size, BPPARAM=bpparam())



ks_test <- function(counts, cexpr, lib.size,
                    formula=NULL, BPPARAM){

  set.seed(0xBEEF)

  #Formulate a simple additive model using all the covariates in 'cexpr'
  covariates <- names(cexpr)
  if(is.null(formula)) {
    message(sprintf("Formulating the additive model..."))
    formula <- 'x ~ 1 '
    if(!identical(covariates, NULL)) {
      for (covar in covariates) {
        formula <- paste(formula, sprintf(' + %s', covar))
      }
    }
  }
  formula <- as.formula(formula)

  #convert data to lists
  if(!identical(covariates, NULL)) {
    gexpr <- apply(counts, 1, function (x) cbind(x,cexpr))
  } else{
    gexpr <- apply(counts, 1, function (x) as.list(as.data.frame(x)))
  }

  #KS test with simulated p-values
  message(sprintf("Performing the KS test..."))
  KS_ZINB <- function(data, formula, lib.size){

    m1 <- tryCatch(pscl::zeroinfl(formula, data, offset=log(lib.size), dist = "negbin"), error = function(e) {print(e$message)})

    if(!is(m1, "character")){
      pi_ML = predict(m1, type = "zero")
      theta_ML = m1$theta
      mean_ML = predict(m1, type = "count")
      var_ML = mean_ML + (mean_ML ^ 2 / theta_ML)
      ccc = rbind(pi_ML, theta_ML, mean_ML)
    }
    else {
      ccc = "NA"
    }


    if(!is(ccc, "character")){

      pp <- tryCatch(VGAM::rzinegbin(n = length(data$x),
                          size = ccc[2,],
                          munb = ccc[3,],
                          pstr0 = ccc[1,]),
                     error = function(e) {print(e$message)})
    }
    else {
      pp <- "NA"
    }

    if(!is(pp, "character")){

      D <- tryCatch(dgof::ks.test(data$x, ecdf(pp), simulate.p.value = TRUE)$p.value, error = function(e) {print(e$message)})
    }
    else {
      D <- "NA"
    }

    if(!is(D, "character")){

      p_value <- D
    }
    else {
      p_value <- "NA"
    }

    return(p_value)
  }


  if(identical(covariates, NULL)) {
    ks.pval.unadj <- BiocParallel::bplapply(gexpr,
                                            FUN = KS_ZINB,
                                            BPPARAM = BPPARAM,
                                            formula <- formula,
                                            lib.size <- lib.size)
  } else {
    ks.pval.unadj <- BiocParallel::bplapply(gexpr,
                                            FUN = KS_ZINB,
                                            BPPARAM = BPPARAM,
                                            formula <- formula,
                                            lib.size <- lib.size)
  }
  return(ks.pval.unadj)

}
