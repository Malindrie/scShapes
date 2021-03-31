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
#' @param workers Number of workers to be used in parallel computation
#' using \code{future.apply}, with argument \code{multisession}.
#'
#' @param seed Seed number to be used in parallel computation
#' using \code{future.apply}, with argument \code{multisession}
#'
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom stats ecdf
#' @importFrom parallelly availableCores
#' @importFrom Matrix Matrix
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom pscl zeroinfl
#' @importFrom VGAM predict
#' @importFrom VGAM rzinegbin
#' @importFrom dgof ks.test
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
#' scData_KS <- ks_test(counts=scData$counts, cexpr=scData$covariates, lib.size=scData$lib_size)


ks_test <- function(counts, cexpr, lib.size,
                    formula=NULL, workers=NULL,
                    seed=NULL){

  if(is.null(workers)) {
    workers <- min(4, parallelly::availableCores())
  }
  if(is.null(seed)) {
    seed <- 0xBEEF
  }

  #Formulate a simple additive model using all the covariates in 'cexpr'
  covariates <- names(cexpr)
  if(is.null(formula)) {
    message(sprintf("Formulating the additive model..."))
    formula <- 'x ~ 1 '
    if(!identical(covariates, character(0))) {
      for (covar in covariates) {
        formula <- paste(formula, sprintf(' + %s', covar))
      }
    }
  }
  formula <- as.formula(formula)

  #set-up multisession
  plan(future::multisession, workers = workers)

  #convert data to lists
  gexpr <- apply(counts, 1, function (x) cbind(x,cexpr))

  #KS test with simulated p-values
  message(sprintf("Performing the KS test..."))
  KS_ZINB <- function(data, formula, lib.size){

    library(pscl)
    m1 <- try(zeroinfl(formula, data, offset=log(lib.size), dist = "negbin"), silent = TRUE)

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
      pp <- try(rzinegbin(n = length(data$x),
                          size = ccc[2,],
                          munb = ccc[3,],
                          pstr0 = ccc[1,]),
                          silent = TRUE)
    }
    else {
      pp <- "NA"
    }

    if(!(class(pp) == "character")){

      library(dgof)
      D <- try(ks.test(data$x, ecdf(pp), simulate.p.value = TRUE)$p.value, silent = TRUE)
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

  ks.pval.unadj <- future.apply::future_lapply(gexpr,
                                               FUN = KS_ZINB,
                                               formula <- formula,
                                               lib.size <- lib.size,
                                               future.seed=TRUE)
  return(ks.pval.unadj)

}
