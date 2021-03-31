#' fit_models
#'
#' This function is used to fit genes with GLM
#'
#' @param counts A non-negative integer matrix of scRNA-seq filtered read counts
#' containing genes belonging to the family of ZINB distributions selected from
#' \code{ks_test}.
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
#' @param model A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#'
#' @export
#'
#' @importFrom stats glm
#' @importFrom stats as.formula
#' @importFrom parallelly availableCores
#' @importFrom Matrix Matrix
#' @importFrom future.apply future_apply
#' @importFrom pscl zeroinfl
#' @importFrom MASS glm.nb
#' @importFrom future plan
#'
#' @return A list of models fitted by 'glm'
#'
#' @examples
#'
#' data(scData)
#'
#' # apply the fit_models function to subset genes belonging to the
#' # family of ZINB distributions, selceted from ks_test function.
#'
#' scData_models <- fit_models(counts=scData$counts, cexpr=scData$covariates,
#' lib.size=scData$lib_size)
#'
#' ## Shut down parallel workers
#' future::plan("sequential")


fit_models <- function(counts, cexpr, lib.size,
                       formula=NULL, workers=NULL,
                       seed=NULL, model=NULL){

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

  f_oth <- as.formula(formula)
  f_nb <- as.formula(paste(formula, ' + offset(log(lib.size))'))

  #set-up multisession
  plan(future::multisession, workers = workers)


  #convert data to lists
  gexpr <- apply(counts, 1, function (x) cbind(x,cexpr))


  #Fit the four distributions
  model_poi <- function(data, formula, lib.size){
    library(stats)
    try(glm(formula, data, offset=log(lib.size), family = "poisson"), silent = TRUE)
  }


  model_nb <- function(data, formula, lib.size){
    library(MASS)
    try(glm.nb(formula, data), silent = TRUE)
  }


  model_zip <- function(data, formula, lib.size){
    library(pscl)
    try(zeroinfl(formula, data, offset=log(lib.size), dist = "poisson"), silent = TRUE)
  }


  model_zinb <- function(data, formula, lib.size){
    library(pscl)
    try(zeroinfl(formula, data, offset=log(lib.size), dist = "negbin"), silent = TRUE)
  }


  fitting <- list()

  if(is.null(model) || 1 %in% model) {
    message('Fitting data with Poisson model...')
    if(identical(covariates, character(0))) {
      fitting[["P"]] <-   future.apply::future_lapply(gexpr,
                                                      FUN = model_poi,
                                                      formula = f_oth,
                                                      lib.size = lib.size,
                                                      future.seed=TRUE)
    } else {
      fitting[["P"]] <- future.apply::future_lapply(gexpr,
                                                    FUN = model_poi,
                                                    formula = f_oth,
                                                    lib.size = lib.size,
                                                    future.seed=TRUE)
    }
  }

  if(is.null(model) || 2 %in% model) {
    message('Fitting data with Negative Binomial model...')
    if(identical(covariates, character(0))) {
      fitting[["NB"]] <-   future.apply::future_lapply(gexpr,
                                                       FUN = model_nb,
                                                       formula = f_nb,
                                                       lib.size = lib.size,
                                                       future.seed=TRUE)
    } else {
      fitting[["NB"]] <- future.apply::future_lapply(gexpr,
                                                     FUN = model_nb,
                                                     formula = f_nb,
                                                     lib.size = lib.size,
                                                     future.seed=TRUE)
    }
  }

  if(is.null(model) || 3 %in% model) {
    message('Fitting data with Zero Inflated Poisson model...')
    if(identical(covariates, character(0))) {
      fitting[["ZIP"]] <-   future.apply::future_lapply(gexpr,
                                                       FUN = model_zip,
                                                       formula = f_oth,
                                                       lib.size = lib.size,
                                                       future.seed=TRUE)
    } else {
      fitting[["ZIP"]] <- future.apply::future_lapply(gexpr,
                                                     FUN = model_zip,
                                                     formula = f_oth,
                                                     lib.size = lib.size,
                                                     future.seed=TRUE)
    }
  }

  if(is.null(model) || 4 %in% model) {
    message('Fitting data with Zero Inflated Negative Binomial model...')
    if(identical(covariates, character(0))) {
      fitting[["ZINB"]] <-   future.apply::future_lapply(gexpr,
                                                       FUN = model_zinb,
                                                       formula = f_oth,
                                                       lib.size = lib.size,
                                                       future.seed=TRUE)
    } else {
      fitting[["ZINB"]] <- future.apply::future_lapply(gexpr,
                                                     FUN = model_zinb,
                                                     formula = f_oth,
                                                     lib.size = lib.size,
                                                     future.seed=TRUE)
    }
  }

  return(fitting)


  future:::ClusterRegistry("stop")

}

