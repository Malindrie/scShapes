#' gof_model
#'
#' This function is used to perform the likelihood ratio test on the models chosen
#' based on the BIC values from \code{best_model} to check for model adequacy.
#'
#' @param lbic A list of genes together with filtered read counts based on the selected
#' distribution from \code{best_model}. Output from \code{best_model}.
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
#' @importFrom stats glm
#' @importFrom stats as.formula
#' @importFrom stats pchisq
#' @importFrom MASS glm.nb
#' @importFrom VGAM logLik
#' @importFrom pscl zeroinfl
#' @importFrom emdbook pchibarsq
#'
#' @return A list of genes with the p-values from performing the GOF tests.
#'
#' @examples
#'
#' data(scData)
#'
#' # apply the gof_model function to perform the likelihood ratio
#' # test on the models selected by using the lbic_model function
#'
#' scData_gof <- gof_model(scData_least.bic, cexpr=scData$covariates, lib.size=scData$lib.size)



gof_model <- function(lbic, cexpr, lib.size,
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

  f_oth <- as.formula(formula)
  f_nb <- as.formula(paste(formula, ' + offset(log(lib.size))'))

  #set-up multisession
  plan(future::multisession, workers = workers)


  #Prepare input data as lists to be inputted to the function
  lbic[["P"]] <- lapply(lbic[["P"]], function (x) cbind(x,cexpr))
  lbic[["NB"]] <- lapply(lbic[["NB"]], function (x) cbind(x,cexpr))
  lbic[["ZIP"]] <- lapply(lbic[["ZIP"]], function (x) cbind(x,cexpr))
  lbic[["ZINB"]] <- lapply(lbic[["ZINB"]], function (x) cbind(x,cexpr))



  #Fit the four distributions
  #Poisson distribution
  model_Chipoi <- function(data, formula, lib.size){

    library(stats)

    1-pchisq(summary(glm(formula, data, offset=log(lib.size), family = "poisson"))$deviance,
             df= summary(glm(formula, data, offset=log(lib.size), family = "poisson"))$df.residual)
  }


  #NB distribution
  model_Chinb <- function(data, formula, lib.size){

    library(MASS)

    1-pchisq(summary(glm.nb(formula, data))$deviance,
             df= summary(glm.nb(formula, data))$df.residual)
  }


  #ZIP distribution
  model_Chizip <- function(data, formula, lib.size){
    library(pscl)
    library(emdbook)

    1-pchibarsq(2*(logLik(zeroinfl(formula, data, offset=log(lib.size), dist = "poisson"))-logLik(glm(formula, data, offset=log(lib.size), family = "poisson"))),df=1,mix=0.5,lower.tail=TRUE, log.p = FALSE)

  }


  #ZINB distribution
  model_Chizinb <- function(data, formula, lib.size){
    library(pscl)
    library(emdbook)

    1-pchibarsq(2*(logLik(zeroinfl(formula, data, offset=log(lib.size), dist = "negbin"))-logLik(glm.nb(formula, data))),df=1,mix=0.5,lower.tail=TRUE, log.p = FALSE)
  }


  fitting_lrt <- list()

  message('Performing LRT for Poisson model...')
  fitting_lrt[["P_lrt"]] <- future.apply::future_lapply(lbic[["P"]],
                                                  FUN = model_Chipoi,
                                                  formula = f_oth,
                                                  lib.size = lib.size,
                                                  future.seed=TRUE)

  message('Performing LRT for Negative Binomial model...')
  fitting_lrt[["NB_lrt"]] <- future.apply::future_lapply(lbic[["NB"]],
                                                      FUN = model_Chinb,
                                                      formula = f_oth,
                                                      lib.size = lib.size,
                                                      future.seed=TRUE)

  message('Performing LRT for Zero Inflated Poisson model...')
  fitting_lrt[["ZIP_lrt"]] <- future.apply::future_lapply(lbic[["ZIP"]],
                                                      FUN = model_Chizip,
                                                      formula = f_oth,
                                                      lib.size = lib.size,
                                                      future.seed=TRUE)

  message('Performing LRT for Zero Inflated Negtaive Binomial model...')
  fitting_lrt[["ZINB_lrt"]] <- future.apply::future_lapply(lbic[["ZINB"]],
                                                      FUN = model_Chizinb,
                                                      formula = f_oth,
                                                      lib.size = lib.size,
                                                      future.seed=TRUE)

  return(fitting_lrt)

}
