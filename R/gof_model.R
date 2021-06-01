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
#' @param BPPARAM configuration parameter related to the method of parallel execution.
#' For further information on how to set-up parallel execution refer to
#' \code{BiocParallel} vignette.
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
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
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
#' scData_models <- fit_models(counts=scData$counts, cexpr=scData$covariates, lib.size=scData$lib_size,
#' BPPARAM=bpparam())
#' scData_bicvals <- model_bic(scData_models)
#' scData_least.bic <- lbic_model(scData_bicvals, scData$counts)
#'
#' scData_gof <- gof_model(scData_least.bic, cexpr=scData$covariates, lib.size=scData$lib_size,
#' BPPARAM=bpparam())


gof_model <- function(lbic, cexpr, lib.size,
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

  f_oth <- as.formula(formula)
  f_nb <- as.formula(paste(formula, ' + offset(log(lib.size))'))


  #Prepare input data as lists to be inputted to the function
  if(identical(covariates, NULL)) {
    lbic[["P"]] <- lapply(lbic[["P"]], function (x) as.list(as.data.frame(x)))
  } else {
    lbic[["P"]] <- lapply(lbic[["P"]], function (x) cbind(x,cexpr))
  }

  if(identical(covariates, NULL)) {
    lbic[["NB"]] <- lapply(lbic[["NB"]], function (x) as.list(as.data.frame(x)))
  } else {
    lbic[["NB"]] <- lapply(lbic[["NB"]], function (x) cbind(x,cexpr))
  }

  if(identical(covariates, NULL)) {
    lbic[["ZIP"]] <- lapply(lbic[["ZIP"]], function (x) as.list(as.data.frame(x)))
  } else {
    lbic[["ZIP"]] <- lapply(lbic[["ZIP"]], function (x) cbind(x,cexpr))
  }

  if(identical(covariates, NULL)) {
    lbic[["ZINB"]] <- lapply(lbic[["ZINB"]], function (x) as.list(as.data.frame(x)))
  } else {
    lbic[["ZINB"]] <- lapply(lbic[["ZINB"]], function (x) cbind(x,cexpr))
  }


  #Fit the four distributions
  #Poisson distribution
  model_Chipoi <- function(data, formula, lib.size){

    1-pchisq(summary(stats::glm(formula, data, offset=log(lib.size), family = "poisson"))$deviance,
             df= summary(stats::glm(formula, data, offset=log(lib.size), family = "poisson"))$df.residual)
  }


  #NB distribution
  model_Chinb <- function(data, formula, lib.size){

    1-pchisq(summary(MASS::glm.nb(formula, data))$deviance,
             df= summary(MASS::glm.nb(formula, data))$df.residual)
  }


  #ZIP distribution
  model_Chizip <- function(data, formula, lib.size){

    1-emdbook::pchibarsq(2*(logLik(pscl::zeroinfl(formula, data, offset=log(lib.size), dist = "poisson"))-logLik(stats::glm(formula, data, offset=log(lib.size), family = "poisson"))),df=1,mix=0.5,lower.tail=TRUE, log.p = FALSE)

  }


  #ZINB distribution
  model_Chizinb <- function(data, formula, lib.size){

    1-emdbook::pchibarsq(2*(logLik(pscl::zeroinfl(formula, data, offset=log(lib.size), dist = "negbin"))-logLik(MASS::glm.nb(formula, data))),df=1,mix=0.5,lower.tail=TRUE, log.p = FALSE)
  }


  fitting_lrt <- list()

  message('Performing LRT for Poisson model...')
  if(identical(covariates, NULL)) {
    fitting_lrt[["P_lrt"]] <- lapply(lbic[["P"]],
                                     FUN = model_Chipoi,
                                     formula = f_oth,
                                     lib.size = lib.size)
  } else {
    fitting_lrt[["P_lrt"]] <- lapply(lbic[["P"]],
                                     FUN = model_Chipoi,
                                     formula = f_oth,
                                     lib.size = lib.size)
  }

  message('Performing LRT for Negative Binomial model...')
  if(identical(covariates, NULL)) {
    fitting_lrt[["NB_lrt"]] <- BiocParallel::bplapply(lbic[["NB"]],
                                                      FUN = model_Chinb,
                                                      BPPARAM = BPPARAM,
                                                      formula = f_oth,
                                                      lib.size = lib.size)
  } else {
    fitting_lrt[["NB_lrt"]] <- BiocParallel::bplapply(lbic[["NB"]],
                                                      FUN = model_Chinb,
                                                      BPPARAM = BPPARAM,
                                                      formula = f_oth,
                                                      lib.size = lib.size)
  }

  message('Performing LRT for Zero Inflated Poisson model...')
  if(identical(covariates, NULL)) {
    fitting_lrt[["ZIP_lrt"]] <- BiocParallel::bplapply(lbic[["ZIP"]],
                                                      FUN = model_Chizip,
                                                      BPPARAM = BPPARAM,
                                                      formula = f_oth,
                                                      lib.size = lib.size)
  } else {
    fitting_lrt[["ZIP_lrt"]] <- BiocParallel::bplapply(lbic[["ZIP"]],
                                                       FUN = model_Chizip,
                                                       BPPARAM = BPPARAM,
                                                       formula = f_oth,
                                                       lib.size = lib.size)
  }

  message('Performing LRT for Zero Inflated Negtaive Binomial model...')
  if(identical(covariates, NULL)) {
    fitting_lrt[["ZINB_lrt"]] <- BiocParallel::bplapply(lbic[["ZINB"]],
                                                      FUN = model_Chizinb,
                                                      BPPARAM = BPPARAM,
                                                      formula = f_oth,
                                                      lib.size = lib.size)
  } else {
    fitting_lrt[["ZINB_lrt"]] <- BiocParallel::bplapply(lbic[["ZINB"]],
                                                        FUN = model_Chizinb,
                                                        BPPARAM = BPPARAM,
                                                        formula = f_oth,
                                                        lib.size = lib.size)
  }

  return(fitting_lrt)

}
