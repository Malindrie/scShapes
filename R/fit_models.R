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
#' @param model A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#'
#' @param BPPARAM configuration parameter related to the method of parallel execution.
#' For further information on how to set-up parallel execution refer to
#' \code{BiocParallel} vignette.
#'
#' @export
#'
#' @importFrom stats glm
#' @importFrom stats as.formula
#' @importFrom Matrix Matrix
#' @importFrom pscl zeroinfl
#' @importFrom MASS glm.nb
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
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
#' library(BiocParallel)
#' scData_models <- fit_models(counts=scData$counts, cexpr=scData$covariates,
#' lib.size=scData$lib_size, BPPARAM=bpparam())


fit_models <- function(counts, cexpr, lib.size,
                       formula=NULL,model=NULL,
                       BPPARAM){


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

  #convert data to lists
  if(!identical(covariates, NULL)) {
    gexpr <- apply(counts, 1, function (x) cbind(x,cexpr))
  } else{
    gexpr <- apply(counts, 1, function (x) as.list(as.data.frame(x)))
  }

  #Fit the four distributions
  model_poi <- function(data, formula, lib.size){
    tryCatch(stats::glm(formula, data, offset=log(lib.size), family = "poisson"), error = function(e) {print(e$message)})
  }


  model_nb <- function(data, formula, lib.size){
    tryCatch(MASS::glm.nb(formula, data), error = function(e) {print(e$message)})
  }


  model_zip <- function(data, formula, lib.size){
    tryCatch(pscl::zeroinfl(formula, data, offset=log(lib.size), dist = "poisson"), error = function(e) {print(e$message)})
  }


  model_zinb <- function(data, formula, lib.size){
    tryCatch(pscl::zeroinfl(formula, data, offset=log(lib.size), dist = "negbin"), error = function(e) {print(e$message)})
  }


  fitting <- list()

  if(is.null(model) || 1 %in% model) {
    message('Fitting data with Poisson model...')
    if(identical(covariates, NULL)) {
      fitting[["P"]] <-   lapply(gexpr,
                                 FUN = model_poi,
                                 formula = f_oth,
                                 lib.size = lib.size)
    } else {
      fitting[["P"]] <- lapply(gexpr,
                               FUN = model_poi,
                               formula = f_oth,
                               lib.size = lib.size)
    }
  }

  if(is.null(model) || 2 %in% model) {
    message('Fitting data with Negative Binomial model...')
    if(identical(covariates, NULL)) {
      fitting[["NB"]] <-   BiocParallel::bplapply(gexpr,
                                                  FUN = model_nb,
                                                  BPPARAM = BPPARAM,
                                                  formula = f_nb,
                                                  lib.size = lib.size)
    } else {
      fitting[["NB"]] <- BiocParallel::bplapply(gexpr,
                                                FUN = model_nb,
                                                BPPARAM = BPPARAM,
                                                formula = f_nb,
                                                lib.size = lib.size)
    }
  }

  if(is.null(model) || 3 %in% model) {
    message('Fitting data with Zero Inflated Poisson model...')
    if(identical(covariates, NULL)) {
      fitting[["ZIP"]] <-   BiocParallel::bplapply(gexpr,
                                                   FUN = model_zip,
                                                   BPPARAM = BPPARAM,
                                                   formula = f_oth,
                                                   lib.size = lib.size)
    } else {
      fitting[["ZIP"]] <- BiocParallel::bplapply(gexpr,
                                                 FUN = model_zip,
                                                 BPPARAM = BPPARAM,
                                                 formula = f_oth,
                                                 lib.size = lib.size)
    }
  }

  if(is.null(model) || 4 %in% model) {
    message('Fitting data with Zero Inflated Negative Binomial model...')
    if(identical(covariates, NULL)) {
      fitting[["ZINB"]] <-   BiocParallel::bplapply(gexpr,
                                                    FUN = model_zinb,
                                                    BPPARAM = BPPARAM,
                                                    formula = f_oth,
                                                    lib.size = lib.size)
    } else {
      fitting[["ZINB"]] <- BiocParallel::bplapply(gexpr,
                                                  FUN = model_zinb,
                                                  BPPARAM = BPPARAM,
                                                  formula = f_oth,
                                                  lib.size = lib.size)
    }
  }

  return(fitting)

}

