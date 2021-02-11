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
#' @param formula A regression formula to fit the covariates in the ZINB GLM.
#'
#' @param workers Number of workers to be used in parallel computation
#' using \code{future.apply}, with argument \code{multisession}.
#'
#' @param seed Seed number to be used in parallel computation
#' using \code{future.apply}, with argument \code{multisession}
#'
#' @param distr A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#'
#' @export
#'
#' @import stats
#' @importFrom parallelly availableCores
#' @importFrom Matrix Matrix
#' @importFrom future.apply future_apply
#' @importFrom pscl zeroinfl
#' @importFrom MASS glm.nb
#'
#' @return List object containing the significant gene indices from the KS test,
#' their adjusted p-values

fit_models <- function(counts, cexpr, formula=NULL, workers=NULL, seed=NULL, distr=NULL){

  if(is.null(workers)) {
    workers <- min(4, parallelly::availableCores())
  }
  if(is.null(seed)) {
    seed <- 0xBEEF
  }

  #Formulate a simple additive model using all the covariates in 'cexpr'
  covariates <- names(cexpr)
  if(is.null(formula)) {
    message(sprintf("Formulating the default additive model..."))
    formula_uni <- 'x ~ 1'
    formula_zi <- 'x ~ 1 | 1'
    if(!identical(covariates, character(0))) {
      for (covar in covariates) {
        formula_uni <- paste(formula_uni, sprintf(' + %s', covar))
        formula_zi <- paste(formula_zi, sprintf(' + %s', covar))
      }
    }
  }

  f_poi <- as.formula(formula_uni)
  f_nb <- as.formula(paste(formula_uni, ' + offset(log(lib.size))'))
  f_zip <- as.formula(formula_zi)
  f_zinb <- as.formula(formula_zi)


  #set-up multisession
  plan(multisession, workers = workers)

  #Calculate library size for each cell
  lib.size <- apply(counts,2, function(x) sum(x))

  #Function to fit the 4 distributions
  model_poi <- function(formula, data, lib.size){

    stats::glm(formula, data, offset=log(lib.size), family = "poisson")
  }


  model_nb <- function(formula, data, lib.size){
    library(MASS)
    try(glm.nb(formula, data), silent = TRUE)
  }


  model_zip <- function(formula, data, lib.size){
    library(pscl)
    try(zeroinfl(formula, data, offset=log(lib.size), dist = "poisson"), silent = TRUE)
  }


  model_zinb <- function(formula, data, lib.size){
    library(pscl)
    try(zeroinfl(formula, data, offset=log(lib.size), dist = "negbin"), silent = TRUE)
  }


  fitting <- list()

  if(is.null(model2fit) || 1 %in% model2fit) {
    message('Fitting data with Poisson model...')
    if(identical(covariates, character(0))) {
      fitting[["P"]] <- glm
    } else {
      fitting[["P"]] <- stan_glmer(f12,
                                   family = poisson,
                                   data = gexpr,
                                   offset = exposure,
                                   cores = nCores,
                                   seed = seed,
                                   refresh = 0)
    }
  }

  if(is.null(model2fit) || 2 %in% model2fit) {
    message('Fitting data with Negative Binomial model...')
    if(identical(covariates, character(0))) {
      fitting[["NB"]] <-   stan_glm(f12,
                                    family = neg_binomial_2,
                                    data = gexpr,
                                    offset = exposure,
                                    cores = nCores,
                                    seed = seed,
                                    refresh = 0)
    } else {
      fitting[["NB"]] <- stan_glmer(f12,
                                    family = neg_binomial_2,
                                    data = gexpr,
                                    offset = exposure,
                                    cores = nCores,
                                    seed = seed,
                                    refresh = 0)
    }
  }

  if(is.null(model2fit) || 3 %in% model2fit) {
    message('Fitting data with Zero-Inflated Poisson model...')
    myprior_3 <- get_prior(bf(f34, zi ~ 1),
                           family = zero_inflated_poisson(),
                           data = gexpr)
    myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
    fitting[["ZIP"]] <- brm(bf(f34, zi ~ 1),
                            family = zero_inflated_poisson(),
                            data = gexpr,
                            prior = myprior_3,
                            control = list(adapt_delta = adapt_delta),
                            cores = nCores,
                            seed = seed,
                            refresh = 500)
  }

  if(is.null(model2fit) || 4 %in% model2fit) {
    message('Fitting data with Zero-Inflated Negative Binomial model...')
    myprior_4 <- get_prior(bf(f34, zi ~ 1),
                           family = zero_inflated_negbinomial(),
                           data = gexpr)
    myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
    fitting[["ZINB"]] <- brm(bf(f34, zi ~ 1),
                             family = zero_inflated_negbinomial(),
                             data = gexpr,
                             control = list(adapt_delta = adapt_delta),
                             prior = myprior_4,
                             cores = nCores,
                             seed = seed,
                             refresh = 500)
  }

  return(fitting)
}

