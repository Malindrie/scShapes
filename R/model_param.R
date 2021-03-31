#' model_param
#'
#' This function returns model parameters based on the best fit distribution
#' as selected by \code{distr_fit} and models fitted by \code{fit_models}
#'
#' @param fit.model A list of models fitted by 'glm' from \code{fit_models}.
#'
#' @param gof.res A list of selected model distributions for genes \code{select_model}.
#'
#' @param model A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#'
#' @export
#'
#' @importFrom VGAM predict
#'
#' @return A list of model parameters estimated. Estimated model parameters include
#' mean (for all 4 models), theta (over-dispersion parameter for NB & ZINB models), pi (zero-inflation parameter
#' for ZIP & ZINB models).
#'
#' @examples
#'
#' data(scData)
#'
#' # apply the model_param function to extract the parameters of the best fit
#' # model obtained by running the select_model function
#'
#' scData_models <- fit_models(counts=scData$counts, cexpr=scData$covariates, lib.size=scData$lib_size)
#' scData_bicvals <- model_bic(scData_models)
#' scData_least.bic <- lbic_model(scData_bicvals, scData$counts)
#' scData_gof <- gof_model(scData_least.bic, cexpr=scData$covariates, lib.size=scData$lib_size)
#' scData_fit <- select_model(scData_gof)
#'
#' scData_params <- model_param(scData_models, scData_fit, model=NULL)

model_param <- function(fit.model, gof.res,
                        model=NULL){

  param_poi <- function(m1){

    mean_mu <- predict(m1, type='response')

    return(mean_mu)
  }


  param_nb <- function(m1){

    mean_mu <- predict(m1, type='response')
    theta <- m1$theta

    theta_mean <- c(theta, mean_mu)
    names(theta_mean)[1] <- "theta"
    return(theta_mean)
  }


  param_zip <- function(m1){

    pi_ML = predict(m1, type = "zero")[1]
    mean_ML = predict(m1, type = "count")

    pi_mean <- c(pi_ML, mean_ML)
    names(pi_mean)[1] <- "pi"
    return(pi_mean)
  }


  param_zinb <- function(m1){

    pi_ML = predict(m1, type = "zero")[1]
    theta_ML = m1$theta
    mean_ML = predict(m1, type = "count")

    pi_theta_mean <- c(pi_ML, theta_ML, mean_ML)
    names(pi_theta_mean)[seq_along(1:2)] <- c("pi", "theta")
    return(pi_theta_mean)
  }


  m1.param <- list()


  if(is.null(model) || 1 %in% model) {

    P_params <- fit.model[["P"]][Reduce(intersect, list(names(fit.model[["P"]]), gof.res[["P_genes"]]))]

    m1.param[["P"]] <- lapply(P_params,  param_poi)
  }


  if(is.null(model) || 2 %in% model) {

    NB_params <- fit.model[["NB"]][Reduce(intersect, list(names(fit.model[["NB"]]), gof.res[["NB_genes"]]))]

    m1.param[["NB"]] <- lapply(NB_params, param_nb)

  }


  if(is.null(model) || 3 %in% model) {

    ZIP_params <- fit.model[["ZIP"]][Reduce(intersect, list(names(fit.model[["ZIP"]]), gof.res[["ZIP_genes"]]))]

    m1.param[["ZIP"]] <- lapply(ZIP_params, param_zip)

  }


  if(is.null(model) || 4 %in% model) {

    ZINB_params <- fit.model[["ZINB"]][Reduce(intersect, list(names(fit.model[["ZINB"]]), gof.res[["ZINB_genes"]]))]

    m1.param[["ZINB"]] <- lapply(ZINB_params, param_zinb)

  }

  return(m1.param )

}
