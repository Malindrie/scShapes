#' model_bic
#'
#' This function is used to calculate the Bayesian information criterion of the models fitted in
#' \code{fit_models}.
#'
#' @param fit_list A list of models fitted from \code{fit_models}
#'
#' @export
#'
#' @return A dataframe containing the BIC values for each distribution type (P, NB, ZIP, ZINB).


model_bic <- function(fit_models){

  model_BIC <- function(z){
    if(!(class(z)=="try-error")){
      stats::BIC(z)
    }
    else{
      "NA"
    }
  }


  bic_p <- t(as.data.frame(lapply(fit_models$P, model_BIC)))
  bic_nb <- t(as.data.frame(lapply(fit_models$NB, model_BIC)))
  bic_zip <- t(as.data.frame(lapply(fit_models$ZIP, model_BIC)))
  bic_zinb <- t(as.data.frame(lapply(fit_models$ZINB, model_BIC)))

  bic <- cbind(bic_p, bic_nb, bic_zip, bic_zinb)
  colnames(bic) <- c("P_bic", "NB_bic", "ZIP_bic", "ZINB_bic")

  return(bic)

}

