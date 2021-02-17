#' lbic_model
#'
#' This function is used to select the best fit model for each gene based on the
#' least BIC value
#'
#' @param bic.value A dataframe of BIC values from fitting GLM using
#' error distributions P, NB, ZIP, ZINB; the output from \code{model_bic}.
#'
#' @param counts A non-negative integer matrix of scRNA-seq filtered read counts
#' containing genes belonging to the family of ZINB distributions selected from
#' \code{ks_test}.
#' The rows of the matrix are genes and columns are samples/cells.
#'
#' @export
#'
#' @import magrittr
#'
#' @return A list of genes chosen to be following one of the 4 distributions
#' P, NB, ZIP, ZINB based on the least BIC value and the corresponding subset
#' of counts from \code{filter_counts}


lbic_model <- function(bic.value, counts){

  #subset counts that passed the KS test and was passed on to calculate the BIC
  counts_bic <- counts[rownames(counts) %in% rownames(bic.value),]
  counts_list <- lapply(as.list(1:dim(counts_bic)[1]), function(x) counts_bic[x[1],])
  names(counts_list) <- rownames(counts_bic)

  #return column name of min value for each row
  BIC_colmax <- colnames(bic.value)[apply(bic.value,1,function(x) which.min(x[x>0]))]

  #create data frame
  BIC_dataframe <- as.data.frame(BIC_colmax)

  #Get the total count of genes follwoing each distribution
  BIC_dataframe %>% group_by(BIC_colmax) %>% tally()


  #Create data frame of minimum BIC values with gene names
  BIC_genename <- as.data.frame(BIC_dataframe)
  row.names(BIC_genename) <- rownames(bic.value)

  #grep the names of genes following each distribution
  poi_genes <- rownames(BIC_genename)[grep(BIC_genename[,1], pattern = "^P_bic")]
  nb_genes <- rownames(BIC_genename)[grep(BIC_genename[,1], pattern = "^NB_bic")]
  zip_genes <- rownames(BIC_genename)[grep(BIC_genename[,1], pattern = "^ZIP_bic")]
  zinb_genes <- rownames(BIC_genename)[grep(BIC_genename[,1], pattern = "^ZINB_bic")]

  least_bic <- list()

  least_bic[["P"]] <- counts_list[names(counts_list) %in% poi_genes]
  least_bic[["NB"]] <- counts_list[names(counts_list) %in% nb_genes]
  least_bic[["ZIP"]] <- counts_list[names(counts_list) %in% zip_genes]
  least_bic[["ZINB"]] <- counts_list[names(counts_list) %in% zinb_genes]

  return(least_bic)
}
