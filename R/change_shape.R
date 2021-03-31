#' change_shape
#'
#' This function returns a list of genes changing shape between conditions and
#' the list of genes changing distribution from a unimodal to distribution to
#' a zero-inflated distribution
#'
#' @param shapes_genes A dataframe consisting of distributions followed by
#' each gene passing the KS test. Where rows are genes and each column the
#' treatment condition.
#'
#' @export
#'
#' @return A list of two lists with genes changing distribution shape between
#' conditions. The list "All" contains all the genes changing distribution.
#' The list "Uni2ZI" contains the genes changing distribution form a unimodal
#' distribution in one condition to a zero-inflated distribution in another
#' condition.


change_shape <- function(shapes_genes){

        genes_CS <- list()

        genes_CS[["All"]] <- c(rownames(shapes_genes)[(shapes_genes[,1] == "Po" & shapes_genes[,2] == "NB")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "Po" & shapes_genes[,2] == "ZIP")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "Po" & shapes_genes[,2] == "ZINB")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "NB" & shapes_genes[,2] == "Po")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "NB" & shapes_genes[,2] == "ZIP")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "NB" & shapes_genes[,2] == "ZINB")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "ZIP" & shapes_genes[,2] == "Po")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "ZIP" & shapes_genes[,2] == "NB")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "ZIP" & shapes_genes[,2] == "ZINB")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "ZINB" & shapes_genes[,2] == "Po")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "ZINB" & shapes_genes[,2] == "NB")],
                               rownames(shapes_genes)[(shapes_genes[,1] == "ZINB" & shapes_genes[,2] == "ZIP")]
        )


        genes_CS[["Uni2ZI"]] <- c(rownames(shapes_genes)[(shapes_genes[,1] == "Po" & shapes_genes[,2] == "ZIP")],
                                  rownames(shapes_genes)[(shapes_genes[,1] == "Po" & shapes_genes[,2] == "ZINB")],
                                  rownames(shapes_genes)[(shapes_genes[,1] == "NB" & shapes_genes[,2] == "ZIP")],
                                  rownames(shapes_genes)[(shapes_genes[,1] == "NB" & shapes_genes[,2] == "ZINB")],
                                  rownames(shapes_genes)[(shapes_genes[,1] == "ZIP" & shapes_genes[,2] == "Po")],
                                  rownames(shapes_genes)[(shapes_genes[,1] == "ZIP" & shapes_genes[,2] == "NB")],
                                  rownames(shapes_genes)[(shapes_genes[,1] == "ZINB" & shapes_genes[,2] == "Po")],
                                  rownames(shapes_genes)[(shapes_genes[,1] == "ZINB" & shapes_genes[,2] == "NB")]
                      )


        return(genes_CS)

}
