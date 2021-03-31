#' filter_counts
#'
#' This function is used to preprocess matrix of read counts to
#' only keep genes with a certain number of nonzero entries.
#'
#' @param counts A non-negative integer matrix of scRNA-seq raw read counts.
#' The rows of the matrix are genes and columns are samples/cells.
#'
#' @param perc.zero A numeric value between 0 and 1 that represents the
#' proportion of zeros per gene in the processed dataset.
#'
#' @export
#'
#' @importFrom Matrix Matrix
#'
#' @return An object of class \code{\link{Matrix}} with genes removed if
#' they have more than \code{perc.zero} zeros.
#'
#' @examples
#'
#' # load toy  example data
#'
#' data(scData)
#'
#' # apply the filter_counts function to filter out genes if they have
#' # more than 10% zero
#'
#' scData_filt <- filter_counts(scData$counts, perc.zero = 0.1)


filter_counts <- function(counts, perc.zero = 0.1){
  # Invalid input control
  if(!is.matrix(counts) & !is.data.frame(counts) & class(counts)[1] != "dgCMatrix")
    stop("Wrong data type of 'counts'")
  if(sum(is.na(counts)) > 0)
    stop("NA detected in 'counts'");gc();
  if(sum(counts < 0) > 0)
    stop("Negative value detected in 'counts'");gc();
  if(all(counts == 0))
    stop("All elements of 'counts' are zero");gc();
  if(any(colSums(counts) == 0))
    warning("Library size of zero detected in 'counts'");gc()

  # Preprocessing
  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"
  if(any(rowSums(counts) == 0))
    message("Removing ", sum(rowSums(counts) == 0), " rows of genes with all zero counts")
  counts <- counts[rowSums(counts) != 0,]

  counts <- counts[apply(counts, 1, function(x){sum(x == 0)}) < ncol(counts)*(1-perc.zero),]

  #convert data matrix to sparse matrix
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  return(counts)
}
