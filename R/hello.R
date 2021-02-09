#' My Hello World Function
#'
#' @param x The name of the person to say hi to
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @examples
#' hello("John")
#' \dontrun{
#' hello("Steve")
#' }
hello <- function(x) {
  print(paste0("Hello ", x, ", this is the world!"))
}
