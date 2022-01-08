#' @useDynLib RCGAL, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

isPositiveNumber <- function(x){
  is.numeric(x) && length(x) == 1L && x > 0 && !is.na(x)
}

isBoolean <- function(x){
  is.logical(x) && length(x) == 1L && !is.na(x)
}
