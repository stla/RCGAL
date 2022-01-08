#' @title Normals for a points could
#' @description Compute some normals for a 3D points cloud.
#'
#' @param points matrix of points, one point per row
#' @param nbNeighbors integer, number of neighbors used to compute the normals
#' @param method one of \code{"pca"} or \code{"jet"}
#'
#' @return A matrix of the same size as the \code{points} matrix, giving one
#'   unit normal for each point.
#' @export
getSomeNormals <- function(points, nbNeighbors, method = "pca"){
  method <- match.arg(method, c("pca", "jet"))
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  dimension <- ncol(points)
  if(dimension != 3L){
    stop("Points must be 3-dimensional.", call. = TRUE)
  }
  if(nrow(points) <= dimension){
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(any(is.na(points))){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  nbNeighbors <- as.integer(nbNeighbors)
  if(nbNeighbors <= 2L){
    stop("There must be at least two neighbors.", call. = TRUE)
  }
  if(method == "pca"){
    pca_normals_cpp(points, nbNeighbors)
  }else{
    jet_normals_cpp(points, nbNeighbors)
  }
}
