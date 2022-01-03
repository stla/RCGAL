#' @title Advancing front surface reconstruction
#' @description Reconstruction of a surface from a cloud of 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#'
#' @return A triangular mesh, of class \code{mesh3d} (ready for plotting
#'   with \strong{rgl}).
#' @export
#' @importFrom rgl tmesh3d
#' @importFrom Rvcg vcgClean
#'
#' @examples library(RCGAL)
#' data(bunny, package = "onion")
#' mesh <- AFSreconstruction(bunny)
#' library(rgl)
#' shade3d(mesh, color = "firebrick")
AFSreconstruction <- function(
  points
){
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  dimension <- ncol(points)
  if(dimension != 3L){
    stop("The dimension must be 3.", call. = TRUE)
  }
  if(nrow(points) <= dimension){
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(any(is.na(points))){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  afs <- AFSreconstruction_cpp(points)
  rglmesh <-
    tmesh3d(afs[["vertices"]], afs[["triangles"]], normals = afs[["normals"]])
  vcgClean(rglmesh, sel = c(0, 7), silent = TRUE)
}
