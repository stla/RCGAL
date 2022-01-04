#' @title Poisson surface reconstruction
#' @description Poisson reconstruction of a surface, from a cloud of 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param normals numeric matrix which stores the normals, one normal per row
#'   (it must have the same size as the \code{points} matrix); if you don't
#'   have normals, set \code{normals=NULL} (the default) and some normals will
#'   be computed with the help of \code{\link[Rvcg]{vcgUpdateNormals}}
#'
#' @return A triangular mesh, of class \code{mesh3d} (ready for plotting
#'   with \strong{rgl}).
#' @export
#' @importFrom rgl tmesh3d addNormals
#' @importFrom Rvcg vcgUpdateNormals
#'
#' @examples library(RCGAL)
#' Psr_mesh <- PoissonReconstruction(SolidMobiusStrip)
#' library(rgl)
#' wire3d(Psr_mesh, color = "black")
PoissonReconstruction <- function(
  points, normals = NULL
){
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  if(!is.null(normals) && (!is.matrix(normals) || !is.numeric(normals))){
    stop(
      "The `normals` argument must be `NULL` or a numeric matrix.",
      call. = TRUE
    )
  }
  dimension <- ncol(points)
  if(dimension != 3L){
    stop("The `points` matrix must have three columns.", call. = TRUE)
  }
  if(!is.null(normals)){
    dimension <- ncol(normals)
    if(dimension != 3L){
      stop("The `normals` matrix must have three columns.", call. = TRUE)
    }
    if(nrow(points) != nrow(normals)){
      stop(
        "The `points` matrix and the `normals` matrix must have the same ",
        "number of rows.", call. = TRUE
      )
    }
  }
  if(nrow(points) <= dimension){
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(any(is.na(points)) || (!is.null(normals) && any(is.na(normals)))){
    stop("Points or normals with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  if(is.null(normals)){
    normals <- vcgUpdateNormals(points, silent = TRUE)
  }else{
    storage.mode(normals) <- "double"
  }
  Psr <- Poisson_reconstruction_cpp(points, normals)
  addNormals(
    tmesh3d(t(Psr[["vertices"]]), t(Psr[["facets"]]), normals = NULL)
  )
}
