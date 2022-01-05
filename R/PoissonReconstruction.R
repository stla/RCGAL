#' @title Poisson surface reconstruction
#' @description Poisson reconstruction of a surface, from a cloud of 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param normals numeric matrix which stores the normals, one normal per row
#'   (it must have the same size as the \code{points} matrix); if you don't
#'   have normals, set \code{normals=NULL} (the default) and some normals will
#'   be computed with the help of \code{\link[Rvcg]{vcgUpdateNormals}}
#' @param spacing size parameter; smaller values increase the precision of the
#'   output mesh at the cost of higher computation time; set to \code{NULL}
#'   (the default) for a reasonable automatic value: an average spacing whose
#'   value will be displayed in a message and that you can also get in the
#'   \code{"spacing"} attribute of the output
#' @param sm_angle bound for the minimum facet angle in degrees
#' @param sm_radius relative bound for the radius of the surface Delaunay balls
#' @param sm_distance relative bound for the center-center distances
#'
#' @return A triangular mesh, of class \code{mesh3d} (ready for plotting
#'   with \strong{rgl}).
#'
#' @details See \href{https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/index.html}{Poisson Surface Recnstruction}.
#'
#' @export
#' @importFrom rgl tmesh3d addNormals
#' @importFrom Rvcg vcgUpdateNormals
#'
#' @examples library(RCGAL)
#' Psr_mesh <- PoissonReconstruction(SolidMobiusStrip)
#' library(rgl)
#' shade3d(Psr_mesh, color= "yellow")
#' wire3d(Psr_mesh, color = "black")
PoissonReconstruction <- function(
  points, normals = NULL, spacing = NULL,
  sm_angle = 20, sm_radius = 30, sm_distance = 0.375
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
    normals <- t(vcgUpdateNormals(points, silent = TRUE)[["normals"]][-4L, ])
  }else{
    storage.mode(normals) <- "double"
  }
  if(is.null(spacing)){
    spacing <- -1
  }else{
    stopifnot(isPositiveNumber(spacing))
  }
  stopifnot(isPositiveNumber(sm_angle))
  stopifnot(isPositiveNumber(sm_radius))
  stopifnot(isPositiveNumber(sm_distance))
  Psr <- Poisson_reconstruction_cpp(
    points, normals, spacing, sm_angle, sm_radius, sm_distance
  )
  out <-  addNormals(
    tmesh3d(t(Psr[["vertices"]]), t(Psr[["facets"]]), normals = NULL)
  )
  if(spacing == -1){
    message(sprintf(
      "Poisson reconstruction using average spacing: %s.",
      formatC(Psr[["spacing"]])
    ))
    attr(out, "spacing") <- Psr[["spacing"]]
  }
  out
}
