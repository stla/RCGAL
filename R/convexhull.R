#' @title Convex hull
#' @description Convex hull of a set of 2D or 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#'
#' @return The convex hull.
#' \itemize{
#'   \item \strong{If the dimension is 2}, the returned value is a list with
#'         five fields:
#'         \describe{
#'           \item{\emph{verticesIds}}{The indices (i.e. the row numbers) of
#'                the points which form the convex hull.}
#'           \item{\emph{vertices}}{A matrix giving the coordinates of the
#'                points which from the convex hull; they are given by rows in
#'                the order defined by \code{verticesIds}.}
#'           \item{\emph{edges}}{A matrix of integers with two columns which
#'                represents the edges; each row provides two indices, the
#'                the indices of the two points which form the edge.}
#'           \item{\emph{surface}}{A number, the surface of the convex hull.}
#'           \item{\emph{perimeter}}{A number, the perimeter of the convex hull.}
#'         }
#'   \item \strong{If the dimension is 3}, the returned value is a list with
#'         five fields:
#'         \describe{
#'           \item{\emph{vertices}}{A list which represents the vertices of the
#'                convex hull. This is a list of lists, each sublist represents
#'                one vertice, by giving its index and its coordinates.}
#'           \item{\emph{edges}}{An integer matrix with three columns,
#'                representing the edges of the convex hull. So each row is
#'                composed of three integers; the two first ones are the indices
#                 of the two points which form the edge, and the third one, in
#'                the column named \code{border}, is a \code{0/1} indicator of
#'                whether the edge is a border edge.}
#'           \item{\emph{faces}}{A matrix of integers with three columns which
#'                represents the faces (these are triangles); each row provides
#'                the indices of the three points which form the face. This
#'                matrix has three attributes: \emph{areas}, which provides the
#'                areas of the faces, \emph{normals}, which provides the normals
#'                of the faces, and \emph{circumcenters}, which provides the
#'                circumcenters of the faces.}
#'           \item{\emph{surface}}{A number, the surface of the convex hull.}
#'           \item{\emph{volume}}{A number, the volume of the convex hull.}
#'         }
#' }
#' @export
#'
#' @examples library(RCGAL)
convexhull <- function(
  points
){
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  dimension <- ncol(points)
  if(!dimension %in% c(2L, 3L)){
    stop("The dimension must be 2 or 3.", call. = TRUE)
  }
  if(!identical(anyDuplicated(points), 0L)){
    stop("There are some duplicated points.", call. = TRUE)
  }
  if(nrow(points) <= dimension){
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(any(is.na(points))){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  if(dimension == 2L){
    out <- cxhull2d_cpp(points)
  }else{
    out <- cxhull3d_cpp(points)
  }
  class(out) <- "cxhull"
  attr(out, "points") <- points
  out
}
