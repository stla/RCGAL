#' @title Delaunay tessellation
#' @description Delaunay tessellation of a set of 2D or 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param elevation if points are three-dimensional and \code{elevation=TRUE},
#'   then the function performs an elevated two-dimensional Delaunay
#'   triangulation, using the \code{z} coordinates of the points for the
#'   elevations; see the example
#'
#' @return The Delaunay tessellation.
#' \itemize{
#'   \item \strong{If the dimension is 2}, the returned value is a list with
#'         two fields: \code{faces} and \code{edges}. The \code{faces} field
#'         contains an integer matrix with three columns; each row represents a
#'         triangle whose each vertex is given by the index (row number) of
#'         this point in the \code{points} matrix. The \code{edges} field
#'         also contains an integer matrix with three columns. The two first
#'         integers of a row are the indices of the two points which form the
#'         edge. The third column, named \code{border}, only contains some
#'         zeros and some ones; a border (exterior) edge is labelled by a
#'         \code{1}.
#'   \item \strong{If the dimension is 3}, the returned value is a list with
#'         four fields: \code{cells}, \code{facets}, \code{edges}, and
#'         \code{volume}. The \code{cells} field represents the tetrahedra
#'         which form the tessellation. The \code{facets} field represents
#'         the faces of these tetrahedra, some triangles. The \code{edges}
#'         field represents the edges of these triangles. The \code{volume}
#'         field provides only one number, the volume of the tessellation,
#'         in other words the volume of the convex hull of the given points.
#'         The \code{cells} field is a list of lists. Each sublist is composed
#'         of three fields: \code{cell} provides the indices of the four
#'         vertices of the corresponding tetrahedron, \code{faces} provides the
#'         indices of the four faces of the tetrahedron, that is to say the row
#'         number of the \code{facets} field which represents this face, and
#'         finally there is a \code{volume} field which provides the volume of
#'         the tetrahedron. The \code{facets} field is an integer matrix with
#'         four columns. The three first integers of a row are the indices of
#'         the points which form the corresponding facet. The fourth column,
#'         named \code{onhull} is composed of zeros and ones only, and a
#'         \code{1} means that the corresponding facet lies on the convex hull
#'         of the points. The \code{edges} field contains an integer matrix
#'         with three columns. Each row represents an edge, given by the two
#'         indices of the points which form this edge, and the third integer,
#'         in the column named \code{onhull} is a \code{0/1} indicator of
#'         whether the edge lies on the convex hull. Finally the \code{volume}
#'         field provides only one number, the volume of the tessellation (i.e.
#'         the volume of the convex hull of the points).
#'   \item \strong{If} \code{elevation=TRUE}, the returned value is a list with
#'         four fields: \code{mesh}, \code{edges}, \code{faceVolumes}, and
#'         \code{volume}. The \code{mesh} field is an object of class
#'         \code{mesh3d}, ready for plotting with the \strong{rgl} package. The
#'         \code{edges} field provides the indices of the edges, given as an
#'         integer matrix with two columns. The \code{faceVolumes} is a numeric
#'         vector, it provides the volumes under the faces that can be found in
#'         the \code{mesh} field. Finally the \code{volume} field provides the
#'         sum of these volumes, that is to say the total volume under the
#'         triangulated surface.
#' }
#' @export
#' @importFrom rgl tmesh3d addNormals
#'
#' @examples library(RCGAL)
#' # elevated Delaunay triangulation
#' f <- function(x, y){
#'   2 * exp(-(x^2 + y^2)) # integrate to 2pi
#' }
#' x <- y <- seq(-4, 4, length.out = 50)
#' grd <- transform(expand.grid(x = x, y = y), z = f(x, y))
#' del <- delaunay(as.matrix(grd), elevation = TRUE)
#' # `del` is a list; its first component is a mesh representing the surface:
#' mesh <- del[["mesh"]]
#' library(rgl)
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(mesh, color = "limegreen")
#' wire3d(mesh)
#' # in `del` you can also found the volume under the surface, which should
#' #   approximate the integral of the function:
#' del[["volume"]]
delaunay <- function(
  points, elevation = FALSE
){
  stopifnot(isBoolean(elevation))
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  if(!identical(anyDuplicated(points), 0L)){
    stop("There are some duplicated rows in the `points` matrix.", call. = TRUE)
  }
  dimension <- ncol(points)
  if(!dimension %in% c(2L, 3L)){
    stop("The dimension must be 2 or 3.", call. = TRUE)
  }
  if(elevation && dimension == 2L){
    stop(
      "If you set `elevation=TRUE`, you must provide three-dimensional points.",
      call. = TRUE
    )
  }
  if(nrow(points) <= dimension){
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(any(is.na(points))){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  if(dimension == 2L){
    out <- del2d_cpp(points)
  }else{
    if(elevation){
      del <- del2d_xy_cpp(points)
      out <- list(
        mesh = addNormals(
          tmesh3d(
            vertices = t(points),
            indices = t(del[["faces"]])
          )
        ),
        edges = del[["edges"]],
        faceVolumes = attr(del[["faces"]], "volumes"),
        volume = del[["volume"]]
      )
    }else{
      out <- del3d_cpp(points)
    }
  }
  class(out) <- "delaunay"
  attr(out, "points") <- points
  out
}


#' @title Plot 2D Delaunay tessellation
#' @description Plot a 2D Delaunay tessellation.
#'
#' @param tessellation the output of \code{\link{delaunay}} used with 2D points
#' @param border the color of the borders of the triangles; \code{NULL} for
#'   no borders
#' @param color controls the filling colors of the triangles, either
#'   \code{FALSE} for no color, \code{"random"} to use
#'   \code{\link[randomcoloR]{randomColor}}, or \code{"distinct"} to use
#'   \code{\link[randomcoloR]{distinctColorPalette}}
#' @param hue,luminosity if \code{color = "random"}, these arguments are passed
#'   to \code{\link[randomcoloR]{randomColor}}
#' @param lty,lwd graphical parameters
#' @param ... arguments passed to \code{\link{plot}}
#'
#' @return No value, just renders a 2D plot.
#' @export
#' @importFrom randomcoloR randomColor distinctColorPalette
#' @importFrom graphics plot polygon par segments
#'
#' @examples # random points in a square
#' set.seed(314)
#' library(tessellation)
#' library(uniformly)
#' square <- rbind(
#'   c(-1, 1), c(1, 1), c(1, -1), c(-1, -1)
#' )
#' ptsin <- runif_in_cube(10L, d = 2L)
#' pts <- rbind(square, ptsin)
#' d <- delaunay(pts)
#' opar <- par(mar = c(0, 0, 0, 0))
#' plotDelaunay2D(
#'   d, xlab = NA, ylab = NA, asp = 1, color = "random", luminosity = "dark"
#' )
#' par(opar)
plotDelaunay2D <- function(
  tessellation, border = "black", color = "distinct", hue = "random",
  luminosity = "light", lty = par("lty"), lwd = par("lwd"), ...
){
  if(!inherits(tessellation, "delaunay")){
    stop(
      "The argument `tessellation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  vertices <- attr(tessellation, "points")
  if(ncol(vertices) != 2L){
    stop(
      sprintf("Invalid dimension (%d instead of 2).", ncol(vertices)),
      call. = TRUE
    )
  }
  plot(vertices, type = "n", ...)
  if(!isFALSE(color)){
    color <- match.arg(color, c("random", "distinct"))
    trianglesIds <- tessellation[["faces"]]
    ntriangles <- nrow(trianglesIds)
    if(color == "random"){
      colors <- randomColor(ntriangles, hue = hue, luminosity = luminosity)
    }else{
      colors <- distinctColorPalette(ntriangles)
    }
    for(i in 1L:ntriangles){
      triangle <- vertices[trianglesIds[i, ], ]
      polygon(triangle, border = NA, col = colors[i])
    }
  }
  if(!is.null(border)){
    edges <- tessellation[["edges"]][, c(1L, 2L)]
    for(i in 1L:nrow(edges)){
      edge <- edges[i, ]
      p0 <- vertices[edge[1L], ]
      p1 <- vertices[edge[2L], ]
      segments(
        p0[1L], p0[2L], p1[1L], p1[2L], col = border, lty = lty, lwd = lwd
      )
    }
  }
}


#' @title Plot 3D Delaunay tessellation
#' @description Plot a 3D Delaunay tessellation with \strong{rgl}.
#'
#' @param tessellation the output of \code{\link{delaunay}} with 3D points
#' @param color controls the filling colors of the tetrahedra, either
#'   \code{FALSE} for no color, \code{"random"} to use
#'   \code{\link[randomcoloR]{randomColor}}, or \code{"distinct"} to use
#'   \code{\link[randomcoloR]{distinctColorPalette}}
#' @param hue,luminosity if \code{color="random"}, these arguments are passed
#'   to \code{\link[randomcoloR]{randomColor}}
#' @param alpha opacity, number between 0 and 1
#'
#' @return No value, just renders a 3D plot.
#' @export
#' @importFrom randomcoloR randomColor distinctColorPalette
#' @importFrom utils combn
#' @importFrom rgl triangles3d
#'
#' @examples library(RCGAL)
#' pts <- rbind(
#'   c(-5, -5,  16),
#'   c(-5,  8,   3),
#'   c(4,  -1,   3),
#'   c(4,  -5,   7),
#'   c(4,  -1, -10),
#'   c(4,  -5, -10),
#'   c(-5,  8, -10),
#'   c(-5, -5, -10)
#' )
#' tess <- delaunay(pts)
#' library(rgl)
#' open3d(windowRect = c(50, 50, 562, 562))
#' plotDelaunay3D(tess)
plotDelaunay3D <- function(
  tessellation, color = "distinct", hue = "random", luminosity = "light",
  alpha = 0.3
){
  if(!inherits(tessellation, "delaunay")){
    stop(
      "The argument `tessellation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  vertices <- attr(tessellation, "points")
  if(ncol(vertices) != 3L){
    stop(
      sprintf("Invalid dimension (%d instead of 3).", ncol(vertices)),
      call. = TRUE
    )
  }
  cells <- tessellation[["cells"]]
  ntetrahedra <- length(cells)
  if(!isFALSE(color)){
    color <- match.arg(color, c("random", "distinct"))
    if(color == "random"){
      colors <- randomColor(ntetrahedra, hue = hue, luminosity = luminosity)
    }else{
      colors <- distinctColorPalette(ntetrahedra)
    }
    triangles <- combn(4L, 3L)
    for(i in 1L:ntetrahedra){
      cellIds <- cells[[i]][["cell"]]
      simplex <- vertices[cellIds, ]
      for(j in 1L:4L){
        triangles3d(simplex[triangles[, j], ], color = colors[i], alpha = alpha)
      }
    }
  }
  edges <- tessellation[["edges"]]
  for(i in 1L:nrow(edges)){
    edge <- edges[i, ]
    p1 <- vertices[edge[1L], ]
    p2 <- vertices[edge[2L], ]
    lines3d(rbind(p1, p2), color = "black")
  }
  invisible(NULL)
}
