makeFaceFamilies <- function(faces, edges){
  Faces <- apply(faces, 1L, function(f){
    f <- sort(f)
    c(
      paste0(f[c(1L, 2L)], collapse = "-"),
      paste0(f[c(1L, 3L)], collapse = "-"),
      paste0(f[c(2L, 3L)], collapse = "-")
    )
  })
  Edges0 <- apply(edges[edges[, 3L]==0L, c(1L, 2L)], 1L, function(e){
    paste0(sort(e), collapse = "-")
  })
  Edges0Faces <- apply(Faces, 2L, function(f){
    sort(match(f, Edges0))
  }, simplify = FALSE)
  nfaces <- nrow(faces)
  families <- vector("list", nfaces)
  fseq <- 1L:nfaces
  for(j in fseq){
    Edges0Faces_j <- Edges0Faces[[j]]
    for(i in fseq[-j]){
      if(any(Edges0Faces_j %in% Edges0Faces[[i]])){
        families[[j]] <- families[[i]] <- union(Edges0Faces[[i]], Edges0Faces_j)
        break
      }
    }
  }
  vapply(families, function(f){
    if(is.null(f)) NA_character_ else paste0(sort(f), collapse = "-")
  }, character(1L))
}

#' @title Convex hull
#' @description Convex hull of a set of 2D or 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param faceFamilies Boolean, for 3D only; faces are always triangular, and
#'   the family of a face is the set of faces adjacent and coplanar with this
#'   face (in other words, the coplanar neighbors of this face); set this
#'   argument to \code{TRUE} if you want the face families; this gives a label
#'   of the family for each face, or \code{NA} is the face is "alone" (has no
#'   coplanar neighbor); this is useful for plotting (see the examples in
#'   \code{\link{plotConvexHull3D}})
#' @param epsilon for 3D only, a small nonnegative number; this number plays
#'   a role in the detection of border edges: an edge is considered as a
#'   non-border edge when there is approximate equality between the unit
#'   normals of the two adjacent faces of this edge, and \code{epsilon}
#'   defines the degree of the approximation (perfect equality corresponds to
#'   \code{epsilon=0})
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
#'                of the faces, \emph{circumcenters}, which provides the
#'                circumcenters of the faces, and \emph{families} if you set
#'                \code{faceFamilies=TRUE}.}
#'           \item{\emph{surface}}{A number, the surface of the convex hull.}
#'           \item{\emph{volume}}{A number, the volume of the convex hull.}
#'         }
#' }
#' @export
#'
#' @examples library(RCGAL)
#' # 2D example ####
#' pts <- rbind(
#'   c(-1, -1),
#'   c(-1,  1),
#'   c( 1, -1),
#'   c( 1,  1),
#'   c( 2,  0),
#'   c( 0,  2),
#'   c(-2,  0),
#'   c( 0, -2)
#' )
#' hull <- convexhull(pts)
#' # it's easy to plot a 2D convex hull:
#' plot(hull[["vertices"]], asp = 1, pch = 19)
#' polygon(hull[["vertices"]], col = "green")
#'
#' # a 3D example ####
#' cube <- rbind(
#'   c(-1, -1, -1),
#'   c(-1, -1,  1),
#'   c(-1,  1, -1),
#'   c(-1,  1,  1),
#'   c( 1, -1, -1),
#'   c( 1, -1,  1),
#'   c( 1,  1, -1),
#'   c( 1,  1,  1),
#'   c( 0,  0,  0)
#' )
#' hull <- convexhull(cube)
#' hull[["vertices"]][[1]]
#' # the non-border edges are the diagonals of the faces:
#' hull[["edges"]]
#' hull[["surface"]]
#' hull[["volume"]]
#' # plot:
#' library(rgl)
#' open3d(windowRect = c(50, 50, 562, 562))
#' plotConvexHull3D(hull)
convexhull <- function(
  points, faceFamilies = FALSE, epsilon = sqrt(.Machine$double.eps)
){
  stopifnot(isNonNegativeNumber(epsilon))
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
    hull <- cxhull2d_cpp(points)
  }else{
    hull <- cxhull3d_cpp(points)
    if(faceFamilies){
      attr(hull[["faces"]], "families") <-
        makeFaceFamilies(hull[["faces"]], hull[["edges"]])
    }
  }
  class(hull) <- "cxhull"
  attr(hull, "points") <- points
  hull
}


#' @title Plot 3D convex hull
#' @description Plot a 3D convex hull with \strong{rgl}.
#'
#' @param hull the output of \code{\link{convexhull}} with 3D points
#' @param color controls the colors of the faces, either
#'   \code{FALSE} for no color, \code{"random"} to use
#'   \code{\link[randomcoloR]{randomColor}}, \code{"distinct"} to use
#'   \code{\link[randomcoloR]{distinctColorPalette}}, or a single color
#' @param hue,luminosity if \code{color="random"}, these arguments are passed
#'   to \code{\link[randomcoloR]{randomColor}}
#' @param alpha opacity, number between 0 and 1
#' @param edgesAsTubes Boolean, whether to plot the edges as tubes
#' @param tubeRadius if \code{edgesAsTubes=TRUE}, the radius of the tubes
#' @param tubeColor if \code{edgesAsTubes=TRUE}, the color of the tubes
#'
#' @return No value, just renders a 3D plot.
#' @export
#' @importFrom randomcoloR randomColor distinctColorPalette
#' @importFrom rgl triangles3d spheres3d cylinder3d shade3d lines3d
#' @examples library(RCGAL)
#' library(rgl)
#' # blue dodecahedron with edges as tubes ####
#' dodecahedron <- t(dodecahedron3d()$vb[-4, ])
#' hull <- convexhull(dodecahedron)
#' open3d(windowRect = c(50, 50, 562, 562))
#' plotConvexHull3D(
#'   hull, color = "navy", edgesAsTubes = TRUE,
#'   tubeRadius = 0.03, tubeColor = "gold"
#' )
#'
#' # the dodecahedron with multiple colors ####
#' hull <- convexhull(dodecahedron, faceFamilies = TRUE)
#' plotConvexHull3D(hull, color = "random", luminosity = "bright")
#'
#' # a strange convex hull ####
#' pt <- function(x){
#'   c(
#'     sin(x) * cos(2 * x),
#'     sin(x) * sin(2 * x),
#'     cos(x)
#'   )
#' }
#' pts <- t(vapply(seq(0, pi, length.out = 50), pt, numeric(3L)))
#' hull <- convexhull(pts)
#' plotConvexHull3D(hull, color = "random", hue = "purple", luminosity = "dark")
plotConvexHull3D <- function(
  hull, color = "distinct", hue = "random", luminosity = "light",
  alpha = 1, edgesAsTubes = FALSE, tubeRadius, tubeColor
){
  if(!inherits(hull, "cxhull")){
    stop(
      "The argument `hull` must be an output of the `convexhull` function.",
      call. = TRUE
    )
  }
  vertices <- attr(hull, "points")
  if(ncol(vertices) != 3L){
    stop(
      sprintf("Invalid dimension (%d instead of 3).", ncol(vertices)),
      call. = TRUE
    )
  }
  triangles <- hull[["faces"]]
  ntriangles <- ncolors <- nrow(triangles)
  faceFamilies <- FALSE
  if(!isFALSE(color)){
    Color <- pmatch(color, c("random", "distinct"))
    if(is.na(Color)){
      colors <- rep(color, ntriangles)
      for(i in 1L:ntriangles){
        triangles3d(vertices[triangles[i, ], ], color = colors[i])
      }
    }else{
      if(!is.null(families <- attr(triangles, "families"))){
        faceFamilies <- TRUE
        NAfamilies <- which(is.na(families))
        families[NAfamilies] <- paste0("NA", NAfamilies)
        distinctFamilies <- unique(families)
        ncolors <- length(distinctFamilies)
      }
      if(Color == 1L){
        colors <- randomColor(ncolors, hue = hue, luminosity = luminosity)
      }else{
        colors <- distinctColorPalette(ncolors)
      }
      if(faceFamilies){
        names(colors) <- distinctFamilies
        for(i in 1L:ntriangles){
          triangles3d(vertices[triangles[i, ], ], color = colors[families[i]])
        }
      }else{
        for(i in 1L:ntriangles){
          triangles3d(vertices[triangles[i, ], ], color = colors[i])
        }
      }
    }
  }
  edges <- hull[["edges"]]
  isborder <- edges[, 3L] == 1L
  edges <- edges[isborder, c(1L, 2L)]
  nedges <- nrow(edges)
  if(edgesAsTubes){
    for(i in 1L:nedges){
      shade3d(
        cylinder3d(vertices[edges[i, ], ], radius = tubeRadius, sides = 90),
        color = tubeColor
      )
    }
    vertices <- do.call(rbind, lapply(hull[["vertices"]], `[[`, "point"))
    spheres3d(
      vertices, radius = 1.5*tubeRadius, color = tubeColor
    )
  }else{
    for(i in 1L:nedges){
      edge <- edges[i, ]
      lines3d(vertices[edge, ])
    }
  }
  invisible(NULL)
}
