delaunayArea <- function(vertices, triangles){
  ntriangles <- nrow(triangles)
  areas <- numeric(ntriangles)
  for(i in 1L:ntriangles){
    points <- vertices[triangles[i, ], ]
    areas[i] <- triangleArea(points[1L, ], points[2L, ], points[3L, ])
  }
  sum(areas)
}

#' @title Delaunay tessellation
#' @description Delaunay tessellation of a set of 2D or 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param elevation if points are three-dimensional and \code{elevation=TRUE},
#'   then the function performs an elevated two-dimensional Delaunay
#'   triangulation, using the \code{z} coordinates of the points for the
#'   elevations; see the example
#' @param constraints \emph{for 2D only}, some edges to perform a constrained
#'   Delaunay triangulation, given as an integer matrix with two columns (each
#'   row provides the indices of the two points forming the edge);
#'   \code{NULL} for no constraint
#'
#' @return The Delaunay tessellation.
#' \itemize{
#'   \item \strong{If the dimension is 2} and \code{constraints=NULL},
#'         the returned value is a list with three fields:
#'         \code{faces}, \code{edges} and \code{area}. The \code{faces} field
#'         contains an integer matrix with three columns; each row represents a
#'         triangle whose each vertex is given by the index (row number) of
#'         this point in the \code{points} matrix. The \code{edges} field
#'         also contains an integer matrix with three columns. The two first
#'         integers of a row are the indices of the two points which form the
#'         edge. The third column, named \code{border}, only contains some
#'         zeros and some ones; a border (exterior) edge is labelled by a
#'         \code{1}. The \code{area} field contains only a number: the area
#'         of the triangulated region (that is, the area of the convex hull of
#'         the points).
#'   \item \strong{If the dimension is 2} and \code{constraints} is not
#'         \code{NULL}, the returned value is a list with
#'         four fields: \code{faces}, \code{edges}, \code{constraints}, and
#'         \code{area}. The \code{faces} field
#'         contains an integer matrix with three columns; each row represents a
#'         triangle whose each vertex is given by the index (row number) of
#'         this point in the \code{points} matrix. The \code{edges} field
#'         is a dataframe with four columns. The first two columns provide
#'         the edges of the triangulation; they are given by row, the two
#'         integers of a row are the indices of the two points which form the
#'         edge. Each integer of the third column is the index of the face
#'         the corresponding edge belongs to. The fourth column,
#'         named \code{border}, only contains some
#'         zeros and some ones; a border edge is labelled by a
#'         \code{1}.
#'         The \code{constraints} field is an integer matrix with
#'         two columns, it represents the constraint edges.
#'         Finally, the \code{area} field contains only a number: the area
#'         of the triangulated region.
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
#'         five fields: \code{mesh}, \code{edges}, \code{faceVolumes},
#'         \code{volume} and \code{area}. The \code{mesh} field is an object of
#'         class \code{mesh3d}, ready for plotting with the \strong{rgl}
#'         package. The \code{edges} field provides the indices of the edges,
#'         given as an integer matrix with two columns. The \code{faceVolumes}
#'         field is a numeric vector, it provides the volumes under the faces
#'         that can be found in the \code{mesh} field. The \code{volume} field
#'         provides the sum of these volumes, that is to say the total volume
#'         under the triangulated surface. Finally, the \code{area} field
#'         provides the sum of the areas of all triangles, thereby
#'         approximating the area of the triangulated surface.
#' }
#' @export
#' @importFrom rgl tmesh3d addNormals
#' @importFrom Rvcg vcgGetEdge
#'
#' @examples library(RCGAL)
#' # elevated Delaunay triangulation ####
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
  points, elevation = FALSE, constraints = NULL
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
  if(!is.null(constraints)){
    if(dimension == 3L){
      stop(
        "If you set some constraints, you must provide two-dimensional points.",
        call. = TRUE
      )
    }
    if(
      !is.matrix(constraints) || !is.numeric(constraints) ||
      ncol(constraints) != 2L
    ){
      stop(
        "The `constraints` argument must be an integer matrix with two columns.",
        call. = TRUE
      )
    }
    if(any(is.na(constraints))){
      stop("Missing values in `constraints` are not allowed.", call. = TRUE)
    }
    storage.mode(constraints) <- "integer"
    stopifnot(all(constraints >= 1L))
    stopifnot(all(constraints <= nrow(points)))
    constraints <- t(apply(constraints, 1L, sort))
    if(anyDuplicated(constraints)){
      stop("There are some duplicated constraints.", call. = TRUE)
    }
    if(any(constraints[, 1L] == constraints[, 2L])){
      stop("There are some invalid constraints.", call. = TRUE)
    }
    triangles <- del2d_constrained_cpp(points, constraints)
    # edges <- apply(triangles, 1L, function(x){
    #   rbind(c(x[1L], x[2L]), c(x[1L], x[3L]), c(x[2L], x[3L]))
    # }, simplify = FALSE)
    # edges <- do.call(rbind, edges)
    # edges <- edges[!duplicated(edges), ]
    mesh <- tmesh3d(
      vertices = t(cbind(points, 0)),
      indices = t(triangles)
    )
    out <- list(
      "faces"       = triangles,
      "edges"       = `colnames<-`(
        as.matrix(vcgGetEdge(mesh))[, c(1L, 2L, 4L)], c("v1", "v2", "border")
      ),
      "constraints" = constraints,
      "area"        = delaunayArea(points, triangles)
    )
    attr(out, "constrained") <- TRUE
  }else if(dimension == 2L && is.null(constraints)){
    out <- del2d_cpp(points)
    out[["area"]] <- delaunayArea(points, out[["faces"]])
  }else{
    if(elevation){
      del <- del2d_xy_cpp(points)
      triangles <- del[["faces"]]
      ntriangles <- nrow(triangles)
      areas <- numeric(ntriangles)
      for(i in 1L:ntriangles){
        triangle <- triangles[i, ]
        vertices <- points[triangle, ]
        areas[i] <- triangleArea(
          vertices[1L, ], vertices[2L, ], vertices[3L, ]
        )
      }
      out <- list(
        "mesh"        = addNormals(
          tmesh3d(
            vertices = t(points),
            indices = t(triangles)
          )
        ),
        "edges"       = del[["edges"]],
        "faceVolumes" = attr(triangles, "volumes"),
        "volume"      = del[["volume"]],
        "area"        = sum(areas)
      )
    }else{
      out <- del3d_cpp(points)
    }
  }
  class(out) <- "delaunay"
  attr(out, "points") <- points
  out
}



#' @title Plot 2D Delaunay triangulation
#' @description Plot a constrained or unconstrained 2D Delaunay triangulation.
#'
#' @param triangulation an output of \code{\link{delaunay}} without constraints
#'   (\code{constraints=NULL}) or with constraints
#' @param col_edges the color of the edges of the triangles which are not
#'   border edges nor constraint edges; \code{NULL} for no color
#' @param col_borders the color of the border edges; note that the border
#'   edges can contain the constraint edges for a constrained
#'   Delaunay tessellation; \code{NULL} for no color
#' @param col_constraints for a constrained Delaunay tessellation, the color
#'   of the constraint edges which are not border edges;
#'   \code{NULL} for no color
#' @param fillcolor controls the filling colors of the triangles, either
#'   \code{NULL} for no color, a single color, \code{"random"} to get multiple
#'   colors with \code{\link[randomcoloR]{randomColor}}, or \code{"distinct"}
#'   get multiple colors with \code{\link[randomcoloR]{distinctColorPalette}}
#' @param hue,luminosity if \code{color = "random"}, these arguments are passed
#'   to \code{\link[randomcoloR]{randomColor}}
#' @param lty_edges,lwd_edges graphical parameters for the edges which are not
#'   border edges nor constraint edges
#' @param lty_borders,lwd_borders graphical parameters for the border edges
#' @param lty_constraints,lwd_constraints in the case of a constrained Delaunay
#'   triangulation, graphical parameters for the constraintt edges which are
#'   not border edges
#' @param ... arguments passed to \code{\link{plot}} for the vertices, such as
#'   \code{type="n"} or \code{asp=1}
#'
#' @return No value, just renders a 2D plot.
#'
#' @seealso \code{\link{mesh2d}} for an interactive plot
#'
#' @export
#' @importFrom randomcoloR randomColor distinctColorPalette
#' @importFrom graphics plot polygon par segments
#' @importFrom gplots col2hex
#'
#' @examples library(RCGAL)
#' # random points in a square ####
#' square <- rbind(
#'   c(-1, 1), c(1, 1), c(1, -1), c(-1, -1)
#' )
#' library(uniformly)
#' set.seed(314)
#' ptsinsquare <- runif_in_cube(10L, d = 2L)
#' pts <- rbind(square, ptsinsquare)
#' d <- delaunay(pts)
#' opar <- par(mar = c(0, 0, 0, 0))
#' plotDelaunay2D(
#'   d, type = "n", xlab = NA, ylab = NA, axes = FALSE, asp = 1,
#'   fillcolor = "random", luminosity = "dark", lwd_borders = 3
#' )
#' par(opar)
#'
#' # a constrained Delaunay triangulation: outer and inner hexagons ####
#' nsides <- 6L
#' angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
#' outer_points <- cbind(cos(angles), sin(angles))
#' inner_points <- outer_points / 2
#' points <- rbind(outer_points, inner_points)
#' # constraint edges
#' indices <- 1L:nsides
#' edges <- cbind(
#'   indices, c(indices[-1L], indices[1L])
#' )
#' edges <- rbind(edges, edges + nsides)
#' # constrained Delaunay triangulation
#' d <- delaunay(points, constraints = edges)
#' opar <- par(mar = c(0, 0, 0, 0))
#' plotDelaunay2D(
#'   d, type = "p", pch = 19, xlab = NA, ylab = NA, axes = FALSE, asp = 1,
#'   fillcolor = "orange", lwd_borders = 3
#' )
#' par(opar)
#'
#' # another constrained Delaunay tesselation: a face ####
#' V <- as.matrix(read.table(
#'   system.file("extdata", "face_vertices.txt", package = "RCGAL")
#' ))[, c(2L, 3L)]
#' E <- as.matrix(read.table(
#'   system.file("extdata", "face_edges.txt", package = "RCGAL")
#' ))[, c(2L, 3L)]
#' d <- delaunay(points = V, constraints = E)
#' opar <- par(mar = c(0, 0, 0, 0))
#' plotDelaunay2D(
#'   d, type = "n", xlab = NA, ylab = NA, axes = FALSE, asp = 1,
#'   fillcolor = "salmon", col_borders = "black",
#'   lwd_borders = 3, lwd_constraints = 2, lty_edges = "dashed"
#' )
#' par(opar)
plotDelaunay2D <- function(
  triangulation,
  col_edges = "black", col_borders = "red", col_constraints = "green",
  fillcolor = "distinct", hue = "random", luminosity = "light",
  lty_edges = par("lty"), lwd_edges = par("lwd"),
  lty_borders = par("lty"), lwd_borders = par("lwd"),
  lty_constraints = par("lty"), lwd_constraints = par("lwd"), ...
){
  if(!inherits(triangulation, "delaunay")){
    stop(
      "The argument `triangulation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  vertices <- attr(triangulation, "points")
  if(ncol(vertices) != 2L){
    stop(
      sprintf("Invalid dimension (%d instead of 2).", ncol(vertices)),
      call. = TRUE
    )
  }
  plot(vertices, ...)
  if(!isFalsy(fillcolor)){
    fillcolor <- tryCatch({
      col2hex(fillcolor)
    }, error = function(e){
      match.arg(fillcolor, c("random", "distinct"))
    })
    triangles <- triangulation[["faces"]]
    ntriangles <- nrow(triangles)
    if(fillcolor == "random"){
      colors <- randomColor(ntriangles, hue = hue, luminosity = luminosity)
    }else if(fillcolor == "distinct"){
      colors <- distinctColorPalette(ntriangles)
    }else{
      colors <- rep(fillcolor, ntriangles)
    }
    for(i in 1L:ntriangles){
      triangle <- makeTriangle(vertices, triangles[i, ])
      polygon(triangle, border = NA, col = colors[i])
    }
  }
  constraintEdges <- triangulation[["constraints"]]
  allEdges <- triangulation[["edges"]]
  borderEdges <- allEdges[allEdges[, "border"] == 1L, c(1L, 2L)]
  allEdges <- allEdges[, c(1L, 2L)]
  specialEdges <- unionEdges(borderEdges, constraintEdges)
  constraintEdges <- subtractEdges(specialEdges, borderEdges)
  otherEdges <- subtractEdges(allEdges, specialEdges)
  if(!isFalsy(col_edges)){
    # if(is.null(constraintEdges)){
    #   edges <- otherEdges
    # }else{
    #   if(!isFalsy(col_constraints)){
    #     edges <- subtractEdges()
    #   }
    # }
    edges <- otherEdges
    for(i in 1L:nrow(edges)){
      edge <- edges[i, ]
      p0 <- vertices[edge[1L], ]
      p1 <- vertices[edge[2L], ]
      segments(
        p0[1L], p0[2L], p1[1L], p1[2L], col = col_edges,
        lty = lty_edges, lwd = lwd_edges
      )
    }
  }
  if(!isFalsy(col_borders)){
    edges <- borderEdges
    for(i in 1L:nrow(edges)){
      edge <- edges[i, ]
      p0 <- vertices[edge[1L], ]
      p1 <- vertices[edge[2L], ]
      segments(
        p0[1L], p0[2L], p1[1L], p1[2L], col = col_borders,
        lty = lty_borders, lwd = lwd_borders
      )
    }
  }
  if(!is.null(constraintEdges) && !isFalsy(col_constraints)){
    edges <- constraintEdges
    for(i in 1L:nrow(edges)){
      edge <- edges[i, ]
      p0 <- vertices[edge[1L], ]
      p1 <- vertices[edge[2L], ]
      segments(
        p0[1L], p0[2L], p1[1L], p1[2L], col = col_constraints,
        lty = lty_constraints, lwd = lwd_constraints
      )
    }
  }
  invisible(NULL)
}

#' @title Convert a 2D Delaunay triangulation to a 'rgl' mesh
#' @description Makes a 'rgl' mesh (\code{\link[rgl]{mesh3d}} object) from
#'   a 2D Delaunay triangulation, unconstrained or constrained.
#'
#' @param triangulation an output of \code{\link{delaunay}} executed with
#'   2D points
#'
#' @return A list with three fields; \code{mesh}, a \code{\link[rgl]{mesh3d}}
#'   object, \code{borderEdges}, a numeric matrix that can be used with
#'   \code{\link[rgl]{segments3d}} to plot the border edges, and
#'   \code{constraintEdges}, a numeric matrix that can be used with
#'   \code{\link[rgl]{segments3d}} to plot the constraint edges which are
#'   not border edges.
#'
#' @seealso \code{\link{plotDelaunay2D}}
#'
#' @export
#' @importFrom rgl tmesh3d
#'
#' @examples library(RCGAL)
#' # outer and inner hexagons ####
#' nsides <- 6L
#' angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
#' outer_points <- cbind(cos(angles), sin(angles))
#' inner_points <- outer_points / 2
#' points <- rbind(outer_points, inner_points)
#' # constraint edges
#' indices <- 1L:nsides
#' edges <- cbind(
#'   indices, c(indices[-1L], indices[1L])
#' )
#' edges <- rbind(edges, edges + nsides)
#' # constrained Delaunay triangulation
#' del <- delaunay(points, constraints = edges)
#' # mesh
#' m2d <- mesh2d(del)
#' mesh <- m2d[["mesh"]]
#' # plot all edges with `wire3d`
#' library(rgl)
#' open3d(windowRect = c(100, 100, 612, 612))
#' shade3d(mesh, color = "red", specular = "orangered")
#' wire3d(mesh, color = "black", lwd = 3, specular = "black")
#' # plot only the border edges
#' open3d(windowRect = c(100, 100, 612, 612))
#' shade3d(mesh, color = "darkred", specular = "firebrick")
#' segments3d(m2d[["borderEdges"]], lwd = 3)
mesh2d <- function(triangulation){
  if(!inherits(triangulation, "delaunay")){
    stop(
      "The argument `triangulation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  vertices <- attr(triangulation, "points")
  if(ncol(vertices) != 2L){
    stop(
      sprintf("Invalid dimension (%d instead of 2).", ncol(vertices)),
      call. = TRUE
    )
  }
  vertices <- cbind(vertices, 0)
  mesh <- tmesh3d(
    vertices = t(vertices),
    indices = t(triangulation[["faces"]])
  )
  constraintEdges <- triangulation[["constraints"]]
  allEdges <- triangulation[["edges"]]
  borderEdges <- allEdges[allEdges[, "border"] == 1L, c(1L, 2L)]
  constraints <- NULL
  if(!is.null(constraintEdges)){
    specialEdges <- unionEdges(borderEdges, constraintEdges)
    constraintEdges <- subtractEdges(specialEdges, borderEdges)
    if(!is.null(constraintEdges)){
      constraints <- do.call(rbind, apply(
        constraintEdges, 1L, function(ij) vertices[ij,], simplify = FALSE
      ))
    }
  }
  borders <- do.call(rbind, apply(
    borderEdges, 1L, function(ij) vertices[ij,], simplify = FALSE
  ))
  list(
    "mesh"            = mesh,
    "borderEdges"     = borders,
    "constraintEdges" = constraints
  )
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
