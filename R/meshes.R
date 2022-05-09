#' @importFrom data.table uniqueN
#' @noRd
checkMesh <- function(vertices, faces){
  if(!is.matrix(vertices) || ncol(vertices) != 3L){
    stop("The `vertices` argument must be a matrix with three columns.")
  }
  storage.mode(vertices) <- "double"
  homogeneousFaces <- FALSE
  isTriangle <- FALSE
  if(is.matrix(faces)){
    if(ncol(faces) < 3L){
      stop("Faces must be given by at least three indices.")
    }
    storage.mode(faces) <- "integer"
    if(anyNA(faces)){
      stop("Found missing values in `faces`.")
    }
    if(any(faces < 1L)){
      stop("Faces cannot contain indices lower than 1.")
    }
    if(any(faces > nrow(vertices))){
      stop("Faces cannot contain indices higher than the number of vertices.")
    }
    homogeneousFaces <- TRUE
    isTriangle <- ncol(faces) == 3L
    faces <- lapply(1L:nrow(faces), function(i) faces[i, ] - 1L)
  }else if(is.list(faces)){
    check <- all(vapply(faces, isAtomicVector, logical(1L)))
    if(!check){
      stop("The `faces` argument must be a list of integer vectors.")
    }
    check <- any(vapply(faces, anyNA, logical(1L)))
    if(check){
      stop("Found missing values in `faces`.")
    }
    faces <- lapply(faces, function(x) as.integer(x) - 1L)
    sizes <- lengths(faces)
    if(any(sizes < 3L)){
      stop("Faces must be given by at least three indices.")
    }
    check <- any(vapply(faces, function(f){
      any(f < 0L) || any(f >= nrow(vertices))
    }, logical(1L)))
    if(check){
      stop(
        "Faces cannot contain indices lower than 1 or higher than the ",
        "number of vertices."
      )
    }
    homogeneousFaces <- uniqueN(sizes) == 1L
    isTriangle <- homogeneousFaces && sizes[1L] == 3L
  }else{
    stop("The `faces` argument must be a list or a matrix.")
  }
  list(
    vertices = t(vertices),
    faces = faces,
    homogeneousFaces = homogeneousFaces,
    isTriangle = isTriangle
  )
}

#' @title Make a 3D mesh
#' @description Make a 3D mesh from given vertices and faces; the returned
#'   faces are coherently oriented and normals are computed if desired.
#'
#' @param vertices a numeric matrix with three columns
#' @param faces either an integer matrix (each row provides the vertex indices
#'   of the corresponding face) or a list of integer vectors, each one
#'   providing the vertex indices of the corresponding face
#' @param triangulate Boolean, whether to triangulate the faces
#' @param merge Boolean, whether to merge duplicated vertices
#' @param normals Boolean, whether to compute the normals
#' @param epsilon if the mesh is triangulated or if \code{triangulate=TRUE},
#'   then \code{epsilon} is used in the detection of exterior edges (see the
#'   \strong{Value} section); the higher value of \code{epsilon}, the lower
#'   number of exterior edges
#'
#' @return A list giving the vertices, the edges, the faces of the mesh, and
#'   optionally the normals. This list has an additional component
#'   \code{edges0} if \code{triangulate=TRUE}, giving the edges before the
#'   triangulation, unless the mesh is already triangulated, in which case
#'   the \code{triangulate} option is ignored. If the mesh is already
#'   triangulated or if \code{triangulate=TRUE}, the output list has an
#'   additional component \code{exteriorEdges}, giving the exterior edges
#'   of the mesh.
#'
#' @export
#'
#' @examples
#' library(RCGAL)
#' library(rgl)
#'
#' # a tetrahedron with ill-oriented faces ####
#' vertices <- rbind(
#'   c(-1, -1, -1),
#'   c(1, 1, -1),
#'   c(1, -1, 1),
#'   c(-1, 1, 1)
#' )
#' faces <- rbind(
#'   c(1, 2, 3),
#'   c(3, 4, 2),
#'   c(4, 2, 1),
#'   c(4, 3, 1)
#' )
#'
#' # plot the tetrahedron, hiding the back of the faces
#' # then some faces do not appear, as their orientation is not correct
#' tmesh1 <- tmesh3d(
#'   vertices = t(vertices),
#'   indices = t(faces),
#'   homogeneous = FALSE
#' )
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(tmesh1, color = "green", back = "cull")
#'
#' # now run the `Mesh` function
#' mesh2 <- Mesh(vertices, faces, normals = FALSE)
#' # plot the tetrahedron, hiding the back of the faces
#' # then all faces appear now
#' tmesh2 <- tmesh3d(
#'   vertices = t(mesh2[["vertices"]]),
#'   indices = t(mesh2[["faces"]]),
#'   homogeneous = FALSE
#' )
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(tmesh2, color = "blue", back = "cull")
#'
#' # illustration of the `merge` option ####
#' # we construct a mesh with a lot of duplicated vertices
#' library(misc3d) # to compute a mesh of an isosurface
#' a <- 0.94; mu <- 0.56; c <- 0.34 # cyclide parameters
#' f <- function(x, y, z, a, c, mu){ # implicit equation of the cyclide
#'   b <- sqrt(a^2 - c^2)
#'   (x^2 + y^2 + z^2 - mu^2 + b^2)^2 - 4*(a*x - c*mu)^2 - 4*b^2*y^2
#' }
#' x <- seq(-c - mu - a, abs(mu - c) + a, length.out = 45)
#' y <- seq(-mu - a, mu + a, length.out = 45)
#' z <- seq(-mu - c, mu + c, length.out = 30)
#' g <- expand.grid(x = x, y = y, z = z)
#' voxel <- array(with(g, f(x, y, z, a, c, mu)), c(45, 45, 30))
#' cont <- computeContour3d(voxel, level = 0, x = x, y = y, z = z)
#' ids <- matrix(1:nrow(cont), ncol = 3, byrow = TRUE)
#' # run the `Mesh` function with `merge=TRUE`
#' mesh <- Mesh(cont, ids, merge = TRUE)
#' # plot the cyclide
#' tmesh <- tmesh3d(
#'   vertices = t(mesh[["vertices"]]),
#'   indices = t(mesh[["faces"]]),
#'   normals = mesh[["normals"]],
#'   homogeneous = FALSE
#' )
#' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
#' shade3d(tmesh, color = "green")
#'
#' # illustration of the `triangulate` option ####
#' # the faces of the truncated icosahedron are hexagonal or pentagonal:
#' truncatedIcosahedron[["faces"]]
#' # so we triangulate them:
#' mesh <- Mesh(
#'   truncatedIcosahedron[["vertices"]],
#'   truncatedIcosahedron[["faces"]],
#'   triangulate = TRUE, normals = FALSE
#' )
#' # now we can plot the truncated icosahedron
#' tmesh <- tmesh3d(
#'   vertices = t(mesh[["vertices"]]),
#'   indices = t(mesh[["faces"]]),
#'   homogeneous = FALSE
#' )
#' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
#' shade3d(tmesh, color = "orange")
Mesh <- function(
  vertices, faces, triangulate = FALSE, merge = FALSE, normals = TRUE,
  epsilon = 0
){
  stopifnot(epsilon >= 0)
  checkedMesh <- checkMesh(vertices, faces)
  vertices <- checkedMesh[["vertices"]]
  faces <- checkedMesh[["faces"]]
  homogeneousFaces <- checkedMesh[["homogeneousFaces"]]
  isTriangle <- checkedMesh[["isTriangle"]]
  mesh <- SurfMesh(
    vertices, faces, isTriangle, triangulate, merge, normals, epsilon
  )
  if(triangulate && isTriangle){
    message(
      "Ignored option `triangulate`, since the mesh is already triangulated."
    )
    triangulate <- FALSE
  }
  # mesh <- if(normals){
  #   SurfMeshWithNormals(t(vertices), faces, merge)
  # }else{
  #   if(triangulate){
  #     SurfTMesh(t(vertices), faces, merge)
  #   }else{
  #     SurfMesh(t(vertices), faces, merge)
  #   }
  # }
  mesh[["vertices"]] <- t(mesh[["vertices"]])
  edges <- unname(t(mesh[["edges"]]))
  if(triangulate || isTriangle){
    mesh[["exteriorEdges"]] <- edges[edges[, 3L] == 1L, c(1L, 2L)]
    mesh[["edges"]] <- edges[, c(1L, 2L)]
  }
  if(triangulate){
    mesh[["edges0"]] <- t(mesh[["edges0"]])
  }
  if(normals){
    mesh[["normals"]] <- t(mesh[["normals"]])
  }
  if(triangulate || homogeneousFaces){
    mesh[["faces"]] <- do.call(rbind, mesh[["faces"]])
  }
  mesh
}

#' @title Meshes intersection
#' @description Computes the intersection of the given meshes.
#'
#' @param meshes a list of \emph{triangular} meshes, each given as a list with
#'   (at least) two fields: \code{vertices} and \code{faces}
#' @param merge Boolean, whether to merge the duplicated vertices of the
#'   input meshes and the output mesh
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#' @param exact Boolean, whether to use exact calculations; this is slower but
#'   more accurate
#'
#' @return A triangular mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges} and \code{normals}
#'   if \code{normals=TRUE}.
#'
#' @importFrom rgl tmesh3d
#' @importFrom Rvcg vcgUpdateNormals
#'
#' @export
#'
#' @examples
#' library(RCGAL)
#' library(rgl)
#'
#' # mesh one: truncated icosahedron; one has to triangulate it
#' mesh1 <- Mesh(
#'   truncatedIcosahedron[["vertices"]],
#'   truncatedIcosahedron[["faces"]],
#'   triangulate = TRUE, normals = FALSE
#' )
#'
#' # mesh two: a cube; one also has to triangulate it
#' cube <- translate3d( # (from the rgl package)
#'   cube3d(), 2, 0, 0
#' )
#' vertices <- t(cube$vb[-4L, ])
#' faces <- t(cube$ib)
#' mesh2 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
#'
#' # compute the intersection
#' inter <- MeshesIntersection(list(mesh1, mesh2))
#'
#' # plot
#' rglmesh1 <- tmesh3d(
#'   vertices = t(mesh1[["vertices"]]),
#'   indices = t(mesh1[["faces"]]),
#'   homogeneous = FALSE
#' )
#' rglinter <- tmesh3d(
#'   vertices = t(inter[["vertices"]]),
#'   indices = t(inter[["faces"]]),
#'   homogeneous = FALSE
#' )
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(rglmesh1, color = "yellow", alpha = 0.2)
#' shade3d(cube, color = "cyan", alpha = 0.2)
#' shade3d(rglinter, color = "red")
#' plotEdges(
#'   vertices = inter[["vertices"]], edges = inter[["exteriorEdges"]],
#'   edgesAsTubes = FALSE, lwd = 3, verticesAsSpheres = FALSE
#' )
MeshesIntersection <- function(
    meshes, merge = FALSE, normals = FALSE, exact = FALSE
){
  stopifnot(is.list(meshes))
  checkMeshes <- lapply(meshes, function(mesh){
    checkMesh(mesh[["vertices"]], mesh[["faces"]])
  })
  areTriangle <- all(vapply(checkMeshes, `[[`, logical(1L), "isTriangle"))
  if(!areTriangle){
    stop("All meshes must be triangular.")
  }
  meshes <- lapply(checkMeshes, `[`, c("vertices", "faces"))
  if(exact){
    inter <- Intersection2_EK(meshes, merge, normals)
  }else{
    inter <- Intersection2_K(meshes, merge, normals)
  }
  inter[["vertices"]] <- t(inter[["vertices"]])
  edges <- unname(t(inter[["edges"]]))
  inter[["exteriorEdges"]] <- edges[edges[, 3L] == 1L, c(1L, 2L)]
  inter[["edges"]] <- edges[, c(1L, 2L)]
  inter[["faces"]] <- do.call(rbind, inter[["faces"]])
  if(normals){
    tmesh <- vcgUpdateNormals(tmesh3d(
      vertices = t(inter[["vertices"]]),
      indices = t(inter[["faces"]]),
      homogeneous = FALSE
    ))
    inter[["normals"]] <- t(tmesh[["normals"]][-4L, ])
  }
  inter
}

#' @title Meshes difference
#' @description Computes the difference between two meshes.
#'
#' @param mesh1,mesh2 two \emph{triangular} meshes, each given as a list with
#'   (at least) two fields: \code{vertices} and \code{faces}
#' @param merge Boolean, whether to merge the duplicated vertices of the
#'   input meshes and the output mesh
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#' @param exact Boolean, whether to use exact calculations; this is slower but
#'   more accurate
#'
#' @return A triangular mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges} and \code{normals}
#'   if \code{normals=TRUE}.
#'
#' @importFrom rgl tmesh3d
#' @importFrom Rvcg vcgUpdateNormals
#'
#' @export
#'
#' @examples
#' library(RCGAL)
#' library(rgl)
#'
#' # mesh one: a cube; one has to triangulate it
#' cube1 <- cube3d() # (from the rgl package)
#' vertices <- t(cube1$vb[-4L, ])
#' faces <- t(cube1$ib)
#' mesh1 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
#'
#' # mesh two: another cube; one also has to triangulate it
#' cube2 <- translate3d( # (from the rgl package)
#'   cube3d(), 1, 1, 0
#' )
#' vertices <- t(cube2$vb[-4L, ])
#' faces <- t(cube2$ib)
#' mesh2 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
#'
#' # compute the difference
#' differ <- MeshesDifference(mesh1, mesh2)
#'
#' # plot
#' rgldiffer <- tmesh3d(
#'   vertices = t(differ[["vertices"]]),
#'   indices = t(differ[["faces"]]),
#'   homogeneous = FALSE
#' )
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(cube1, color = "yellow", alpha = 0.2)
#' shade3d(cube2, color = "cyan", alpha = 0.2)
#' shade3d(rgdiffer, color = "red")
#' plotEdges(
#'   vertices = differ[["vertices"]], edges = differ[["exteriorEdges"]],
#'   edgesAsTubes = FALSE, lwd = 3, verticesAsSpheres = FALSE
#' )
MeshesDifference <- function(
    mesh1, mesh2, merge = FALSE, normals = FALSE, exact = FALSE
){
  stopifnot(is.list(mesh1))
  stopifnot(is.list(mesh2))
  checkMesh1 <- checkMesh(mesh1[["vertices"]], mesh1[["faces"]])
  if(!checkMesh1[["isTriangle"]]){
    stop("The first mesh is not triangular.")
  }
  checkMesh2 <- checkMesh(mesh2[["vertices"]], mesh2[["faces"]])
  if(!checkMesh2[["isTriangle"]]){
    stop("The second mesh is not triangular.")
  }
  mesh1 <- checkMesh1[c("vertices", "faces")]
  mesh2 <- checkMesh2[c("vertices", "faces")]
  if(exact){
    differ <- Difference_EK(mesh1, mesh2, merge, normals)
  }else{
    differ <- Difference_K(mesh1, mesh2, merge, normals)
  }
  differ[["vertices"]] <- t(differ[["vertices"]])
  edges <- unname(t(differ[["edges"]]))
  differ[["exteriorEdges"]] <- edges[edges[, 3L] == 1L, c(1L, 2L)]
  differ[["edges"]] <- edges[, c(1L, 2L)]
  differ[["faces"]] <- do.call(rbind, differ[["faces"]])
  if(normals){
    tmesh <- vcgUpdateNormals(tmesh3d(
      vertices = t(differ[["vertices"]]),
      indices = t(differ[["faces"]]),
      homogeneous = FALSE
    ))
    differ[["normals"]] <- t(tmesh[["normals"]][-4L, ])
  }
  differ
}

#' @title Plot some edges
#' @description Plot the given edges with \strong{rgl}.
#'
#' @param vertices a three-columns matrix giving the coordinates of the vertices
#' @param edges a two-columns integer matrix giving the edges by pairs of
#'   vertex indices
#' @param color a color for the edges
#' @param lwd line width, a positive number, ignored if \code{edgesAsTubes=TRUE}
#' @param edgesAsTubes Boolean, whether to draw the edges as tubes
#' @param tubesRadius the radius of the tubes when \code{edgesAsTubes=TRUE}
#' @param verticesAsSpheres Boolean, whether to draw the vertices as spheres
#' @param spheresRadius the radius of the spheres when
#'   \code{verticesAsSpheres=TRUE}
#' @param spheresColor the color of the spheres when
#'   \code{verticesAsSpheres=TRUE}
#'
#' @return No value.
#'
#' @importFrom rgl cylinder3d shade3d lines3d spheres3d
#' @export
#'
#' @examples
#' library(RCGAL)
#' library(rgl)
#' mesh <- Mesh(
#'   truncatedIcosahedron[["vertices"]],
#'   truncatedIcosahedron[["faces"]],
#'   triangulate = TRUE, normals = FALSE
#' )
#' # now we can plot the truncated icosahedron
#' tmesh <- tmesh3d(
#'   vertices = t(mesh[["vertices"]]),
#'   indices = t(mesh[["faces"]]),
#'   homogeneous = FALSE
#' )
#' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
#' shade3d(tmesh, color = "gold")
#' plotEdges(mesh[["vertices"]], mesh[["edges0"]], color = "navy")
plotEdges <- function(
  vertices,
  edges,
  color = "black",
  lwd = 2,
  edgesAsTubes = TRUE,
  tubesRadius = 0.03,
  verticesAsSpheres = TRUE,
  spheresRadius = 0.05,
  spheresColor = color
){
  for(i in 1L:nrow(edges)){
    edge <- edges[i, ]
    if(edgesAsTubes){
      tube <- cylinder3d(
        vertices[edge, ], radius = tubesRadius, sides = 90
      )
      shade3d(tube, color = color)
    }else{
      lines3d(vertices[edge, ], color = color, lwd = lwd)
    }
  }
  if(verticesAsSpheres){
    spheres3d(vertices, radius = spheresRadius, color = spheresColor)
  }
  invisible(NULL)
}

