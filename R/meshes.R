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
#' @importFrom data.table uniqueN
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
  mesh <- SurfMesh(
    t(vertices), faces, isTriangle, triangulate, merge, normals, epsilon
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
