#' @title Rational approximation of square roots
#' @description Returns a rational approximation of the square root of
#'   an integer.
#'
#' @param x the positive integer whose square root is desired
#' @param n a positive integer, the higher the better approximation
#' @param ... ignored
#'
#' @return The \code{qsqrt} function returns a \strong{gmp} rational number
#'   (class \code{\link[gmp]{bigq}}) approximating the square root of
#'   \code{x}. The \code{qsqrt2}, \code{qsqrt3}, and \code{qsqrtPhi} functions
#'   return a \strong{gmp} rational number approximating the square root of
#'   \code{2}, \code{3}, and \code{phi} (the golden number) respectively.
#'   Their value converge more fastly than the value obtained with \code{qsqrt}.
#'
#' @importFrom gmp as.bigz as.bigq asNumeric matrix.bigq add.bigq factorialZ `%*%`
#'
#' @aliases qsqrt qsqrt2 qsqrt3 qsqrtPhi print.qsqrt
#' @rdname qsqrt
#' @export
#'
#' @examples
#' library(RCGAL)
#' qsqrt(2, 7)
#' qsqrt2(7)
#' qsqrt3(22)
#' qsqrtPhi(17)
qsqrt <- function(x, n){
  stopifnot(isPositiveInteger(x))
  stopifnot(isStrictPositiveInteger(n))
  zero <- as.bigz(0L)
  one <- as.bigz(1L)
  A <- matrix.bigq(
    c(zero, as.bigz(x)-1L, one, as.bigz(2L)),
    nrow = 2L, ncol = 2L
  )
  zs <- c(gmp::`%*%`(A %^% n, c(zero, one)))
  out <- as.bigq(zs[2L], zs[1L]) - 1L
  attr(out, "error") <- abs(asNumeric(out) - sqrt(x))
  class(out) <- c("qsqrt", class(out))
  out
}

#' @rdname qsqrt
#' @export
qsqrt2 <- function(n){
  stopifnot(isStrictPositiveInteger(n) && n >= 2)
  out <- as.bigq(99L, 70L) + Reduce(add.bigq, sapply(2L:n, function(i){
    numer <- 10L * prod(1L - 2L*(0L:(i-1L)))
    denom <- 7L * factorialZ(i) * as.bigz(-100L)^i
    as.bigq(numer, denom)
  }))
  attr(out, "error") <- abs(asNumeric(out) - sqrt(2))
  class(out) <- c("qsqrt", class(out))
  out
}

#' @rdname qsqrt
#' @export
qsqrt3 <- function(n){
  stopifnot(isStrictPositiveInteger(n) && n >= 2)
  out <- as.bigq(7L, 4L) + Reduce(add.bigq, sapply(2L:n, function(i){
    as.bigq(
      2L * prod(1L - 2L*(0L:(i-1L))),
      factorialZ(i) * as.bigz(-8L)^i
    )
  }))
  attr(out, "error") <- abs(asNumeric(out) - sqrt(3))
  class(out) <- c("qsqrt", class(out))
  out
}

#' @rdname qsqrt
#' @export
qsqrtPhi <- function(n){
  stopifnot(isStrictPositiveInteger(n) && n >= 2)
  out <- as.bigq(13L, 8L) + Reduce(add.bigq, sapply(2L:n, function(i){
    as.bigq(
      10L * prod(1L - 2L*(0L:(i-1L))),
      as.bigz(8L) * factorialZ(i) * as.bigz(-10L)^(i)
    )
  }))
  attr(out, "error") <- abs(asNumeric(out) - (1+sqrt(5))/2)
  class(out) <- c("qsqrt", class(out))
  out
}

#' @rdname qsqrt
#' @exportS3Method print qsqrt
print.qsqrt <- function(x, ...){
  print(as.bigq(x))
  cat('attr("error")')
  print(attr(x, "error"))
  invisible(NULL)
}

#' @importFrom gmp is.bigq is.matrixZQ
#' @importFrom data.table uniqueN
#' @noRd
checkMesh <- function(vertices, faces, gmp){
  if(gmp){
    if(!is.matrixZQ(vertices) || ncol(vertices) != 3L){
      stop("The `vertices` argument must be a matrix with three columns.")
    }
    stopifnot(is.bigq(vertices))
    vertices <- as.character(vertices)
  }else{
    if(!is.matrix(vertices) || ncol(vertices) != 3L){
      stop("The `vertices` argument must be a matrix with three columns.")
    }
    stopifnot(is.numeric(vertices))
    storage.mode(vertices) <- "double"
  }
  if(anyNA(vertices)){
    stop("Found missing values in `vertices`.")
  }
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
#' @param clean Boolean, whether to clean the mesh (merging duplicated
#'   vertices, duplicated faces, removed isolated vertices)
#' @param normals Boolean, whether to compute the normals
#' @param numbersType the type of the numbers used in C++ for the
#'   computations; must be one of \code{"double"}, \code{"lazyExact"}
#'   (a type provided by CGAL for exact computations), or \code{"gmp"}
#'   (exact computations with rational numbers); using exact computations can
#'   improve the detection of the exterior edges
#' @param epsilon if the mesh is triangulated or if \code{triangulate=TRUE},
#'   then \code{epsilon} is used in the detection of exterior edges (see the
#'   \strong{Value} section); the higher value of \code{epsilon}, the lower
#'   number of exterior edges
#'
#' @return A list giving the vertices, the edges, the faces of the mesh, the
#'   exterior edges, the exterior vertices and optionally the normals. This
#'   list has two additional components \code{edges0} and \code{normals0} if
#'   \code{triangulate=TRUE}, giving the edges and the normals before the
#'   triangulation, unless the mesh is already triangulated, in which case
#'   the \code{triangulate} option is ignored.
#'
#' @export
#'
#' @importFrom gmp as.bigq asNumeric
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
#' # illustration of the `clean` option ####
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
#' # run the `Mesh` function with `clean=TRUE`
#' mesh <- Mesh(cont, ids, clean = TRUE)
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
    vertices, faces, triangulate = FALSE, clean = FALSE, normals = FALSE,
    numbersType = "double", epsilon = 0
){
  numbersType <- match.arg(numbersType, c("double", "lazyExact", "gmp"))
  gmp <- numbersType == "gmp"
  stopifnot(epsilon >= 0)
  checkedMesh <- checkMesh(vertices, faces, gmp = gmp)
  vertices <- checkedMesh[["vertices"]]
  faces <- checkedMesh[["faces"]]
  homogeneousFaces <- checkedMesh[["homogeneousFaces"]]
  isTriangle <- checkedMesh[["isTriangle"]]
  rmesh <- list("vertices" = vertices, "faces" = faces)
  if(numbersType == "double"){
    mesh <- SurfMesh(
      rmesh, isTriangle, triangulate, clean, normals, epsilon
    )
  }else if(numbersType == "lazyExact"){
    mesh <- SurfEMesh(
      rmesh, isTriangle, triangulate, clean, normals, epsilon
    )
  }else{
    mesh <- SurfQMesh(
      rmesh, isTriangle, triangulate, clean, normals, epsilon
    )
  }
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
  if(gmp){
    vertices <- as.bigq(t(mesh[["vertices"]]))
    mesh[["gmpVertices"]] <- vertices
    vertices <- asNumeric(vertices)
  }else{
    vertices <- t(mesh[["vertices"]])
  }
  mesh[["vertices"]] <- vertices
  edges <- unname(t(mesh[["edges"]]))
  exteriorEdges <- edges[edges[, 3L] == 1L, c(1L, 2L)]
  mesh[["exteriorEdges"]] <- exteriorEdges
  mesh[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  mesh[["edges"]] <- edges[, c(1L, 2L)]
  if(normals){
    mesh[["normals"]] <- t(tmesh[["normals"]])
  }
  if(triangulate){
    mesh[["edges0"]] <- t(mesh[["edges0"]])
    if(normals){
      mesh[["normals0"]] <- t(mesh[["normals0"]])
    }
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
#'   (at least) two fields: \code{vertices} and \code{faces}; the `vertices`
#'   matrix must have the \code{bigq} class if \code{numberTypes="gmp"},
#'   otherwise it must be numeric
#' @param clean Boolean, whether to clean the input meshes (merging
#'   duplicated vertices, duplicated faces, removed isolated vertices)
#'   as well as the output mesh
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#' @param numbersType the type of the numbers used in C++ for the
#'   computations; must be one of \code{"double"}, \code{"lazyExact"}
#'   (a type provided by CGAL for exact computations), or \code{"gmp"}
#'   (exact computations with rational numbers); of course using
#'   exact computations is slower but more accurate
#'
#' @return A triangular mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges}, \code{gmpvertices}
#'   if \code{numberTypes="gmp"}, and \code{normals} if \code{normals=TRUE}.
#'
#' @importFrom gmp as.bigq asNumeric
#' @importFrom rgl tmesh3d
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
#'
#' # other example, with 'gmp' rational numbers ####
#' library(RCGAL)
#' library(gmp)
#' library(rgl)
#'
#' cube <- cube3d()
#' vertices <- t(cube$vb[-4L, ])
#' faces <- t(cube$ib)
#'
#' rglmesh1 <- cube
#' mesh1 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
#' mesh1$vertices <- as.bigq(mesh1$vertices)
#'
#' rotMatrix <- t(cbind( # pi/3 around a great diagonal
#'   as.bigq(c(2, -1, 2), c(3, 3, 3)),
#'   as.bigq(c(2, 2, -1), c(3, 3, 3)),
#'   as.bigq(c(-1, 2, 2), c(3, 3, 3))
#' ))
#' mesh2 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
#' mesh2$vertices <- as.bigq(vertices) %*% rotMatrix
#' rglmesh2 <- rotate3d(cube, pi/3, 1, 1, 1)
#'
#' inter <- MeshesIntersection(list(mesh1, mesh2), numbersType = "gmp")
#' # perfect vertices:
#' inter[["vertices"]]
#' rglinter <- tmesh3d(
#'   vertices = t(inter[["vertices"]]),
#'   indices = t(inter[["faces"]]),
#'   homogeneous = FALSE
#' )
#'
#' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
#' bg3d("#363940")
#' shade3d(rglmesh1, color = "yellow", alpha = 0.2)
#' shade3d(rglmesh2, color = "orange", alpha = 0.2)
#' shade3d(rglinter, color = "hotpink")
#' plotEdges(
#'   inter[["vertices"]], inter[["exteriorEdges"]],
#'   only = inter[["exteriorVertices"]],
#'   color = "firebrick",
#'   tubesRadius = 0.05, spheresRadius = 0.07
#' )
MeshesIntersection <- function(
    meshes, clean = FALSE, normals = FALSE, numbersType = "double"
){
  numbersType <- match.arg(numbersType, c("double", "lazyExact", "gmp"))
  gmp <- numbersType == "gmp"
  stopifnot(is.list(meshes))
  checkMeshes <- lapply(meshes, function(mesh){
    checkMesh(mesh[["vertices"]], mesh[["faces"]], gmp)
  })
  areTriangle <- all(vapply(checkMeshes, `[[`, logical(1L), "isTriangle"))
  if(!areTriangle){
    stop("All meshes must be triangular.")
  }
  meshes <- lapply(checkMeshes, `[`, c("vertices", "faces"))
  if(numbersType == "double"){
    inter <- Intersection2_K(meshes, clean, normals)
  }else if(numbersType == "lazyExact"){
    inter <- Intersection2_EK(meshes, clean, normals)
  }else{
    inter <- Intersection_Q(meshes, clean, normals)
  }
  if(gmp){
    vertices <- as.bigq(t(inter[["vertices"]]))
    inter[["gmpVertices"]] <- vertices
    vertices <- asNumeric(vertices)
  }else{
    vertices <- t(inter[["vertices"]])
  }
  inter[["vertices"]] <- vertices
  edges <- unname(t(inter[["edges"]]))
  exteriorEdges <- edges[edges[, 3L] == 1L, c(1L, 2L)]
  inter[["exteriorEdges"]] <- exteriorEdges
  inter[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  inter[["edges"]] <- edges[, c(1L, 2L)]
  inter[["faces"]] <- do.call(rbind, inter[["faces"]])
  if(normals){
    inter[["normals"]] <- t(tmesh[["normals"]])
  }
  # if(normals){
  #   tmesh <- vcgUpdateNormals(tmesh3d(
  #     vertices = vertices,
  #     indices = t(inter[["faces"]]),
  #     homogeneous = FALSE
  #   ))
  #   inter[["normals"]] <- t(tmesh[["normals"]][-4L, ])
  # }
  inter
}

#' @title Meshes difference
#' @description Computes the difference between two meshes.
#'
#' @param mesh1,mesh2 two \emph{triangular} meshes, each given as a list with
#'   (at least) two fields: \code{vertices} and \code{faces}; the `vertices`
#'   matrix must have the \code{bigq} class if \code{numberTypes="gmp"},
#'   otherwise it must be numeric
#' @param clean Boolean, whether to clean the input mesh (merging duplicated
#'   vertices, duplicated faces, removed isolated vertices) as well as the
#'   output mesh
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#' @param numbersType the type of the numbers used in C++ for the
#'   computations; must be one of \code{"double"}, \code{"lazyExact"}
#'   (a type provided by CGAL for exact computations), or \code{"gmp"}
#'   (exact computations with rational numbers); of course using
#'   exact computations is slower but more accurate
#'
#' @return A triangular mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges}, \code{gmpvertices}
#'   if \code{numberTypes="gmp"}, and \code{normals} if \code{normals=TRUE}.
#'
#' @importFrom gmp as.bigq asNumeric
#' @importFrom rgl tmesh3d
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
#' shade3d(rgldiffer, color = "red")
#' plotEdges(
#'   vertices = differ[["vertices"]], edges = differ[["exteriorEdges"]],
#'   edgesAsTubes = TRUE, verticesAsSpheres = TRUE
#' )
MeshesDifference <- function(
    mesh1, mesh2, clean = FALSE, normals = FALSE, numbersType = "double"
){
  stopifnot(is.list(mesh1))
  stopifnot(is.list(mesh2))
  numbersType <- match.arg(numbersType, c("double", "lazyExact", "gmp"))
  gmp <- numbersType == "gmp"
  checkMesh1 <- checkMesh(mesh1[["vertices"]], mesh1[["faces"]], gmp)
  if(!checkMesh1[["isTriangle"]]){
    stop("The first mesh is not triangular.")
  }
  checkMesh2 <- checkMesh(mesh2[["vertices"]], mesh2[["faces"]], gmp)
  if(!checkMesh2[["isTriangle"]]){
    stop("The second mesh is not triangular.")
  }
  mesh1 <- checkMesh1[c("vertices", "faces")]
  mesh2 <- checkMesh2[c("vertices", "faces")]
  if(numbersType == "double"){
    differ <- Difference_K(mesh1, mesh2, clean, normals)
  }else if(numbersType == "lazyExact"){
    differ <- Difference_EK(mesh1, mesh2, clean, normals)
  }else{
    differ <- Difference_K(mesh1, mesh2, clean, normals)
  }
  if(gmp){
    vertices <- as.bigq(t(differ[["vertices"]]))
    differ[["gmpVertices"]] <- vertices
    vertices <- asNumeric(vertices)
  }else{
    vertices <- t(differ[["vertices"]])
  }
  differ[["vertices"]] <- vertices
  edges <- unname(t(differ[["edges"]]))
  exteriorEdges <- edges[edges[, 3L] == 1L, c(1L, 2L)]
  differ[["exteriorEdges"]] <- exteriorEdges
  differ[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  differ[["edges"]] <- edges[, c(1L, 2L)]
  differ[["faces"]] <- do.call(rbind, differ[["faces"]])
  if(normals){
    differ[["normals"]] <- t(tmesh[["normals"]])
  }
  differ
}

#' @title Meshes union
#' @description Computes the union of the given meshes.
#'
#' @param meshes a list of \emph{triangular} meshes, each given as a list with
#'   (at least) two fields: \code{vertices} and \code{faces}; the `vertices`
#'   matrix must have the \code{bigq} class if \code{numberTypes="gmp"},
#'   otherwise it must be numeric
#' @param clean Boolean, whether to clean the input meshes (merging duplicated
#'   vertices, duplicated faces, removed isolated vertices) as well as the
#'   output mesh
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#' @param numbersType the type of the numbers used in C++ for the
#'   computations; must be one of \code{"double"}, \code{"lazyExact"}
#'   (a type provided by CGAL for exact computations), or \code{"gmp"}
#'   (exact computations with rational numbers); of course using
#'   exact computations is slower but more accurate
#'
#' @return A triangular mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges}, \code{gmpvertices}
#'   if \code{numberTypes="gmp"}, and \code{normals} if \code{normals=TRUE}.
#'
#' @importFrom gmp as.bigq asNumeric
#' @importFrom rgl tmesh3d
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
#'   cube3d(), 1, 1, 1
#' )
#' vertices <- t(cube2$vb[-4L, ])
#' faces <- t(cube2$ib)
#' mesh2 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
#'
#' # compute the union
#' umesh <- MeshesUnion(list(mesh1, mesh2))
#'
#' # plot
#' rglumesh <- tmesh3d(
#'   vertices = t(umesh[["vertices"]]),
#'   indices = t(umesh[["faces"]]),
#'   homogeneous = FALSE
#' )
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(rglumesh, color = "red")
#' plotEdges(
#'   vertices = umesh[["vertices"]], edges = umesh[["exteriorEdges"]],
#'   edgesAsTubes = TRUE, verticesAsSpheres = TRUE
#' )
MeshesUnion <- function(
    meshes, clean = FALSE, normals = FALSE, numbersType = "double"
){
  stopifnot(is.list(meshes))
  numbersType <- match.arg(numbersType, c("double", "lazyExact", "gmp"))
  gmp <- numbersType == "gmp"
  checkMeshes <- lapply(meshes, function(mesh){
    checkMesh(mesh[["vertices"]], mesh[["faces"]], gmp)
  })
  areTriangle <- all(vapply(checkMeshes, `[[`, logical(1L), "isTriangle"))
  if(!areTriangle){
    stop("All meshes must be triangular.")
  }
  meshes <- lapply(checkMeshes, `[`, c("vertices", "faces"))
  if(numbersType == "double"){
    umesh <- Union_K(meshes, clean, normals)
  }else if(numbersType == "lazyExact"){
    umesh <- Union_EK(meshes, clean, normals)
  }else{
    umesh <- Intersection_Q(meshes, clean, normals)
  }
  if(gmp){
    vertices <- as.bigq(t(umesh[["vertices"]]))
    umesh[["gmpVertices"]] <- vertices
    vertices <- asNumeric(vertices)
  }else{
    vertices <- t(umesh[["vertices"]])
  }
  umesh[["vertices"]] <- vertices
  edges <- unname(t(umesh[["edges"]]))
  exteriorEdges <- edges[edges[, 3L] == 1L, c(1L, 2L)]
  umesh[["exteriorEdges"]] <- exteriorEdges
  umesh[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  umesh[["edges"]] <- edges[, c(1L, 2L)]
  umesh[["faces"]] <- do.call(rbind, umesh[["faces"]])
  if(normals){
    umesh[["normals"]] <- t(tmesh[["normals"]])
  }
  umesh
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
#' @param only integer vector made of the indices of the vertices you want
#'   to plot (as spheres), or \code{NULL} to plot all vertices
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
    only = NULL,
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
    if(!is.null(only)){
      vertices <- vertices[only, ]
    }
    spheres3d(vertices, radius = spheresRadius, color = spheresColor)
  }
  invisible(NULL)
}
