#' @title Make a 3D mesh
#' @description Make a 3D mesh from given vertices and faces; the returned
#'   faces are coherently oriented and normals are computed if desired.
#'
#' @param vertices a numeric matrix with three columns
#' @param faces either an integer matrix (each row provides the vertex indices
#'   of the corresponding face) or a list of integer vectors, each one
#'   providing the vertex indices of the corresponding face
#' @param normals Boolean, whether to compute the normals
#'
#' @return A list giving the vertices, the edges, the faces of the mesh, and
#'   optionally the normals.
#'
#' @importFrom data.table uniqueN
#' @export
#'
#' @examples
#' # a tetrahedron with ill-oriented faces:
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
#' Mesh(vertices, faces, normals = FALSE)
Mesh <- function(vertices, faces, normals = TRUE){
  if(!is.matrix(vertices) || ncol(vertices) != 3L){
    stop("The `vertices` argument must be a matrix with three columns.")
  }
  storage.mode(vertices) <- "double"
  homogeneousFaces <- FALSE
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
    faces <- lapply(1L:nrow(faces), function(i) faces[i, ])
  }else if(is.list(faces)){
    check <- all(vapply(faces, isAtomicVector, logical(1L)))
    if(!check){
      stop("The `faces` argument must be a list of integer vectors.")
    }
    check <- any(vapply(faces, anyNA, logical(1L)))
    if(check){
      stop("Found missing values in `faces`.")
    }
    faces <- lapply(faces, as.integer)
    sizes <- lengths(faces)
    if(any(sizes < 3L)){
      stop("Faces must be given by at least three indices.")
    }
    check <- any(vapply(faces, function(f){
      any(f < 1L) || any(f > nrow(vertices))
    }, logical(1L)))
    if(check){
      stop(
        "Faces cannot contain indices lower than 1 or higher than the ",
        "number of vertices."
      )
    }
    homogeneousFaces <- uniqueN(sizes) == 1L
  }else{
    stop("The `faces` argument must be a list or a matrix.")
  }
  mesh <- if(normals){
    SurfMeshWithNormals(t(vertices), faces)
  }else{
    SurfMesh(t(vertices), faces)
  }
  mesh[["vertices"]] <- t(mesh[["vertices"]])
  mesh[["edges"]] <- t(mesh[["edges"]])
  if(homogeneousFaces){
    mesh[["faces"]] <- do.call(rbind, mesh[["faces"]])
  }
  if("normals" %in% names(mesh)){
    mesh[["normals"]] <- t(mesh[["normals"]])
  }
  mesh
}
