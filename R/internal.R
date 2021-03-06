distance <- function(A, B){
  sqrt(c(crossprod(A-B)))
}

triangleArea <- function(A, B, C){
  a <- distance(B, C)
  b <- distance(A, C)
  c <- distance(A, B)
  s <- (a + b + c) / 2
  sqrt(s*(s-a)*(s-b)*(s-c))
}

isFalsy <- function(x){
  isFALSE(x) || is.null(x) || is.na(x)
}

makeTriangle <- function(vertices, indices){
  vertices[indices, ]
}

subtractEdges <- function(Edges, edges){
  if(is.null(edges)){
    return(Edges)
  }
  if(nrow(Edges) == nrow(edges)){
    return(NULL)
  }
  Strings <- paste0(Edges[, 1L], "-", Edges[, 2L])
  strings <- paste0(edges[, 1L], "-", edges[, 2L])
  rownames(Edges) <- Strings
  keep <- setdiff(Strings, strings)
  Edges[keep, ]
}

unionEdges <- function(edges1, edges2){
  if(is.null(edges2)){
    return(edges1)
  }
  Edges <- rbind(edges1, edges2)
  Edges[!duplicated(Edges), ]
  # strings1 <- paste0(edges1[, 1L], "-", edges1[, 2L])
  # strings2 <- paste0(edges2[, 1L], "-", edges2[, 2L])
  # strings <- union(strings1, strings2)
}

isAtomicVector <- function(x){
  is.atomic(x) && is.vector(x)
}

isPositiveNumber <- function(x){
  is.numeric(x) && length(x) == 1L && x > 0 && !is.na(x)
}

isNonNegativeNumber <- function(x){
  is.numeric(x) && length(x) == 1L && x >= 0 && !is.na(x)
}

isPositiveInteger <- function(x){
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x
}

isStrictPositiveInteger <- function(x){
  isPositiveInteger(x) && x > 0
}

isBoolean <- function(x){
  is.logical(x) && length(x) == 1L && !is.na(x)
}

#' @importFrom gmp `%*%`
#' @noRd
`%^%` <- function(A, n){
  Reduce(gmp::`%*%`, replicate(n, A, simplify = FALSE))
}
