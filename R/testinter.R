testinter <- function() {
  phi <- (1+sqrt(5))/2
  a <- 1/sqrt(3)
  b <- a/phi
  c <- a*phi
  
  vertices <-
    1000 * rbind(
      c( a,  a,  a),
      c( a,  a, -a),
      c( a, -a,  a),
      c(-a, -a,  a),
      c(-a,  a, -a),
      c(-a,  a,  a),
      c( 0,  b, -c),
      c( 0, -b, -c),
      c( 0, -b,  c),
      c( c,  0, -b),
      c(-c,  0, -b),
      c(-c,  0,  b),
      c( b,  c,  0),
      c( b, -c,  0),
      c(-b, -c,  0),
      c(-b,  c,  0),
      c( 0,  b,  c),
      c( a, -a, -a),
      c( c,  0,  b),
      c(-a, -a, -a)
    )
vertices <- round(vertices, 11)

  vs1 <- vertices[c(17, 14, 2, 11), ]
  vs2 <- vertices[c(18, 1, 4, 5), ]
  vs3 <- vertices[c(19, 6, 15, 7), ]
  vs4 <- vertices[c(3, 13, 12, 8), ]
  vs5 <- vertices[c(20, 16, 10, 9), ]
  
  faces <- rbind(
    c(1, 2, 3),
    c(3, 2, 4),
    c(4, 2, 1),
    c(1, 3, 4)
  )
  lfaces <- lapply(1:nrow(faces), function(i) as.integer(faces[i, ]-1))
  lfacescopy = lapply(lfaces, identity)
  
  mesh1 <- list(
    vertices = t(vs1),
    faces = lfaces
  )
  mesh2 <- list(
    vertices = t(vs2),
    faces = lfaces
  )
  mesh3 <- list(
    vertices = t(vs3),
    faces = lfacescopy
  )
  mesh4 <- list(
    vertices = t(vs4),
    faces = lfaces
  )
  mesh5 <- list(
    vertices = t(vs5),
    faces = lfaces
  )
  
  ii <- Intersection2_K(
    list(mesh1, mesh2, mesh3),
    FALSE, FALSE
  )
  vv = ii$vertices
  ff = lapply(ii$faces, function(x) x-1L)
  mm = list(vertices = vv, faces = ff)
  iii <- Intersection2_K(
    list(
      mesh3,
      mesh4,
      mesh5
    ),
    FALSE, FALSE
  )
  vv = iii$vertices
  ff = lapply(iii$faces, function(x) x-1L)
  mmm = list(vertices = vv, faces = ff)
  
  iiii <- Intersection2_K(
    list(
      mm,
      mmm
    ),
    FALSE, FALSE
  )
  iiii
}
