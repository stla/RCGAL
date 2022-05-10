phi <- (1+sqrt(5))/2
a <- 1/sqrt(3)
b <- a/phi
c <- a*phi

vertices <-
  rbind(
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

faces <- rbind(
  c(1L, 2L, 3L),
  c(3L, 2L, 4L),
  c(4L, 2L, 1L),
  c(1L, 3L, 4L)
)

tetrahedraCompound <- list(
  list(
    vertices = vertices[c(17, 14, 2, 11), ],
    faces = faces
  ),
  list(
    vertices = vertices[c(18, 1, 4, 5), ],
    faces = faces
  ),
  list(
    vertices = vertices[c(19, 6, 15, 7), ],
    faces = faces
  ),
  list(
    vertices = vertices[c(3, 13, 12, 8), ],
    faces = faces
  ),
  list(
    vertices = vertices[c(20, 16, 10, 9), ],
    faces = faces
  )
)
