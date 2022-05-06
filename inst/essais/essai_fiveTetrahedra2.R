library(rgl)

plotEdges <- function(vertices, edges, ...){
  for(i in 1L:nrow(edges)){
    edge <- edges[i, ]
    lines3d(rbind(vertices[edge[1L], ], vertices[edge[2L], ]), ...)
  }
  invisible(NULL)
}

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

ii <- RCGAL:::Intersection2(
  list(mesh1, mesh2),#, mesh3, mesh4, mesh5),
  FALSE, FALSE
)
vv = ii$vertices
ff = lapply(ii$faces, function(x) x-5L)
mm = list(vertices = vv, faces = ff)
iii <- RCGAL:::Intersection2(
  list(
    mm,
    mesh3
  ),
  TRUE, FALSE
)

tmesh <- tmesh3d(
  vertices = ii$vertices,
  indices = do.call(cbind, lapply(ii$faces, function(x) x-4L)),
  homogeneous = FALSE
)
t2 <- tmesh3d(
  vertices = mesh2$vertices,
  indices = t(faces),
  homogeneous = FALSE
)
t3 <- tmesh3d(
  vertices = mesh3$vertices,
  indices = t(faces),
  homogeneous = FALSE
)
t4 <- tmesh3d(
  vertices = mesh4$vertices,
  indices = t(faces),
  homogeneous = FALSE
)


open3d()
shade3d(t3, color = "cyan", alpha = 0.2)
shade3d(t2, color = "palegreen", alpha = 0.4)
shade3d(t4, color = "yellow", alpha = 0.2)

shade3d(tmesh, color="red")
#plotEdges(t(ii$vertices), t(ii$edges0))
wire3d(tmesh, lwd=3)
