library(gmp)
library(RCGAL)
library(rgl)

Vertices <- function(a, b, c, phi, O){
  t(cbind(
    c( a,  a,  a),
    c( a,  a, -a),
    c( a, -a,  a),
    c(-a, -a,  a),
    c(-a,  a, -a),
    c(-a,  a,  a),
    c( O,  b, -c),
    c( O, -b, -c),
    c( O, -b,  c),
    c( c,  O, -b),
    c(-c,  O, -b),
    c(-c,  O,  b),
    c( b,  c,  O),
    c( b, -c,  O),
    c(-b, -c,  O),
    c(-b,  c,  O),
    c( O,  b,  c),
    c( a, -a, -a),
    c( c,  O,  b),
    c(-a, -a, -a)
  ))
}

phi <- (1+sqrt(5))/2
a <- 1 / sqrt(3)
b <- a / phi
c <- a * phi
O <- 0
vertices <- Vertices(a, b, c, phi, O)

phi <- qsqrtPhi(17)
a <- 1L/qsqrt3(22)
b <- a / phi
c <- a * phi
O <- as.bigq(0)
gmpVertices <- Vertices(a, b, c, phi, O)
# check:
asNumeric(gmpVertices) - vertices

faces <- rbind(
  c(1L, 2L, 3L),
  c(3L, 2L, 4L),
  c(4L, 2L, 1L),
  c(1L, 3L, 4L)
)

mesh1 <- list(
  "vertices" = vertices[c(17, 14, 2, 11), ],
  "faces" = faces
)
mesh2 <- list(
  "vertices" = vertices[c(18, 1, 4, 5), ],
  "faces" = faces
)
mesh3 <- list(
  "vertices" = vertices[c(19, 6, 15, 7), ],
  "faces" = faces
)
mesh4 <- list(
  "vertices" = vertices[c(3, 13, 12, 8), ],
  "faces" = faces
)
mesh5 <- list(
  "vertices" = vertices[c(20, 16, 10, 9), ],
  "faces" = faces
)

meshes <- list(
  mesh1, mesh2, mesh3, mesh4, mesh5
)

tetrahedraCompound <- list(
  "meshes" = lapply(meshes, function(mesh){
    list("vertices" = as.bigq(mesh[["vertices"]]), faces = faces)
  }),
  "rglmeshes" = list(
    tmesh3d(
      "vertices" = t(mesh1[["vertices"]]),
      "indices" = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices" = t(mesh2[["vertices"]]),
      "indices" = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices" = t(mesh3[["vertices"]]),
      "indices" = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices" = t(mesh4[["vertices"]]),
      "indices" = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices" = t(mesh5[["vertices"]]),
      "indices" = t(faces),
      "homogeneous" = FALSE
    )
  )
)



inter <- MeshesIntersection(tetrahedraCompound$meshes, numbersType = "gmp")

open3d(windowRect = c(50, 50, 562, 562))
bg3d("#363940")
invisible(lapply(
  tetrahedraCompound[["rglmeshes"]], shade3d,
  color = "yellow", alpha = 0.2
))
rglinter <- tmesh3d(
  "vertices" = t(inter[["vertices"]]),
  "indices" = t(inter[["faces"]]),
  "homogeneous" = FALSE
)
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
shade3d(rglinter, color = "gold")
plotEdges(
  inter[["vertices"]], inter[["exteriorEdges"]],
  only = inter[["exteriorVertices"]], color = "navy"
)
