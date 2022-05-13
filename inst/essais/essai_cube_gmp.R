library(rgl)
library(RCGAL)
library(gmp)

cube <- cube3d()

vertices <- t(cube$vb[-4L, ])
faces <- t(cube$ib)
mesh1 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
mesh1$vertices <- as.bigq(mesh1$vertices)

# cube <- rotate3d( # (from the rgl package)
#   cube3d(), pi/3, 1, 1, 1
# )
vertices <- t(cube$vb[-4L, ])
faces <- t(cube$ib)

rotMatrix <- t(cbind(
  as.bigq(c(2, -1, 2), c(3, 3, 3)),
  as.bigq(c(2, 2, -1), c(3, 3, 3)),
  as.bigq(c(-1, 2, 2), c(3, 3, 3))
))
mesh2 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
mesh2$vertices <- as.bigq(vertices) %*% rotMatrix

inter <- MeshesIntersection(list(mesh1, mesh2), numbersType = "gmp")

inter2 <- Mesh(vertices = inter[["vertices"]], faces = inter[["faces"]], epsilon = 1e-13)
rglinter <- tmesh3d(
  vertices = t(inter[["vertices"]]),
  indices = t(inter[["faces"]]),
  homogeneous = FALSE
)
shade3d(rglinter, color = "red")
plotEdges(
  vertices = inter2[["vertices"]], edges = inter2[["exteriorEdges"]],
  edgesAsTubes = FALSE, lwd = 3, verticesAsSpheres = FALSE
)

hull <- convexhull(inter$vertices, epsilon = 1e-17)
plotConvexHull3D(hull)
