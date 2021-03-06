library(RCGAL)
library(rgl)

# mesh one: truncated icosahedron; one has to triangulate it
mesh1 <- Mesh(
  truncatedIcosahedron[["vertices"]],
  truncatedIcosahedron[["faces"]],
  triangulate = TRUE, normals = FALSE
)

# mesh two: a cube; one also has to triangulate it
cube <- translate3d( # (from the rgl package)
  cube3d(), 2, 0, 0
)
vertices <- t(cube$vb[-4L, ])
faces <- t(cube$ib)
mesh2 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)

# compute the intersection
inter <- MeshesIntersection(list(mesh1, mesh2))

# plot
rglmesh1 <- tmesh3d(
  vertices = t(mesh1[["vertices"]]),
  indices = t(mesh1[["faces"]]),
  homogeneous = FALSE
)
rglinter <- tmesh3d(
  vertices = t(inter[["vertices"]]),
  indices = t(inter[["faces"]]),
  homogeneous = FALSE
)
open3d(windowRect = c(50, 50, 562, 562))
shade3d(rglmesh1, color = "yellow", alpha = 0.2)
shade3d(cube, color = "cyan", alpha = 0.2)
shade3d(rglinter, color = "red")
plotEdges(
  vertices = inter[["vertices"]], edges = inter[["exteriorEdges"]],
  edgesAsTubes = FALSE, lwd = 3, verticesAsSpheres = FALSE
)

