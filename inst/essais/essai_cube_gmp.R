# other example, with 'gmp' rational numbers ####
library(RCGAL)
library(gmp)
library(rgl)

cube <- cube3d()
vertices <- t(cube$vb[-4L, ])
faces <- t(cube$ib)

rglmesh1 <- cube
mesh1 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
mesh1$vertices <- as.bigq(mesh1$vertices)

rotMatrix <- t(cbind( # pi/3 around a great diagonal
  as.bigq(c(2, -1, 2), c(3, 3, 3)),
  as.bigq(c(2, 2, -1), c(3, 3, 3)),
  as.bigq(c(-1, 2, 2), c(3, 3, 3))
))
mesh2 <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
mesh2$vertices <- as.bigq(vertices) %*% rotMatrix
rglmesh2 <- rotate3d(cube, pi/3, 1, 1, 1)

inter <- MeshesIntersection(list(mesh1, mesh2), numbersType = "gmp")
# perfect vertices:
inter[["vertices"]]
rglinter <- tmesh3d(
  vertices = t(inter[["vertices"]]),
  indices = t(inter[["faces"]]),
  homogeneous = FALSE
)

open3d(windowRect = c(50, 50, 562, 562), zoom = 0.7)
bg3d("#363940")
shade3d(rglmesh1, color = "yellow", alpha = 0.2)
shade3d(rglmesh2, color = "orange", alpha = 0.2)
shade3d(rglinter, color = "hotpink")
plotEdges(
  inter[["vertices"]], inter[["exteriorEdges"]],
  only = inter[["exteriorVertices"]],
  color = "firebrick",
  tubesRadius = 0.05, spheresRadius = 0.07
)


stop()

# animation ####
movie3d(spin3d(axis = c(0, 1, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE,
        startTime = 1/10,
        webshot = FALSE)


command <- "gifski --fps=10 --frames=zzpic*.png -o interCubeRotatedCube.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)
















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
