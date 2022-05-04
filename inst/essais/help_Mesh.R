library(RCGAL)
library(rgl)

# a tetrahedron with ill-oriented faces ####
vertices <- rbind(
  c(-1, -1, -1),
  c(1, 1, -1),
  c(1, -1, 1),
  c(-1, 1, 1)
)
faces <- rbind(
  c(1, 2, 3),
  c(3, 4, 2),
  c(4, 2, 1),
  c(4, 3, 1)
)

# plot the tetrahedron, hiding the back of the faces
# then some faces do not appear, as their orientation is not correct
tmesh1 <- tmesh3d(
  vertices = t(vertices),
  indices = t(faces),
  homogeneous = FALSE
)
open3d(windowRect = c(50, 50, 562, 562))
shade3d(tmesh1, color = "green", back = "cull")

# now run the `Mesh` function
mesh2 <- Mesh(vertices, faces, normals = FALSE)
# plot the tetrahedron, hiding the back of the faces
# then all faces appear now
tmesh2 <- tmesh3d(
  vertices = t(mesh2[["vertices"]]),
  indices = t(mesh2[["faces"]]),
  homogeneous = FALSE
)
open3d(windowRect = c(50, 50, 562, 562))
shade3d(tmesh2, color = "blue", back = "cull")

# illustration of the `merge` option ####
# we construct a mesh with a lot of duplicated vertices
library(misc3d) # to compute a mesh of an isosurface
a <- 0.94; mu <- 0.56; c <- 0.34 # cyclide parameters
f <- function(x, y, z, a, c, mu){ # implicit equation of the cyclide
  b <- sqrt(a^2 - c^2)
  (x^2 + y^2 + z^2 - mu^2 + b^2)^2 - 4*(a*x - c*mu)^2 - 4*b^2*y^2
}
x <- seq(-c - mu - a, abs(mu - c) + a, length.out = 45)
y <- seq(-mu - a, mu + a, length.out = 45)
z <- seq(-mu - c, mu + c, length.out = 30)
g <- expand.grid(x = x, y = y, z = z)
voxel <- array(with(g, f(x, y, z, a, c, mu)), c(45, 45, 30))
cont <- computeContour3d(voxel, level = 0, x = x, y = y, z = z)
ids <- matrix(1:nrow(cont), ncol = 3, byrow = TRUE)
# run the `Mesh` function with `merge=TRUE`
mesh <- Mesh(cont, ids, merge = TRUE)
# plot the cyclide
tmesh <- tmesh3d(
  vertices = t(mesh[["vertices"]]),
  indices = t(mesh[["faces"]]),
  normals = mesh[["normals"]],
  homogeneous = FALSE
)
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
shade3d(tmesh, color = "green")

# illustration of the `triangulate` option ####
# the faces of the truncated icosahedron are hexagonal or pentagonal:
truncatedIcosahedron[["faces"]]
# so we triangulate them:
mesh <- Mesh(
  truncatedIcosahedron[["vertices"]],
  truncatedIcosahedron[["faces"]],
  triangulate = TRUE, normals = FALSE
)
# now we can plot the truncated icosahedron
tmesh <- tmesh3d(
  vertices = t(mesh[["vertices"]]),
  indices = t(mesh[["faces"]]),
  homogeneous = FALSE
)
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
shade3d(tmesh, color = "orange")
