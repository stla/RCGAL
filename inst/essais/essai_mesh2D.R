library(RCGAL)
library(rgl)

# vertices
R <- sqrt((5-sqrt(5))/10)     # outer radius
r <- sqrt((25-11*sqrt(5))/10) # circumradius of the inner pentagon
k <- pi/180 # factor to convert degrees to radians
X <- R * vapply(0L:4L, function(i) cos(k * (90+72*i)), numeric(1L))
Y <- R * vapply(0L:4L, function(i) sin(k * (90+72*i)), numeric(1L))
x <- r * vapply(0L:4L, function(i) cos(k * (126+72*i)), numeric(1L))
y <- r * vapply(0L:4L, function(i) sin(k * (126+72*i)), numeric(1L))
vertices <- rbind(
  c(X[1L], Y[1L]),
  c(x[1L], y[1L]),
  c(X[2L], Y[2L]),
  c(x[2L], y[2L]),
  c(X[3L], Y[3L]),
  c(x[3L], y[3L]),
  c(X[4L], Y[4L]),
  c(x[4L], y[4L]),
  c(X[5L], Y[5L]),
  c(x[5L], y[5L])
)
# edge indices
edges <- cbind(1L:10L, c(2L:10L, 1L))
# constrained Delaunay triangulation
del <- delaunay(vertices, constraints = edges)

m <- mesh2d(del)
mesh <- m[["mesh"]]

open3d(windowRect = c(100, 100, 612, 612))
shade3d(mesh, color = "red", specular = "red")
wire3d(mesh, color = "black")

open3d(windowRect = c(100, 100, 612, 612))
shade3d(mesh, color = "red", specular = "red")
segments3d(m$borderEdges, lwd = 3)

# hexagon ####
nsides <- 6L
angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
outer_points <- cbind(cos(angles), sin(angles))
inner_points <- outer_points / 2
points <- rbind(outer_points, inner_points)
# constraint edges
indices <- 1L:nsides
edges <- cbind(
  indices,
  c(indices[-1L], indices[1L])
)
edges <- rbind(edges, edges + nsides)
# constrained Delaunay triangulation
del <- delaunay(points, constraints = edges)

m <- mesh2d(del)
mesh <- m[["mesh"]]

open3d(windowRect = c(100, 100, 612, 612))
shade3d(mesh, color = "red", specular = "red")
wire3d(mesh, color = "black", lwd = 3, specular = "black")

open3d(windowRect = c(100, 100, 612, 612))
shade3d(mesh, color = "red", specular = "red")
segments3d(m$borderEdges, lwd = 3)
