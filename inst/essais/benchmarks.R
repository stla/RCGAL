library(RCGAL)
library(cxhull)
library(geometry)
library(tessellation)
library(microbenchmark)
#library(uniformly)

data(bunny, package="onion")

# Convexhull:
microbenchmark(
  RCGAL = convexhull(bunny),
  cxhull = cxhull(bunny, triangulate = TRUE),
  geometry = convhulln(bunny, output.options = TRUE),
  times = 2
)

# Delaunay:
microbenchmark(
  RCGAL = RCGAL::delaunay(bunny),
  tessellation = tessellation::delaunay(bunny), # plante
  geometry = delaunayn(bunny, output.options = TRUE),
  times = 2
)
