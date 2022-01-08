library(cxhull)
library(geometry)
library(tessellation)
library(RCGAL)
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

# Elevated Delaunay:
f <- function(x, y){
  2 * exp(-(x^2 + y^2)) # integrate to 2pi
}
x <- y <- seq(-4, 4, length.out = 30)
grd <- transform(expand.grid(x = x, y = y), z = f(x, y))
microbenchmark(
  RCGAL = delaunay(as.matrix(grd), elevation = TRUE),
  deldir = deldir::triang.list(deldir::deldir(grd[["x"]], grd[["y"]], z = grd[["z"]])),
  times = 2
)
