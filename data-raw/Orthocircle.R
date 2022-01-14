library(misc3d)
library(rgl)

f <- function(x, y, z, a, b){
  x2 <- x*x
  y2 <- y*y
  z2 <- z*z
  xy2 <- x2 + y2 - 1
  yz2 <- y2 + z2 - 1
  zx2 <- z2 + x2 - 1
  (xy2*xy2 + z2) * (yz2*yz2 + x2) * (zx2*zx2 + y2) - a*a*(1 + b*(x2 + y2 + z2))
}

a = 0.075; b = 3

nx <- 100; ny <- 100; nz <- 100
x <- seq(-1.3, 1.3, length.out = nx)
y <- seq(-1.3, 1.3, length.out = ny)
z <- seq(-1.3, 1.3, length.out = nz)
G <- expand.grid(x = x, y = y, z = z)
voxel <- array(with(G, f(x, y, z, a, b)), c(nx, ny, nz))

surf <- computeContour3d(voxel, level = 0, x = x, y = y, z = z)

mesh0 <- tmesh3d(
  vertices = t(surf),
  indices = matrix(1L:nrow(surf), nrow = 3L)
)

mesh <- Rvcg::vcgUniformRemesh(mesh0)
Orthocircle <- t(mesh$vb[-4L, ])
