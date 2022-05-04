a <- 0.94
mu <- 0.56
c <- 0.34
f <- function(x, y, z, a, c, mu) {
  b <- sqrt(a^2 - c^2)
  (x^2 + y^2 + z^2 - mu^2 + b^2)^2 - 4 * (a * x - c * mu)^2 - 4 * b^2 * y^2
}
x <- seq(-c - mu - a, abs(mu - c) + a, length.out = 45)
y <- seq(-mu - a, mu + a, length.out = 45)
z <- seq(-mu - c, mu + c, length.out = 30)
g <- expand.grid(x = x, y = y, z = z)
voxel <- array(with(g, f(x, y, z, a, c, mu)), c(45, 45, 30))

library(misc3d)
cont <- computeContour3d(voxel, level = 0, x = x, y = y, z = z)
idx <- matrix(1:nrow(cont), ncol = 3, byrow = TRUE)

library(RCGAL)
mesh <- Mesh(cont, idx, merge = TRUE)

library(rgl)
tmesh <- tmesh3d(
  vertices = t(mesh[["vertices"]]),
  indices = t(mesh[["faces"]]),
  normals = mesh[["normals"]],
  homogeneous = FALSE
)

shade3d(tmesh, color = "green")
