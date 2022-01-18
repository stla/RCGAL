library(RCGAL)
library(rgl)
library(viridisLite)

xy <- as.matrix(expand.grid(x = 1L:nrow(volcano), y = 1L:ncol(volcano)))
z <-  c(volcano)
pts <- cbind(5 * xy, z)

del <- delaunay(pts, elevation = TRUE)
mesh <- del[["mesh"]]

squaredNorms <- apply(mesh[["vb"]][-4L, ], 2L, crossprod)
normalizedSquaredNorms <- squaredNorms /  max(squaredNorms)

palette <- function(x) {
  RGB <- colorRamp(turbo(256))(x)
  rgb(RGB, maxColorValue = 255)
}

mesh$material <- list(color = palette(normalizedSquaredNorms))

open3d(windowRect = c(50, 50, 562, 562))
view3d(0, -40, zoom = 0.9)
wire3d(mesh, color = "black")
shade3d(mesh)
