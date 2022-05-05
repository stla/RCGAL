library(rgl)

a <- 1/sqrt(3)
vertices2 <- rbind(
  c(a, -a, -a),
  c(a, a, a),
  c(-a, -a, a),
  c(-a, a, -a)
)
faces <- rbind(
  c(1, 2, 3),
  c(3, 2, 4),
  c(4, 2, 1),
  c(1, 3, 4)
)

t1 <- scale3d(Rvcg::vcgSphere(subdivision = 4), 0.5, 0.5, 0.5)
t2 <- tmesh3d(
  vertices = t(vertices2),
  indices = t(faces),
  homogeneous = FALSE
)
open3d()
shade3d(t1, color = "cyan", alpha = 0.2)
shade3d(t2, color = "palegreen", alpha = 0.4)

lfaces <- lapply(1:nrow(faces), function(i) as.integer(faces[i, ]-1))
mesh1 <- list(
  vertices = t1$vb[1:3, ],
  faces = lapply(1:ncol(t1$it), function(i) t1$it[, i]-1L)
)
mesh2 <- list(
  vertices = t(vertices2),
  faces = lfaces
)
ii <- RCGAL:::Intersection(list(mesh1, mesh2), FALSE, TRUE)

tmesh <- tmesh3d(
  vertices = ii$vertices,
  indices = do.call(cbind, ii$faces),
  normals = t(ii$vertices),
  homogeneous = FALSE
)
Morpho::plotNormals(tmesh, long = 0.05)
shade3d(tmesh, color="red")
wire3d(tmesh, lwd=3)
