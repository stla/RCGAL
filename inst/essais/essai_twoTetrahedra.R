library(rgl)



phi <- (1+sqrt(5))/2
a <- 1/sqrt(3)
b <- a/phi
c <- a*phi;

vertices1 <- rbind(
  c(0, b, c),
  c(b, -c, 0),
  c(a, a, -a),
  c(-c, 0, -b)
)

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

t1 <- tmesh3d(
  vertices = t(vertices1),
  indices = t(faces),
  homogeneous = FALSE
)
t2 <- tmesh3d(
  vertices = t(vertices2),
  indices = t(faces),
  homogeneous = FALSE
)
open3d()
shade3d(t1, color = "cyan", alpha = 0.4)
shade3d(t2, color = "palegreen", alpha = 0.6)

lfaces <- lapply(1:nrow(faces), function(i) as.integer(faces[i, ]-1))
mesh1 <- list(
  vertices = vertices1,
  faces = lfaces
)
mesh2 <- list(
  vertices = vertices2,
  faces = lfaces
)
ii <- RCGAL:::Intersection(list(mesh1, mesh2), FALSE, TRUE)
