library(RCGAL)
library(Rvcg)
library(rgl)

cap <- vcgSphericalCap(angleRad = pi/2, subdivision = 3, normals = TRUE)
shade3d(cap, color = "green")
wire3d(cap)

R <- 1
h <- R*(1-sin(pi/2/2))
2*pi*R*h
vcgArea(cap)

points <- t(cap$vb[-4,])
del <- delaunay(points, elevation = TRUE)
del$area


# volume
pi*h^2/3 * (3*R-h)
del$volume
