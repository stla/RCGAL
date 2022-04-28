library(RCGAL)
library(misc3d)
library(rgl)

R <- function(u, v){
  cos(v)^2 * pmax(abs(sin(4*u)), 0.9 - 0.2*abs(cos(8*u)))
}
fx <- function(u, v){
  R(u, v) * cos(u) * cos(v)
}
fy <- function(u, v){
  R(u, v) * sin(u) * cos(v)
}
fz <- function(u, v){
  R(u, v) * sin(v) / 2
}

tris <- parametric3d(
  fx, fy, fz, umin = 0, umax = 2*pi, vmin = 0, vmax = 2*pi,
  engine = "none"
)

# mesh00 <- misc3d:::t2ve(tris)
# mesh0 <- tmesh3d(mesh00$vb, mesh00$ib)
# mesh0 <- Rvcg::vcgClean(mesh0, sel = 7)
mesh <- Rvcg::vcgUniformRemesh(mesh0, discretize = TRUE)
flower <- t(mesh0$vb[-4L, ])

psr <- PoissonReconstruction(flower, spacing = 0.005)
shade3d(psr, color = "hotpink")
wire3d(psr)

