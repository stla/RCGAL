a = 0.94; d = 0.56; c = 0.34; b = sqrt(a^2-c^2)
fx <- function(u, v){
  (d*(c-a*cos(u)*cos(v)) + b^2*cos(u)) / (a-c*cos(u)*cos(v))
}
fy <- function(u, v){
  (b*sin(u)*(a-d*cos(v))) / (a-c*cos(u)*cos(v))
}
fz <- function(u, v){
  (b*sin(v)*(c*cos(u)-d)) / (a-c*cos(u)*cos(v))
}
library(misc3d)
tris <- parametric3d(fx, fy, fz,
                     umin = 0, umax = 2*pi, vmin = 0, vmax = 2*pi,
                     n = 50, engine = "none")
mesh00 <- misc3d:::t2ve(tris)
mesh0 <- rgl::tmesh3d(mesh00$vb, mesh00$ib)
mesh <- Rvcg::vcgUniformRemesh(mesh0)
cyclide <- t(mesh$vb[-4L, ])
