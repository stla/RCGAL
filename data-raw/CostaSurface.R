library(elliptic)
e1 <- Re(P(1/2, Omega = c(1/2, 1i/2)))
c <- 4*e1^2
fx <- function(u, v){
  w <- u + 1i*v
  Re(pi*(u + pi/4/e1) - zeta(w, c(c, 0)) +
       pi/2/e1*(zeta(w-1/2, c(c, 0)) - zeta(w-1i/2, c(c, 0))))/2
}
fy <- function(u, v) fx(v, u)
fz <- function(u, v){
  w <- u + 1i*v
  p <- P(w, c(c, 0))
  sqrt(pi/2)*log(Mod((p-e1)/(p+e1)))/2
}

library(misc3d)
library(rgl)
tris <- parametric3d(Vectorize(fx), Vectorize(fy), Vectorize(fz),
             umin = 0.05, umax = 0.95, vmin = 0.05, vmax = 0.95, n = 50,
             engine = "none")
mesh00 <- misc3d:::t2ve(tris)
mesh0 <- tmesh3d(mesh00$vb, mesh00$ib)
mesh <- Rvcg::vcgUniformRemesh(mesh0)
CostaSurface <- t(mesh$vb[-4L, ])
