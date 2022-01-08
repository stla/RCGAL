library(rgl)
library(Rvcg)

HT <- function(u, v, nlobes, A){
  B <- pi/2 - (pi/2-A)*cos(u*nlobes)
  C <- u + A*sin(2*u*nlobes)
  p1 <- cos(B)
  p2 <- sin(B) * cos(C)
  p3 <- sin(B) * sin(C)
  y1 <- 1 + p1
  y2 <- p2
  y3 <- p3
  cos_v <- cos(v)
  sin_v <- sin(v)
  x1 <- cos_v*y3 + sin_v*y2
  x2 <- cos_v*y2 - sin_v*y3
  x3 <- sin_v * y1
  x4 <- cos_v * y1
  yden <- sqrt(2*(1+p1))
  c(x1, x2, x3) / (yden-x4)
}

hmesh <- function(nu, nv, nlobes = 3, A = 0.44){
  vs <- matrix(NA_real_, nrow=3, ncol=nu*nv)
  u_ <- seq(0, 2*pi, length.out = nu+1)[-1L]
  v_ <- seq(0, 2*pi, length.out = nv+1)[-1L]
  for(i in 1:nu){
    for(j in 1:nv){
      vs[, (i-1L)*nv+j] <- HT(u_[i], v_[j], nlobes, A)
    }
  }
  tris1 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
  tris2 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
  for(i in 1L:nu){
    ip1 <- ifelse(i == nu, 1L, i+1L)
    for(j in 1L:nv){
      jp1 <- ifelse(j == nv, 1L, j+1L)
      tris1[, (i-1L)*nv+j] <- c((i-1L)*nv+j, (i-1L)*nv+jp1, (ip1-1L)*nv+j)
      tris2[, (i-1L)*nv+j] <- c((i-1L)*nv+jp1, (ip1-1L)*nv+jp1,(ip1-1L)*nv+j)
    }
  }
  addNormals(tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    normals = NULL
  ))
}

HopfTorus_Mesh <- vcgUniformRemesh(hmesh(200, 150))
HopfTorus <- t(HopfTorus_Mesh[["vb"]][-4L, ])
