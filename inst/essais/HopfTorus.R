library(RCGAL)
library(rgl)

HT <- function(u, v, nlobes, A){
  p2 <- sin(pi/2 - (pi/2-A)*cos(u*nlobes)) * cos(u+A*sin(2*u*nlobes))
  p3 <- sin(pi/2 - (pi/2-A)*cos(u*nlobes)) * sin(u+A*sin(2*u*nlobes))
  p1 <- cos(pi/2 - (pi/2-A)*cos(u*nlobes))
  yden <- sqrt(2*(1+p1))
  y1 <- (1+p1)/yden
  y2 <- p2/yden
  y3 <- p3/yden
  x4 <- cos(v) * y1
  x3 <- sin(v) * y1
  x2 <- cos(v)*y2 - sin(v)*y3
  x1 <- cos(v)*y3 + sin(v)*y2
  c(x1, x2, x3) / (1-x4)
}

xprod <- function(v, w){
  c(v[2]*w[3] - v[3]*w[2], v[3]*w[1] - v[1]*w[3], v[1]*w[2] - v[2]*w[1])
}

hmesh <- function(nu, nv, normals = TRUE, nlobes = 3, A = 0.44, ...){
  vs <- matrix(NA_real_, nrow=3, ncol=nu*nv)
  u_ <- seq(0, 2*pi, length.out = nu+1)[-1]
  v_ <- seq(0, 2*pi, length.out = nv+1)[-1]
  for(i in 1:nu){
    for(j in 1:nv){
      vs[,(i-1)*nv+j] <- HT(u_[i], v_[j], nlobes, A)
    }
  }
  if(normals){
    Normals <- matrix(NA_real_, nrow=4, ncol=nu*nv)
    for(i in 1L:nu){
      im1 <- ifelse(i==1L, nu, i-1L)
      ip1 <- ifelse(i==nu, 1L, i+1L)
      for(j in 1L:nv){
        jm1 <- ifelse(j==1L, nv, j-1L)
        jp1 <- ifelse(j==nv, 1L, j+1L)
        n1 <- xprod(vs[,(i-1)*nv+jp1]-vs[,(i-1)*nv+j],
                    vs[,(ip1-1)*nv+j]-vs[,(i-1)*nv+j])
        n2 <- -xprod(vs[,(i-1)*nv+jm1]-vs[,(i-1)*nv+j],
                     vs[,(ip1-1)*nv+jm1]-vs[,(i-1)*nv+j])
        n3 <- xprod(vs[,(im1-1)*nv+j]-vs[,(i-1)*nv+j],
                    vs[,(im1-1)*nv+jp1]-vs[,(i-1)*nv+j])
        n4 <- xprod(vs[,(ip1-1)*nv+j]-vs[,(i-1)*nv+j],
                    vs[,(ip1-1)*nv+jm1]-vs[,(i-1)*nv+j])
        n5 <- -xprod(vs[,(im1-1)*nv+j]-vs[,(i-1)*nv+j],
                     vs[,(i-1)*nv+jm1]-vs[,(i-1)*nv+j])
        n6 <- xprod(vs[,(im1-1)*nv+jp1]-vs[,(i-1)*nv+j],
                    vs[,(i-1)*nv+jp1]-vs[,(i-1)*nv+j])
        Normals[,(i-1)*nv+j] <- c(c(n1+n2+n3+n4+n5+n6)/6,1)
      }
    }
  }
  tris1 <- matrix(NA_integer_, nrow=3, ncol=nu*nv)
  tris2 <- matrix(NA_integer_, nrow=3, ncol=nu*nv)
  for(i in 1L:nu){
    ip1 <- ifelse(i==nu, 1L, i+1L)
    for(j in 1L:nv){
      jp1 <- ifelse(j==nv, 1L, j+1L)
      tris1[,(i-1)*nv+j] <- c((i-1L)*nv+j,(i-1L)*nv+jp1, (ip1-1L)*nv+j)
      tris2[,(i-1)*nv+j] <- c((i-1L)*nv+jp1,(ip1-1L)*nv+jp1,(ip1-1L)*nv+j)
    }
  }
  tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    normals = if(normals) Normals else NULL
  )
}

hopf <- addNormals(hmesh(200, 100, normals = FALSE))
hopf <- Rvcg::vcgIsotropicRemeshing(hopf, TargetLen = 0.5)
points <- t(hopf$vb[-4L, ])
points3d(points)
mesh <- PoissonReconstruction(points, getSomeNormals(points, 10))#, spacing = 0.05, sm_distance = 0.1)
shade3d(mesh, color = "red")
wire3d(mesh, color = "black")


normals <- t(hopf$normals[-4L, ])

indices <- sample.int(40000, 10000)
mesh <- PoissonReconstruction(points[indices, ], normals[indices, ])#spacing = 0.05, sm_distance = 0.1)
shade3d(mesh, color = "red")

set.seed(666)
points <- points[sample.int(40000, 10000), ]

mesh <- PoissonReconstruction(points)#, spacing = 0.05, sm_distance = 0.1)
shade3d(mesh, color = "red")
wire3d(mesh, color = "black")

mesh <- PoissonReconstruction(points, getSomeNormals(points, 10))
shade3d(mesh, color = "red")
wire3d(mesh, color = "black")

mesh <- PoissonReconstruction(points, compute_normals_cpp2(points, 10))
shade3d(mesh, color = "red")
wire3d(mesh, color = "black")

mesh <- AFSreconstruction(points)
shade3d(mesh, color = "red")


##############################################################################
A <- 0.44
n <- 3
Gamma <- function(t){
  alpha <- pi/2 - (pi/2-A)*cos(n*t)
  beta <- t + A*sin(2*n*t)
  c(
    sin(alpha) * cos(beta),
    sin(alpha) * sin(beta),
    cos(alpha)
  )
}
HopfInverse <- function(p, phi){
  c(
    (1+p[3])*cos(phi),
    p[1]*sin(phi) - p[2]*cos(phi),
    p[1]*cos(phi) + p[2]*sin(phi),
    (1+p[3])*sin(phi)
  ) / sqrt(2*(1+p[3]))
}
Stereo <- function(q){
  2*q[1:3] / (1-q[4])
}
F <- function(t, phi){
  Stereo(HopfInverse(Gamma(t), phi))
}

fx <- Vectorize(function(u,v) F(u,v)[1])
fy <- Vectorize(function(u,v) F(u,v)[2])
fz <- Vectorize(function(u,v) F(u,v)[3])
library(misc3d)
tgls <- parametric3d(fx, fy, fz, umin = 0, umax = 2*pi, vmin = 0, vmax = 2*pi,
             n = 200, engine = "none")
