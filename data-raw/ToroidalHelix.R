library(rgl)

# helix curve
helix <- function(t, R, r, w){
  c(
    (R + r*cos(t)) * cos(t/w),
    (R + r*cos(t)) * sin(t/w),
    r*sin(t)
  )
}
# derivative (tangent)
dhelix <- function(t, R, r, w){
  v <- c(
    -r*sin(t)*cos(t/w) - (R+r*cos(t))/w*sin(t/w),
    -r*sin(t)*sin(t/w) + (R+r*cos(t))/w*cos(t/w),
    r*cos(t)
  )
  v / sqrt(c(crossprod(v)))
}
# second derivative (normal)
ddhelix <- function(t, R, r, w){
  v <- c(
    -r*cos(t)*cos(t/w) + r*sin(t)/w*sin(t/w) +
      r*sin(t)/w*sin(t/w) - (R+r*cos(t))/w^2*cos(t/w),
    -r*cos(t)*sin(t/w) - r*sin(t)/w*cos(t/w) -
      r*sin(t)/w*cos(t/w) - (R+r*cos(t))/w^2*sin(t/w),
    -r*sin(t)
  )
  v / sqrt(c(crossprod(v)))
}
# binormal
bnrml <- function(t, R, r, w){
  v <- rgl:::xprod(dhelix(t, R, r, w), ddhelix(t, R, r, w))
  v / sqrt(c(crossprod(v)))
}

# mesh maker
scos <- function(x,alpha) sign(cos(x)) * abs(cos(x))^alpha
ssin <- function(x,alpha) sign(sin(x)) * abs(sin(x))^alpha
TubularHelixMesh <- function(R, r, w, a, nu, nv, alpha = 1, twists = 2){
  vs <- matrix(NA_real_, nrow = 3L, ncol = nu*nv)
  u_ <- seq(0, w*2*pi, length.out = nu+1)[-1L]
  v_ <- seq(0, 2*pi, length.out = nv+1)[-1L]
  for(i in 1:nu){
    u <- u_[i]
    for(j in 1:nv){
      v <- v_[j]
      h <- helix(u, R, r, w)
      vs[,(i-1)*nv+j] <-
        h +
        a*(scos(v,alpha) *
             (cos(twists*u)*ddhelix(u, R, r, w) +
                sin(twists*u)*bnrml(u, R, r, w)) +
             ssin(v,alpha) *
             (-sin(twists*u)*ddhelix(u, R, r, w) +
                cos(twists*u)*bnrml(u, R, r, w)))
    }
  }
  tris1 <- matrix(NA_integer_, nrow=3, ncol=nu*nv)
  tris2 <- matrix(NA_integer_, nrow=3, ncol=nu*nv)
  nv <- as.integer(nv)
  for(i in 1L:nu){
    ip1 <- ifelse(i==nu, 1L, i+1L)
    for(j in 1L:nv){
      jp1 <- ifelse(j==nv, 1L, j+1L)
      tris1[,(i-1L)*nv+j] <- c((i-1L)*nv+j,(i-1L)*nv+jp1, (ip1-1L)*nv+j)
      tris2[,(i-1L)*nv+j] <- c((i-1L)*nv+jp1,(ip1-1L)*nv+jp1,(ip1-1L)*nv+j)
    }
  }
  out <- tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    homogeneous = FALSE
  )
  addNormals(out)
}


mesh <- TubularHelixMesh(
  R = 5, r = 1.35, w = 10, a = 0.7, 10*20, 30, alpha = 1, twists = 1
)

ToroidalHelix <- t(mesh$vb[-4L, ])
