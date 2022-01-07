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

# mesh: astrotoroidal twisted tubular helix
scos <- function(x,alpha) sign(cos(x)) * abs(cos(x))^alpha
ssin <- function(x,alpha) sign(sin(x)) * abs(sin(x))^alpha
TubularHelixMesh <- function(R, r, w, a, nu, nv, alpha=1, twists=2){
  vs <- matrix(NA_real_, nrow=3, ncol=nu*nv)
  colors <- matrix(NA_character_, nrow = 3L, ncol = nu*nv)
  u_ <- seq(0, w*2*pi, length.out = nu+1)[-1]
  v_ <- seq(0, 2*pi, length.out = nv+1)[-1]
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
      colors[,(i-1)*nv+j] <-
        ifelse(findInterval(v, pi*(1/16+seq(0,2,len=17))) %% 2 == 0,
               "#440154FF", "#FDE725FF")
    }
  }
  tris1 <- matrix(NA_integer_, nrow=3, ncol=nu*nv)
  tris2 <- matrix(NA_integer_, nrow=3, ncol=nu*nv)
  nv <- as.integer(nv)
  for(i in 1L:nu){
    ip1 <- ifelse(i==nu, 1L, i+1L)
    for(j in 1L:nv){
      jp1 <- ifelse(j==nv, 1L, j+1L)
      tris1[,(i-1)*nv+j] <- c((i-1L)*nv+j,(i-1L)*nv+jp1, (ip1-1L)*nv+j)
      tris2[,(i-1)*nv+j] <- c((i-1L)*nv+jp1,(ip1-1L)*nv+jp1,(ip1-1L)*nv+j)
    }
  }
  out <- tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    homogeneous = FALSE,
    material = list(color = colors)
  )
  addNormals(out)
}

# draw ####
m <- TubularHelixMesh(R=5, r=1.35, w=10, a=0.7, 10*20, 30, alpha=1, twists=1)
open3d(windowRect = c(50,50,550,550))
bg3d(rgb(54,57,64, maxColorValue = 255))
view3d(0,-55)
wire3d(m, color="black")

pts <- t(m$vb[-4L,])
mesh <- PoissonReconstruction(pts)
open3d(windowRect = c(50,50,562,562))
view3d(0,-55)
shade3d(mesh, color = "yellow")
wire3d(mesh, color="black")

mesh <- PoissonReconstruction(pts, getSomeNormals(pts, 10))
open3d(windowRect = c(50,50,562,562))
view3d(0,-55)
shade3d(mesh, color = "orange")
wire3d(mesh, color="black")

# draw ####
m <- TubularHelixMesh(R=4, r=1, w=20, a=0.25, 20*60, 60, alpha=0.5, twists=2)
open3d(windowRect = c(50,50,550,550))
bg3d(rgb(54,57,64, maxColorValue = 255))
view3d(0,-55)
shade3d(m)

# animation ####
movie3d(spin3d(axis = c(0, 0, 1), rpm = 60),
        duration = 1, fps = 60,
        movie = "anim00", dir = ".",
        convert = "convert -dispose previous -loop 0 -delay 1x%d %s*.png %s.%s",
        startTime = 1/60)
