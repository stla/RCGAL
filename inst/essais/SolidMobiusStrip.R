# http://data.imaginary-exhibition.com/IMAGINARY-Moebiusband-Stephan-Klaus.pdf
# ####

library(RCGAL)
library(rgl)

a = 0.4
b = 0.1
f <- function(x,y,z,a,b){
  ((x*x+y*y+1)*(a*x*x+b*y*y)+z*z*(b*x*x+a*y*y)-2*(a-b)*x*y*z-a*b*(x*x+y*y))^2 -
    4*(x*x+y*y)*(a*x*x+b*y*y-x*y*z*(a-b))^2
}

library(spray) # define f as a polynomial (to get the gradient later)
P <- f(lone(1,3), lone(2,3), lone(3,3), a, b)

# run the marching cubes algorithm ####
nx <- 120; ny <- 120; nz <- 120
x <- seq(-1.4, 1.4, length=nx)
y <- seq(-1.7, 1.7, length=ny)
z <- seq(-0.7, 0.7, length=nz)
g <- expand.grid(x=x, y=y, z=z)
voxel <- array(with(g, f(x,y,z,a,b)), c(nx,ny,nz))
library(misc3d)
surf <- computeContour3d(voxel, maxvol=max(voxel), level=0,
                         x=x, y=y, z=z)


set.seed(666L)
SolidMobiusStrip <- surf[sample.int(nrow(surf), 10000L), ]
mesh <- AFSreconstruction(SolidMobiusStrip)
shade3d(mesh, color = "maroon")





# build the rgl mesh ####
dfx <- as.function(deriv(P, 1))
dfy <- as.function(deriv(P, 2))
dfz <- as.function(deriv(P, 3))
gradient <- function(xyz){
  cbind(dfx(xyz), dfy(xyz), dfz(xyz))
}
library(rgl)
mesh <- tmesh3d(vertices = t(surf),
                indices = matrix(1:nrow(surf), nrow=3),
                homogeneous = FALSE,
                normals = -gradient(surf))
mesh$normals <- rbind(mesh$normals,1)

# plot ####
open3d(windowRect = c(50, 50, 562, 562), zoom= 0.75)
bg3d(rgb(54, 57, 64, maxColorValue = 255))
shade3d(mesh, color=rgb(1,0,1))
