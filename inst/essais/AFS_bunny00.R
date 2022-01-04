library(RCGAL)
library(rgl)

data(bunny, package = "onion")
mesh <- AFSreconstruction(bunny)
shade3d(mesh, color = "darkred")

pts <- uniformly::runif_on_sphere(500, 3)
mesh <- AFSreconstruction(pts)
shade3d(mesh, color = "red")

tp <- tessellation::teapot()
mesh <- AFSreconstruction(tp)
shade3d(mesh, color = "red")

data(dummyhead, package = "Rvcg")
dummyhead <- t(dummyhead.mesh$vb[-4L, ])
mesh <- AFSreconstruction(dummyhead)
shade3d(mesh, color = "red")

setwd("C:/SL/MyPackages/RCGAL/inst/essais")

tp <- readLines("teapot.obj", n=3644)
writeLines(tp, "teapotVertices.txt")
dat <- read.table("teapotVertices.txt")
pts <- as.matrix(dat[, c(2, 3, 4)])
mesh <- AFSreconstruction(pts)
shade3d(mesh, color = "red")

