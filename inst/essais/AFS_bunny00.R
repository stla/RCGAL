library(RCGAL)
library(rgl)

data(bunny, package = "onion")
mesh <- AFSreconstruction(bunny)

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

