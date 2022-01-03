data(bunny, package = "onion")
afs <- RCGAL:::AFSreconstruction(bunny)

pts <- uniformly::runif_on_sphere(1500, 3)
afs <- RCGAL:::AFSreconstruction(pts)

library(rgl)
mesh <- tmesh3d(afs$vertices, afs$triangles, normals=-afs$normals)
