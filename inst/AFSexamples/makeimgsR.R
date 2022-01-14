library(rgl)
library(RCGAL)

setwd("C:/SL/MyPackages/RCGAL/inst/AFSexamples")

data(bunny, package = "onion")
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.8)
points3d(bunny)
rgl.snapshot("BunnyCloud.png")

afs <- AFSreconstruction(bunny)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.8)
shade3d(afs, color = "darkred")
rgl.snapshot("BunnyMesh.png")

command <-
  "magick convert BunnyCloud.png BunnyMesh.png +append Bunny.png"
system(command)
