library(rgl)
library(RCGAL)

setwd("C:/SL/MyPackages/RCGAL/inst/AFSexamples")


open3d(windowRect = c(50, 50, 450, 450))
view3d(-30, 0, zoom = 0.75)
points3d(StanfordDragon)
rgl.snapshot("StanfordDragonCloud.png")

afs <- AFSreconstruction(StanfordDragon)
open3d(windowRect = c(50, 50, 450, 450))
view3d(-30, 0, zoom = 0.75)
shade3d(afs, color = "forestgreen")
rgl.snapshot("StanfordDragonMesh.png")

command <-
  "magick convert StanfordDragonCloud.png StanfordDragonMesh.png +append StanfordDragon.png"
system(command)



data(dummyhead, package = "Rvcg")
dummyhead <- t(dummyhead.mesh$vb[-4L, ])
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.8)
points3d(dummyhead)
rgl.snapshot("DummyHeadCloud.png")

afs <- AFSreconstruction(dummyhead)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.8)
shade3d(afs, color = "darksalmon")
rgl.snapshot("DummyHeadMesh.png")

command <-
  "magick convert DummyHeadCloud.png DummyHeadMesh.png +append DummyHead.png"
system(command)




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
