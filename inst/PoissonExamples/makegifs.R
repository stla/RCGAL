library(rgl)
library(RCGAL)

setwd("C:/SL/MyPackages/RCGAL/inst/makegifs")

open3d(windowRect = c(50, 50, 450, 450))
view3d(0, -50, zoom = 0.7)
points3d(SpiderCage)
rgl.snapshot("SpiderCageCloud.png")

psr <- PoissonReconstruction(SpiderCage)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, -50, zoom = 0.7)
shade3d(psr, color = "orangered")
wire3d(psr, color = "black")
rgl.snapshot("SpiderCageMesh.png")

command <-
  "magick convert SpiderCageCloud.png SpiderCageMesh.png +append SpiderCage.png"
system(command)


open3d(windowRect = c(50, 50, 450, 450))
view3d(0, -50, zoom = 0.7)
points3d(SolidMobiusStrip)
rgl.snapshot("SolidMobiusStripCloud.png")

psr <- PoissonReconstruction(SolidMobiusStrip)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, -50, zoom = 0.7)
shade3d(psr, color = "orangered")
wire3d(psr, color = "black")
rgl.snapshot("SolidMobiusStripMesh.png")

command <-
  "magick convert SolidMobiusStripCloud.png SolidMobiusStripMesh.png +append SolidMobiusStrip.png"
system(command)


open3d(windowRect = c(50, 50, 450, 450))
view3d(0, -50, zoom = 0.7)
points3d(ToroidalHelix)
rgl.snapshot("ToroidalHelixCloud.png")

psr <- PoissonReconstruction(ToroidalHelix)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, -50, zoom = 0.7)
shade3d(psr, color = "orangered")
wire3d(psr, color = "black")
rgl.snapshot("ToroidalHelixMesh.png")

command <-
  "magick convert ToroidalHelixCloud.png ToroidalHelixMesh.png +append ToroidalHelix.png"
system(command)


open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.7)
points3d(HopfTorus)
rgl.snapshot("HopfTorusCloud.png")

psr <- PoissonReconstruction(HopfTorus)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.7)
shade3d(psr, color = "orangered")
wire3d(psr, color = "black")
rgl.snapshot("HopfTorusMesh.png")

command <-
  "magick convert HopfTorusCloud.png HopfTorusMesh.png +append HopfTorus.png"
system(command)

psr <- PoissonReconstruction(HopfTorus, spacing = 0.2)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.7)
shade3d(psr, color = "orangered")
wire3d(psr, color = "black")
rgl.snapshot("HopfTorusMesh_spacing02.png")
