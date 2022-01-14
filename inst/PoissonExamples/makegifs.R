library(rgl)
library(RCGAL)

setwd("C:/SL/MyPackages/RCGAL/inst/PoissonExamples")

data(bunny, package = "onion")
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.8)
points3d(bunny)
rgl.snapshot("StanfordBunnyCloud.png")

psr <- PoissonReconstruction(bunny, spacing = 0.0005)
open3d(windowRect = c(50, 50, 450, 450))
view3d(0, 0, zoom = 0.8)
shade3d(psr, color = "yellowgreen")
wire3d(psr, color = "black")
rgl.snapshot("StanfordBunnyMesh.png")

command <-
  "magick convert StanfordBunnyCloud.png StanfordBunnyMesh.png +append StanfordBunny.png"
system(command)


open3d(windowRect = c(50, 50, 450, 450))
view3d(10, -60, zoom = 0.9)
points3d(CostaSurface)
rgl.snapshot("CostaSurfaceCloud.png")

psr <- PoissonReconstruction(CostaSurface)
open3d(windowRect = c(50, 50, 450, 450))
view3d(10, -60, zoom = 0.9)
shade3d(psr, color = "orangered")
wire3d(psr, color = "black")
rgl.snapshot("CostaSurfaceMesh.png")

command <-
  "magick convert CostaSurfaceCloud.png CostaSurfaceMesh.png +append CostaSurface.png"
system(command)



open3d(windowRect = c(50, 50, 450, 450))
view3d(40, -80, zoom = 0.7)
par3d(userMatrix =
        matrix(c(0.9, -0.5, 0.06, 0,
                 0.36, 0.72, 0.6, 0,
                 -0.33, -0.5, 0.8, 0,
                 0, 0, 0, 1), nrow = 4))
points3d(cyclide)
rgl.snapshot("cyclideCloud.png")

psr <- PoissonReconstruction(cyclide)
open3d(windowRect = c(50, 50, 450, 450))
view3d(40, -80, zoom = 0.7)
par3d(userMatrix =
        matrix(c(0.9, -0.5, 0.06, 0,
                 0.36, 0.72, 0.6, 0,
                 -0.33, -0.5, 0.8, 0,
                 0, 0, 0, 1), nrow = 4))
shade3d(psr, color = "orangered")
wire3d(psr, color = "black")
rgl.snapshot("cyclideMesh.png")

command <-
  "magick convert cyclideCloud.png cyclideMesh.png +append cyclide.png"
system(command)


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
