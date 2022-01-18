library(rgl)
library(misc3d)

a <- 2.4

f <- function(x, y, z, a){
  (sqrt((x^2 - a^2)^2/(x^2 + a^2) + y^2) - 3)^2 +
    (sqrt((x*a)^2/(x^2 + a^2) + z^2) - 1.5)^2
}

# run the marching cubes algorithm ####
nx <- 200; ny <- 200; nz <- 200
x <- seq(-5, 5, length.out = nx)
y <- seq(-4, 4, length.out = ny)
z <- seq(-2.5, 2.5, length.out = nz)
G <- expand.grid(x = x, y = y, z = z)
voxel <- array(with(G, f(x, y, z, a)), dim = c(nx, ny, nz))
surface <- computeContour3d(
  voxel, maxvol = max(voxel), level = 0.5, x = x, y = y, z = z
)

mesh0 <- misc3d:::t2ve(makeTriangles(surface))
mesh <- addNormals(tmesh3d(
  vertices = mesh0$vb,
  indices = mesh0$ib
))

open3d(windowRect = c(50, 50, 562, 562))
bg3d(rgb(54, 57, 64, maxColorValue = 255))
view3d(30, -30, zoom= 0.75)
shade3d(mesh, color = "tomato")

# animation ####
M <- par3d("userMatrix")
movie3d(
  par3dinterp(
    time = seq(0, 1, len = 9),
    userMatrix = list(
      M,
      rotate3d(M, pi, 1, 0, 0),
      rotate3d(M, pi, 1, 1, 0),
      rotate3d(M, pi, 1, 1, 1),
      rotate3d(M, pi, 0, 1, 1),
      rotate3d(M, pi, 0, 1, 0),
      rotate3d(M, pi, 1, 0, 1),
      rotate3d(M, pi, 0, 0, 1),
      M
    )
  ),
  fps = 120,
  duration = 1,
  dir = ".",
  frames = "zzzpic",
  convert = "echo \"%d %s %s %s\"",
  clean = FALSE,
  webshot = FALSE
)

pngs <- list.files(".", pattern = "^zzzpic", full.names = TRUE)
library(gifski)
gifski(pngs, "ICN5D_eight.gif",
       width = 512, height = 512, delay = 1/9)

file.remove(pngs)


# https://thumbs.gfycat.com/MajorCarefreeHuman-size_restricted.gif
