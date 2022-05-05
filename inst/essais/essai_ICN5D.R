library(rgl)
library(Rvcg)
library(rmarchingcubes)

f <- function(x, y, z, a, cosb, sinb){
    (sqrt((sqrt(x*x + (y*sinb + a*cosb)^2) - 2)^2) - 1)^2 +
      (sqrt((sqrt(z*z + (y*cosb - a*sinb)^2) - 2)^2) - 1)^2
}

a <- 0.6
b <- 0.785
cosb <- cos(b)
sinb <- sin(b)

x <- z <- seq(-3.5, 3.5, len = 150L)
y <- seq(-4.2, 4.2, len = 150L)
g <- expand.grid(X = x, Y = y, Z = z)

voxel <- array(
  with(g, f(X, Y, Z, a, cosb, sinb)),
  dim = c(150L, 150L, 150L)
)

contour_shape <- contour3d(
  griddata = voxel,
  level = 0.1,
  x = x,
  y = y,
  z = z
)

tmesh <- tmesh3d(
  vertices = t(contour_shape[["vertices"]]),
  indices = t(contour_shape[["triangles"]]),
  normals = contour_shape[["normals"]],
  homogeneous = FALSE
)

open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
#view3d(55, -80)
shade3d(tmesh, color = "orangered")

edges <- as.matrix(vcgGetEdge(tmesh)[, c("vert1", "vert2")])

library(cc)
cc <- cc_cpp(edges)

ncc <- attr(cc, "ncomponents")
colors <- hcl.colors(ncc)

open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
for(i in 1L:ncc){
  comp_i <- cc[2L, ] == i
  mesh <- tmesh3d(
    vertices = tmesh$vb,#[, comp_i],
    normals = t(tmesh$normals),
    indices = tmesh$it[, tmesh$it[1L, ] %in% which(comp_i)]
  )
  shade3d(mesh, color = colors[i])
}

# animation ####
movie3d(spin3d(axis = c(0, 1, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE,
        startTime = 1/10,
        webshot = FALSE)


command <- "gifski --fps=9 --frames=zzpic*.png -o ICN5D_viridis.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)
