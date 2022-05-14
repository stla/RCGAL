# the compound of five tetrahedra ####
library(gmp) # we use rational numbers for the vertex coordinates
library(RCGAL)
library(rgl)

# all vertices
Vertices <- function(a, b, c, phi, O){
  t(cbind(
    c( a,  a,  a),
    c( a,  a, -a),
    c( a, -a,  a),
    c(-a, -a,  a),
    c(-a,  a, -a),
    c(-a,  a,  a),
    c( O,  b, -c),
    c( O, -b, -c),
    c( O, -b,  c),
    c( c,  O, -b),
    c(-c,  O, -b),
    c(-c,  O,  b),
    c( b,  c,  O),
    c( b, -c,  O),
    c(-b, -c,  O),
    c(-b,  c,  O),
    c( O,  b,  c),
    c( a, -a, -a),
    c( c,  O,  b),
    c(-a, -a, -a)
  ))
}

# 'double' vertex coordinates
phi <- (1 + sqrt(5)) / 2
a <- 1 / sqrt(3)
b <- a / phi
c <- a * phi
O <- 0
vertices <- Vertices(a, b, c, phi, O)

# 'gmp' vertex coordinates
phi <- qsqrtPhi(17)
a <- 1L / qsqrt3(22)
b <- a / phi
c <- a * phi
O <- as.bigq(0)
gmpVertices <- Vertices(a, b, c, phi, O)
# check the difference is close to 0:
max(abs(asNumeric(gmpVertices) - vertices))

# the face indices are common to all tetrahedra
faces <- rbind(
  c(1L, 2L, 3L),
  c(3L, 2L, 4L),
  c(4L, 2L, 1L),
  c(1L, 3L, 4L)
)

# define the five tetrahedra meshes
mesh1 <- list(
  "vertices" = vertices[c(17, 14, 2, 11), ],
  "faces" = faces
)
mesh2 <- list(
  "vertices" = vertices[c(18, 1, 4, 5), ],
  "faces" = faces
)
mesh3 <- list(
  "vertices" = vertices[c(19, 6, 15, 7), ],
  "faces" = faces
)
mesh4 <- list(
  "vertices" = vertices[c(3, 13, 12, 8), ],
  "faces" = faces
)
mesh5 <- list(
  "vertices" = vertices[c(20, 16, 10, 9), ],
  "faces" = faces
)

# put them in a list to apply the intersection algorithm
meshes <- list(
  mesh1, mesh2, mesh3, mesh4, mesh5
)

tetrahedraCompound <- list(

  "meshes" = lapply(meshes, function(mesh){
    list("vertices" = mesh[["vertices"]], faces = faces)
  }),

  "rglmeshes" = list(

    tmesh3d(
      "vertices"    = t(mesh1[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(mesh2[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(mesh3[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(mesh4[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(mesh5[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    )
  )
)

library(RCGAL)
library(rgl)
# plot
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.75)
bg3d("#363940")
#
palette <- trekcolors::trek_pal("klingon")
colors <- colorRampPalette(palette)(5)
invisible(lapply(1:5, function(i){
  shade3d(tetrahedraCompound[["rglmeshes"]][[i]], color = colors[i])
}))

# animation ####
movie3d(spin3d(axis = c(0, 1, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE,
        startTime = 1/10,
        webshot = FALSE)


command <- "gifski --fps=10 --frames=zzpic*.png -o tetrahedraCompound.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)


stop("")

# compute the intersection of the 'gmp' meshes
inter <- MeshesIntersection(
  tetrahedraCompound[["meshes"]], numbersType = "lazy", clean = TRUE
)

# plot
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.75)
bg3d("#363940")
# first the five tetrahedra with transparency
invisible(lapply(
  tetrahedraCompound[["rglmeshes"]], shade3d,
  color = "yellow", alpha = 0.1
))
# now the intersection
rglinter <- tmesh3d(
  "vertices"    = t(inter[["vertices"]]),
  "indices"     = t(inter[["faces"]]),
  "homogeneous" = FALSE
)
shade3d(rglinter, color = "gainsboro")
# and finally the edges
plotEdges(
  inter[["vertices"]], inter[["exteriorEdges"]],
  only = inter[["exteriorVertices"]], color = "darkmagenta"
)
# this is an icosahedron

# animation ####
movie3d(spin3d(axis = c(0, 1, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE,
        startTime = 1/10,
        webshot = FALSE)


command <- "gifski --fps=10 --frames=zzpic*.png -o tetrahedraCompound.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)


# we can recognize some exact values in the 'gmpVertices', for instance
#   the value at entry \code{(1, 1)} is \code{1/sqrt(3)}:
asNumeric(1L / gmpVertices[1L, 1L]^2L)


library(microbenchmark)

tetrahedraCompound <- list(

  "double_meshes" = lapply(meshes, function(mesh){
    list("vertices" = mesh[["vertices"]], faces = faces)
  }),

  "gmp_meshes" = lapply(meshes, function(mesh){
    list("vertices" = as.bigq(mesh[["vertices"]]), faces = faces)
  })

)

microbenchmark(
  exactKernel = MeshesIntersection(
    tetrahedraCompound[["double_meshes"]], numbersType = "lazyExact"
  ),
  gmpKernel = MeshesIntersection(
    tetrahedraCompound[["gmp_meshes"]], numbersType = "gmp"
  ),
  times = 5
)

rglmeshes <- tetrahedraCompound[["rglmeshes"]]
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.75)
bg3d("#363940")
colors <- hcl.colors(5, palette = "Spectral")
invisible(lapply(
  1:5, function(i) shade3d(rglmeshes[[i]], color = colors[i])
))
