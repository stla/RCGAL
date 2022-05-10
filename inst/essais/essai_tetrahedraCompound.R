library(RCGAL)
library(rgl)

mesh <- MeshesUnion(tetrahedraCompound, exact = TRUE)
tmesh <- tmesh3d(
  vertices = t(mesh$vertices),
  indices = t(mesh$faces),
  homogeneous = FALSE
)
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.8)
shade3d(tmesh, color = "chartreuse")
plotEdges(mesh$vertices, mesh$exteriorEdges,
          edgesAsTubes = FALSE, verticesAsSpheres = FALSE)

# animation ####
movie3d(spin3d(axis = c(0, 1, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE,
        startTime = 1/10,
        webshot = FALSE)


command <- "gifski --fps=9 --frames=zzpic*.png -o tetrahedraCompound.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)
