dat <- read.table("rat.txt")
vertices <- subset(dat, V1 == "v")[, c(2, 3, 4)]
normals <- subset(dat, V1 == "vn")[, c(2, 3, 4)]

psr <- PoissonReconstruction(
  as.matrix(vertices), as.matrix(normals), spacing = 0.08,
  sm_angle = 50, sm_radius = 3, sm_distance = 0.1
)

shade3d(psr, color = "gray")
wire3d(psr)

afs <- AFSreconstruction(as.matrix(vertices))
