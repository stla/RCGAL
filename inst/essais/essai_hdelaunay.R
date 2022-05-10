library(RCGAL)
library(uniformly)
library(rsvg)

path <- "C:/SL/MyPackages/RCGAL/inst/trash/"

points <- runif_in_sphere(11L, d = 2)
hh <- RCGAL:::htest(t(points))
hedges <- t(hh$edges)
points <- t(hh$vertices)



# svg(filename = file.path(path, "hdelaunay.svg"), onefile = TRUE)
opar <- par(mar = c(0, 0, 0, 0))
plot(
  NULL, type = "n", asp = 1, xlim = c(-1, 1), ylim = c(-1, 1),
  xlab = NA, ylab = NA, axes = FALSE
)
plotrix::draw.circle(0, 0, radius = 1, border = "black")
points(points, pch = 19)
for(i in 1L:nrow(hedges)){
  hedge <- hedges[i, ] + 1L
  hsegment <- RCGAL:::gyrosegment(points[hedge[1L], ], points[hedge[2L], ])
  lines(hsegment, lty = "dashed", col = "red")
}
par(opar)
dev.off()

rsvg_png(
  file.path(path, "hdelaunay.svg"),
  file = file.path(path, "hdelaunay.png")
)

