library(RCGAL)
library(uniformly)
library(rsvg)

path <- "C:/SL/MyPackages/RCGAL/inst/trash/"

points <- runif_in_sphere(11L, d = 2)
hh <- RCGAL:::htest(t(points))
hedges <- t(hh$edges)
hfaces <- t(hh$faces)
points <- t(hh$vertices)



svg(filename = file.path(path, "hdelaunay.svg"), onefile = TRUE)
opar <- par(mar = c(0, 0, 0, 0))
plot(
  NULL, type = "n", asp = 1, xlim = c(-1, 1), ylim = c(-1, 1),
  xlab = NA, ylab = NA, axes = FALSE
)
plotrix::draw.circle(0, 0, radius = 1, border = "black")
points(points, pch = 19, cex = 0.9)
colors <- hcl.colors(nrow(hfaces))
for(i in 1L:nrow(hfaces)){
  hface <- hfaces[i, ]
  hpolygon <- rbind(
    RCGAL:::gyrosegment(points[hface[1L], ], points[hface[3L], ])[-1L, ],
    RCGAL:::gyrosegment(points[hface[3L], ], points[hface[2L], ])[-1L, ],
    RCGAL:::gyrosegment(points[hface[2L], ], points[hface[1L], ])[-1L, ]
  )
  polygon(hpolygon, border = NA, col = colors[i])
}
for(i in 1L:nrow(hedges)){
  hedge <- hedges[i, ]
  hsegment <- RCGAL:::gyrosegment(points[hedge[1L], ], points[hedge[2L], ])
  lines(hsegment, lty = "solid", col = "black", lwd = 2)
}
par(opar)
dev.off()

rsvg_png(
  file.path(path, "hdelaunay.svg"),
  file = file.path(path, "hdelaunay.png")
)

