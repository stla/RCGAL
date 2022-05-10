library(RCGAL)
library(uniformly)

points <- runif_in_sphere(10L, d = 2)
hedges <- t(RCGAL:::htest(t(points)))

plot(
  NULL, type = "n", asp = 1, xlim = c(-1, 1), ylim = c(-1, 1),
  xlab = NA, ylab = NA, axes = FALSE
)
plotrix::draw.circle(0, 0, radius = 1, border = "black")
points(points, pch = 19)
for(i in 1L:nrow(hedges)){
  hedge <- hedges[i, ]
  hsegment <- RCGAL:::gyrosegment(points[hedge[1L], ], points[hedge[2L], ])
  lines(hsegment, lty = "dashed", col = "red")
}

