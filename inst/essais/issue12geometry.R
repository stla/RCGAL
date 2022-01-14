setwd("C:/SL/MyPackages/RCGAL/inst/essais")

library(RCGAL)

X <- as.matrix(read.table("points.txt"))

d2 <- delaunay(X)

plotDelaunay2D(d2)

plot(X, pch = 19, cex = 0.5)
for(i in 1L:nrow(d2$edges)){
  edge <- d2$edges[i, ]
  segments(X[edge[1], 1], X[edge[1], 2], X[edge[2], 1], X[edge[2], 2])
}
