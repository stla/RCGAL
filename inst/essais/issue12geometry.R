setwd("C:/SL/MyPackages/RCGAL/inst/essais")

X <- as.matrix(read.table("points.txt"))

cxh <- cxhull2d(X)

d2 <- del2d(X)

plot(X, pch = 19, cex = 0.5)
for(i in 1L:nrow(d2$edges)){
  edge <- d2$edges[i,]
  segments(X[edge[1],1], X[edge[1],2], X[edge[2],1], X[edge[2],2])
}
