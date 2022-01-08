library(misc3d)

a <- 0.9
f <- function(x,y,z){
  (sqrt((x^2 - y^2)^2/(x^2 + y^2) + 3*(z*sin(a))^2) - 3)^2 +
    6*(sqrt((x*y)^2/(x^2 + y^2) + (z*cos(a))^2) - 1.5)^2
}

n <- 250
x <- y <- seq(-4.7, 4.7, length.out = n)
z <- seq(-2.8, 2.8, length.out = n)
g <- expand.grid(x = x, y = y, z = z)
v <- array(with(g, f(x, y, z)), dim = c(n, n, n))

surf <- computeContour3d(v, max(v), 0.5, x = x, y = y, z = z)

set.seed(666L)
SpiderCage <- surf[sample.int(nrow(surf), 10000L), ]
