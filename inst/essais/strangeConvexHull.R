library(RCGAL)

#' pt <- function(x){
#'   c(
#'     sin(x) * cos(2 * x),
#'     sin(x) * sin(2 * x),
#'     cos(x)
#'   )
#' }
#' pts <- t(vapply(seq(0, pi, length.out = 50), pt, numeric(3L)))
#' hull <- convexhull(pts)
#' plotConvexHull3D(hull, color = "random", luminosity = "dark")
