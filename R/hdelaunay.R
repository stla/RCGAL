dotprod <- function(x, y = NULL){
  c(crossprod(x, y))
}

gyroadd <- function(X, Y, s = 1){
  x <- dotprod(X) / s^2
  y <- dotprod(Y) / s^2
  xy <- 2 * dotprod(X, Y) / s^2
  ((1 + xy + y) * X + (1 - x) * Y) / (1 + xy + x*y)
}

gyroscalar <- function(r, X, s = 1){
  Xnorm <- sqrt(dotprod(X))
  s / Xnorm * tanh(r * atanh(Xnorm / s)) * X
}

gyrovector <- function(A, B, s = 1){
  gyroadd(-A, B, s = s)
}

gyroABt <- function(A, B, t, s = 1){
  gyroadd(A, gyroscalar(t, gyrovector(A, B, s = s), s = s), s = s)
}

gyrosegment <- function(A, B, n = 50, s = 1){
  t(vapply(seq(0, 1, length.out = n), function(t){
    gyroABt(A, B, t, s = s)
  }, numeric(3L)))
}
