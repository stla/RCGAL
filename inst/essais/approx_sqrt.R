# sqrt(2)
u <- function(n){
  dbls <- (0L:(n-1L))
  prod(1L/2 - dbls) / factorial(n) / (-50)^n
}
10/7 * (1 - 1/100 + sum(sapply(2:10, u)))

# sqrt(5)
u <- function(n){
  dbls <- (0L:(n-1L))
  prod(1L/2 - dbls) / factorial(n) / (-5)^n
}
10/4 * (1 - 1/10  + sum(sapply(2:10, u)))

library(gmp)

u <- function(n){
  dbls <- 2*(0L:(n-1L))
  as.bigq(10*prod(1L - dbls), factorialZ(n) * (as.bigz(-10))^(n))
}
asNumeric(as.bigq(1L, 4L) * (as.bigq(9L)  + Reduce("+", sapply(2:10, u), 0L)))
sqrt(5)

as.bigq(1L, 2L) * # phi
  (1L + as.bigq(1L, 4L) * (as.bigq(9L)  + Reduce("+", sapply(2:10, u), 0L)))

# sqrt(3)
u <- function(n){
  dbls <- 2*(0L:(n-1L))
  prod(1L - dbls) / factorial(n) / 2^n/ (-4)^n
}
10/5 * (1 - 5/40  + sum(sapply(2:10, u)))

library(gmp)

u <- function(n){
  dbls <- 2*(0L:(n-1L))
  as.bigq(8*prod(1L - dbls), factorialZ(n) * (as.bigz(-8))^(n))
}
asNumeric(as.bigq(1L, 4L) * (as.bigq(7L)  + Reduce("+", sapply(2:11, u), 0L)))


# u(1L) + u(2L) + u(3L) + u(4L) + u(5L) + u(6L) + u(7L)
# 1 + sum(sapply(1:30, u))


`%^%` <- function(A, n){
  Reduce(`%*%`, replicate(n, A, simplify = FALSE))
}

library(gmp)

qsqrt <- function(x, n){
  # stopifnot(isPositiveInteger(x) || x == 0)
  # stopifnot(isPositiveInteger(n))
  zero <- as.bigz(0L)
  one <- as.bigz(1L)
  A <- matrix(c(zero, as.bigz(x)-1L, one, as.bigz(2L)), nrow = 2L, ncol = 2L)
  zs <- c((A %^% n) %*% c(zero, one))
  out <- as.bigq(zs[2L], zs[1L]) - 1L
  class(out) <- c("qsqrt", class(out))
  attr(out, "error") <- abs(asNumeric(out) - sqrt(x))
  out
}

qsqrt(3, 7)

print.qsqrt <- function(x, ...){
  print(x)
  cat(attr(x, "error"))
}

print.qsqrt(qsqrt(3, 7))
