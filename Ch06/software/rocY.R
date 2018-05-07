# rocY.R
# x is rocX, i.e., FPF
rocY <- function (x, a, b) {
  y <- pnorm(a + b*qnorm(x))
  return(y)
}
