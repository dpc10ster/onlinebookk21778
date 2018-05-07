CbmRocY <- function (x, mu, alpha) {
  y <- (1-alpha)*(1-pnorm(qnorm(1-x))) + alpha*(1-pnorm(qnorm(1-x)-mu))
  return(y)
}
