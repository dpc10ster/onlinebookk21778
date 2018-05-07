F1 <- function(x, y, rho1){
  sigma <- rbind(c(1, rho1), c(rho1, 1))
  return(pmvnorm(c(-Inf, -Inf), c(x, y), sigma = sigma))
}

F2 <- function(x, y, muX, muY, alpha, rho2){
  return((1-alpha) * F1(x, y, rho2) + alpha * F1(x - muX, y - muY, rho2))
}

F1dx <- function(x, y, rho1){
  return(dnorm(x)* pnorm((y - x*rho1)/sqrt(1 - rho1^2)))
}

F1dy <- function(x, y, rho1){
  return(dnorm(y)* pnorm((x - y*rho1)/sqrt(1 - rho1^2)))
}

F1drho1 <- function(x, y, rho1){
  sigma <- rbind(c(1, rho1), c(rho1, 1))
  return(dmvnorm(c(x, y), sigma = sigma))
}

F2dx <- function(x, y, muX, muY, alpha, rho2){
  if (x == -Inf || x == Inf){
    return(0)
  }else{
    d <- (1- alpha) * dnorm(x) * pnorm((y - rho2 * x)/sqrt(1 - rho2^2)) + (1- alpha) * dnorm(x - muX) * pnorm(((y - muY) - rho2 * (x - muX))/sqrt(1 - rho2^2))
  }
  return(d)
}

F2dy <- function(x, y, muX, muY, alpha, rho2){
  if (y == -Inf || y == Inf){
    return(0)
  }else{
    d <- (1- alpha) * dnorm(y) * pnorm((x - rho2 * y)/sqrt(1 - rho2^2)) + (1- alpha) * dnorm(y - muY) * pnorm(((x - muX) - rho2 * (y - muY))/sqrt(1 - rho2^2))
    return(d)
  }
}

F2dmuX <- function(x, y, muX, muY, alpha, rho2){
  if (x == -Inf || x == Inf){
    return(0)
  }else{
    d <- -alpha * dnorm(x - muX) * pnorm((y - muY - rho2 * (x - muX)) / sqrt( 1- rho2^2))
    return(d)
  }
}

F2dmuXdzetar <- function(x, y, muX, muY, alpha, rho2){
  delta <- 1e-7
  d <- (F2dmuX(x + delta, y, muX, muY, alpha, rho2) - F2dmuX(x, y, muX, muY, alpha, rho2)) / delta
  return(d)
}

F2dmuY <- function(x, y, muX, muY, alpha, rho2){
  if (y == -Inf || y == Inf){
    return(0)
  }else{
    d <- -alpha * dnorm(y - muY) * pnorm((x - muX - rho2 * (y - muY)) / sqrt( 1- rho2^2))
    return(d)
  }
}

F2dalpha <- function(x, y, muX, muY, rho2){
  sigma <- rbind(c(1, rho2), c(rho2, 1))
  mu <- c(muX, muY)
  d <- pmvnorm(c(-Inf, -Inf), c(x, y) - mu, sigma = sigma) - pmvnorm(c(-Inf, -Inf), c(x, y), sigma = sigma)
  return(d)
}

F2drho2 <- function(x, y, muX, muY, alpha, rho2){
  if (any(c(x, y) == -Inf) || any(c(x,y) == Inf)){
    return (0)
  }else{
    sigma <- rbind(c(1, rho2), c(rho2, 1))
    mu <- c(muX, muY)
    d <- (1 - alpha) * dmvnorm(c(x, y), sigma = sigma) + alpha * dmvnorm(c(x, y) - mu, sigma = sigma)
    return(d)
  }
}
