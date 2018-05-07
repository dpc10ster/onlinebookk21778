deePrimeToAz <- function(deePrime) {
  return(pnorm(deePrime/sqrt(2)))
}

azToDeePrime <- function(Az) {
  return(sqrt(2)*qnorm(Az))
}

effectSizeDeePrime <- function(effectSizeAz, Az)
{
  deePrime1 <-  azToDeePrime (Az)
  temp <- (Az + effectSizeAz) < 1.0
  Az <- Az[temp]
  deePrime1 <-  deePrime1[temp]
  deePrime2 <-  azToDeePrime (Az + effectSizeAz) 
  return (deePrime2-deePrime1)
}

