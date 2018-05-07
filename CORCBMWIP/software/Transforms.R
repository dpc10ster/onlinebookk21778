ThetaPrime <- function (Theta)
{
  # the vector of Theta, starting with a, b, zeta2,..., zeta_R-1
  # R is the number of ratings bins; 5 in the example below
  R <- length(Theta) - 1 
  ThetaPrime <- Theta # allocates the output vector   
  ThetaPrime[ 1 ] <- log(-log((Theta[1]-a_min)/(a_max-a_min))) # a to a'
	ThetaPrime[ 2 ] <- log(-log((Theta[2]-b_min)/(b_max-b_min))) # b to b'
  ThetaPrime[ 3 ] <- Theta[3];		# zeta1 to zeta1' no change needed  
  for (r in 4:(R+1)) ThetaPrime[r] <- log(Theta[r] - Theta[r-1])
  return (ThetaPrime)
}
  


Theta <- function (ThetaPrime)
{
  # the vector of transformed ThetaPrime, starting with a', b', zeta2', zeta2',..., zeta_R-1'
  # R is the number of ratings bins; 5 in the example below
  R <- length(ThetaPrime) - 1  
  Theta <- ThetaPrime # allocates the output vector  
  Theta[ 1 ] <- a_min+(a_max-a_min)*exp(-exp(ThetaPrime[1]))# a' --> a
  Theta[ 2 ] <- b_min+(b_max-b_min)*exp(-exp(ThetaPrime[2]))# b' --> b
  Theta[ 3 ] <- ThetaPrime[ 3 ];	# # zeta1' to zeta1 no change needed
  for (r in 4:(R+1)) Theta[r] <- exp(ThetaPrime[r]) + Theta[r-1]
  return (Theta)
}

ForwardZetas <- function(zetas){
  zetasFwd <- zetas
  if (length(zetas) > 1){
    zetasFwd[2:length(zetasFwd)] <- log(zetasFwd[2:length(zetasFwd)] - zetasFwd[1:(length(zetasFwd) - 1)])
  }
  zetasFwd[1] <- ForwardValue(zetasFwd[1], -3, 3)
  zetasFwd[-1] <- ForwardValue(zetasFwd[-1], -7, 3)
  return(zetasFwd)
}

InverseZetas <- function(zetasFwd){
  zetasFwd[1] <- InverseValue(zetasFwd[1], -3, 3)
  zetasFwd[-1] <- InverseValue(zetasFwd[-1], -7, 3)
  zetas <- zetasFwd
  if (length(zetasFwd) > 1) {
    for (i in 2:length(zetasFwd)) zetas[i] <- exp(zetasFwd[i]) + zetas[i - 1]
  }
  return(zetas)
}

ForwardValue <- function(value, valueLower, valueUpper){
  maxValue <- valueUpper
  minValue <- valueLower
  return(log(-log((value - minValue)/(maxValue - minValue))))
}

InverseValue <- function(valueFwd, valueLower, valueUpper){
  maxValue <- valueUpper
  minValue <- valueLower
  return(minValue + (maxValue - minValue) * (exp(-exp(valueFwd))))
}