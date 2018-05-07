ForwardTransform <- function (parameters)
{
  # the vector of parameters, starting with a, b, zeta1, zeta2,...
  # R is the number of ratings bins; 5 in the example below
  R <- length(parameters) - 1 
  parameters_prime <- parameters # allocates the output vector   
  Y_SAVE <- parameters
  parameters_prime[ 1 ] <- log(-log((Y_SAVE[1]-a_min)/(a_max-a_min))) # a to a'
	parameters_prime[ 2 ] <- log(-log((Y_SAVE[2]-b_min)/(b_max-b_min))) # b to b'
  parameters_prime[ 3 ] <- Y_SAVE[3];		# zeta1 to zeta1' no change needed  
  for (r in 4:(R+1)) parameters_prime[r] <- log(Y_SAVE[r] - Y_SAVE[r-1])
  return (parameters_prime)
}
  


InverseTransform <- function (parameters_prime)
{
  # the vector of transformed parameters parameters_prime, starting with a', b', zeta1', zeta2',...
  # R is the number of ratings bins; 5 in the example below
  R <- length(parameters_prime) - 1  
  parameters <- parameters_prime  
  Z_SAVE <- parameters_prime 
  parameters[ 1 ] <- a_min+(a_max-a_min)*exp(-exp(Z_SAVE[1]))# a' --> a
  parameters[ 2 ] <- b_min+(b_max-b_min)*exp(-exp(Z_SAVE[2]))# b' --> b
  parameters[ 3 ] <- Z_SAVE[ 3 ];	# # zeta1' to zeta1 no change needed
  for (r in 4:(R+1)) parameters[r] <- exp(Z_SAVE[r]) + parameters[r-1]
  return (parameters)
}

 