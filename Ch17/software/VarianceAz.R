VarianceAz <- function (a, b,Cov)
{
  derivWrtA <- dnorm (a/sqrt(1+b^2))/sqrt(1+b^2)
  derivWrtB <- dnorm (a/sqrt(1+b^2))*(-a*b*(1+b^2)^(-1.5)) 
  VarAz <- (derivWrtA)^2*Cov[1,1]+(derivWrtB)^2*Cov[2,2] +
    2 * derivWrtA*derivWrtB*Cov[1,2]
  return (VarAz)  
}

# using values from Eng website makes little difference
#   VarAz <- (derivWrtA)^2*0.0656 +(derivWrtB)^2*0.0254 +
#   2 * derivWrtA*derivWrtB*0.0259
