EffectSize <- function (muNH, sigmaNH, muAH, sigmaAH)
{
  
  ES <- pnorm(muAH/sqrt(1+sigmaAH^2)) - pnorm(muNH/sqrt(1+sigmaNH^2))
  return (ES)
}