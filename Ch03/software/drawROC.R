drawROC <- function (mu)
{
  zeta <- seq(-3, mu + 3, 0.01)
  lines(pnorm(-zeta),pnorm(mu-zeta))
  text(pnorm(-mu/2),pnorm(mu/2), paste ("mu =", toString(mu)))
}