#MainNumericalIntegration.R
rm(list = ls()); set.seed(1); #options(digits=3) 
# numerically integrate pdf of N(0,1) from -infinity to infinty; should give unity
res <- integrate(dnorm, -Inf, Inf)
cat("integral from -inf to inf = ",res$value, "\n")
# numerically integrate the pdf of N(0,1) from 1.96 to infinty; should give 0.025
res <- integrate(dnorm, 1.96, Inf)
cat("integral from 1.96 to inf = ",res$value, "\n")

