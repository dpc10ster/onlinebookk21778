# mainLarsenMarx455.R

# this code illustrates theorem that estimated variance times (N-1)
# divided by true variance follows the chi-square distribution with (N-1) df
# see Larsen and Marx, page 455??
rm(list = ls());seed <- 1;set.seed(seed)
require(stats)
library(ggplot2)
sigma <- 0.4#true value
N <- 10#number of samples over which variance will be estimated

S <- 1000 #number of independent simulations over which variance of variance will be estimated

Var <- array(dim = S)# to hold individual variance estimates
scaledS2 <- array(dim = S)#to hold (N-1)*Var[s]/sigma^2 estimates

for (s in 1:S) {
  Var[s] <- var(rnorm(N, sd =  sigma))# sample variance for this simulation
  scaledS2[s] <- (N-1)*Var[s]/sigma^2
}

# preferred way of testing a model
x <- scaledS2
qqData <- data.frame(x = sort(x), quantile = sort(qchisq(ppoints(x), df = N-1)))
qqPlot <- ggplot(data = qqData, mapping = aes(x = x, y = quantile)) + 
  geom_point() + geom_abline(slope = 1, color = "red", linetype = 2)
print(qqPlot)

# if you insist on using histogram...
pdfData <- data.frame(x = x, y = dchisq(x, df = N-1))
histogram <- ggplot(data = qqData, mapping = aes(x = x)) + geom_histogram(mapping = aes(y = ..density..), color = "black", fill = "grey") + 
  geom_line(mapping = aes(y = y), data = pdfData, color = "red", linetype = 2)
print(histogram)
