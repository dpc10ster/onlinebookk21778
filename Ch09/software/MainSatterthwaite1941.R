#this illustrates Satterthewaites' approximation as described in his 1941 paper
# mainSatterthwaite1941.R
rm(list = ls());seed <- 1;set.seed(seed)
library(ggplot2)
sigma1 <- 1;sigma2 <- 4
N1 <- 15;N2 <- 20 # expected to give bad approximation; use N1 <- 15;N2 <- 20 for better fit
cat("sigma1 = ", sigma1, "sigma2 = ", sigma2, "N1 = ", N1, "N2 = ", N2, "\n")
r1 <- N1-1;r2 <- N2-1
# s1Sqd and s2Sqd are what Satterthwaite calls "simple estimates of variance"
s1Sqd <- var(rnorm(N1, sd = sigma1)) # same as sum((xx1-mean(xx1))^2)/(N1-1)
s2Sqd <- var(rnorm(N2, sd = sigma2))
rSat <- (s1Sqd+s2Sqd)^2/(s1Sqd^2/r1+s2Sqd^2/r2) # Satterthwaite approximation
cat("Satterthwaite df = ", rSat, ", approx. accuracy (should be close to unity) = ", 
    r1*s2Sqd/(r2*s1Sqd), "\n") #if second qnty is close to one, approximation should work

S <- 1000
v1 <- array(dim = S);v2 <- array(dim = S)
for (s in 1:S) {
  v1[s] <- var(rnorm(N1, sd = sigma1))# simple estimate of variance
  v2[s] <- var(rnorm(N2, sd = sigma2))# do 
}
scaledS2 <- rSat*(v1 + v2)/(sigma1^2+sigma2^2)# scaled values of v1+v2

x <- scaledS2
qqData <- data.frame(scaledS2 = sort(x), quantile = sort(qchisq(ppoints(v1), df = rSat-1)))
qqPlot <- ggplot(data = qqData, mapping = aes(x = scaledS2, y = quantile)) + geom_point() + geom_abline(slope = 1, color = "red", linetype = 2) + 
  xlab(label = "Scaled S2") + ylab(label = "Chi-square Quantile")
print(qqPlot)
