rm( list = ls()) # MainEffectMuOnSd.R
source('GenerateCaseSamples.R')
source('VarJack.R')
source('Wilcoxon.R')
FinalParameters <- c(1.320455, 0.6074974, 0.007676989, 0.8962713, 1.515645, 2.396711)
a  <- FinalParameters[1];b  <-  FinalParameters[2];zetas <- FinalParameters[3:length((FinalParameters))]
mu <- a/b;sigma <- 1/b;K <- c(60, 50)#;mu <- 1.5;sigma <- 1.3
seed <- 1;cat("K1 = ", K[1], ", K2 = ", K[2], ", sigma = ", sigma, "\n")
P1 <- 20
muArr <- mu*seq(0.2,2,length.out = 5)
for (i in 1:length(muArr)) {
  set.seed( seed );VJK<- array(dim = P1); {for (p in 1 : P1) VJK[p] <- VarJack(K, muArr[i], sigma, zetas, bin = FALSE)}
  sdAuc <- mean(sqrt(VJK))
  sdMu <- sqrt(2*sdAuc^2/(dnorm(mu/sqrt(2)))^2)
  cat("mu = ", muArr[i], ", Mean Sd Jack Sampling = ", sdAuc, ", sdMu/mu = ", sdMu/mu, ", sdMu = ", sdMu, "\n")
}

