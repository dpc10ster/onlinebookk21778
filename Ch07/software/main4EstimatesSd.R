rm( list = ls()) # Main4EstimatesSd.R
source('GenerateCaseSamples.R')
source('VarPopSampling.R')
source('VarBootstrap.R')
source('VarJack.R')
source('Wilcoxon.R')
source("VarDeLong.R")
FinalParameters <- c(1.320455, 0.6074974, 0.007676989, 0.8962713, 1.515645, 2.396711)
a  <- FinalParameters[1];b  <-  FinalParameters[2];zetas <- FinalParameters[3:length((FinalParameters))]
mu <- a/b;sigma <- 1/b;K <- c(600, 500)#;mu <- 1.5;sigma <- 1.3
seed <- 1;cat("seed = ", seed, "K1 = ", K[1], ", K2 = ", K[2], ", mu = ", mu, ", sigma = ", sigma, "\n")

P <- 2000;B <- 2000;P1 <- 20

set.seed( seed );VPS<- array(dim = P1); {for (p in 1 : P1) VPS[p] <- VarPopSampling(K, mu, sigma, zetas, P)}
set.seed( seed );VBS<- array(dim = P1); {for (p in 1 : P1) VBS[p] <- VarBootstrap(K, mu, sigma, zetas, B)}
set.seed( seed );VJK<- array(dim = P1); {for (p in 1 : P1) VJK[p] <- VarJack(K, mu, sigma, zetas)} 
set.seed( seed );VDL <- array(dim = P1); {for (p in 1 : P1) VDL[p] <- VarDeLong(K, mu, sigma, zetas)}
cat("Mean Sd Pop Sampling = ", mean(sqrt(VPS)),"\n")
cat("Mean Sd Boot Sampling = ", mean(sqrt(VBS)),"\n")
cat("Mean Sd Jack Sampling = ", mean(sqrt(VJK)),"\n")
cat("Mean Sd DeLong  = ", mean(sqrt(VDL)),"\n")

