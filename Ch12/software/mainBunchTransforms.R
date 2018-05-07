# mainBunchTransforms.R
rm(list=ls())
a <- 1.3204 # see book section 6.2.7.1
b <- 0.6075
mu <- a/b; sigma <- 1/b
zeta = seq(-2,mu+3, 0.1)
FPF <- pnorm(-zeta)
TPF <- pnorm(a-b*zeta)
plot(FPF, TPF, type = "l", xlim=c(0,1), ylim=c(0,1))
LLF <- (TPF-FPF)/(1-FPF)
NLF <- -log(1-FPF)
plot(NLF, LLF, type = "l", xlim=c(0,max(NLF)), ylim=c(0,1))
