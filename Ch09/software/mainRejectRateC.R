# mainRejectRateC.R
# Simulating testing of analysis method; CLUSTER version
rm( list = ls(all = TRUE))
library(RJafroc)
library(foreach)
library(doParallel)
source( 'RoeMetzZSamples.R' )

nomAlpha <- 0.05;seed <- NULL;set.seed( seed )
I <- 2;J <- 5;K1 <- 52;K2 <- 50;VarStrString <- "LH1"
VarStr <- RoeMetzVarStr( VarStrString )
cat("Var Str = ", unlist(VarStr), "\n")
mu2 <- sqrt(2)*qnorm(VarStr$auc);tau22 <- 0.0*mu2 # zero multipler sets NH condition
eSize <- pnorm((mu2+tau22)/sqrt(2)) - pnorm((mu2)/sqrt(2)) # effect size in AUC units
isBinned <- TRUE;NumBins <- 5
failed <- 0
for (rep in 1:1) {
  S <- 2000
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  pVal <- foreach(s = 1 : S, .combine = "c", .packages = c("stringr", "RJafroc")) %dopar% {
    # following line creates z-samples 
    zSamplesInit <- InitRandomSamples( I, J, K1, K2 ) # assuming all components of VarStr are unity
    # following line includes the effect of the variance components and bins data, if needed 
    zSamplesRaw <- RoeMetzZSamples( 
      I, J, K1, K2, mu2, tau22, VarStr, zSamplesInit)
    NL <- array(-Inf,dim=c(I,J,(K1+K2),1));NL[,,1:K1,1] <- zSamplesRaw$FP
    LL <- array(-Inf,dim=c(I,J,K2,1));LL[,,1:K2,1] <- zSamplesRaw$TP
    dataset <- Df2RJafrocDataset(NL, LL)
    if (isBinned) {
      dataset <- DfBinDataset(dataset,desiredNumBins = NumBins)
    }
    ret <- StSignificanceTesting(dataset, fom = "Wilcoxon", option = "RRRC")
    ret$pRRRC
  }
  stopCluster(cl)
  empAlpha <- sum(pVal < nomAlpha )/S
  ciWidth <- qnorm(1-nomAlpha/2)*sqrt(empAlpha*(1-empAlpha)/S)
  CI <- c(empAlpha-ciWidth, empAlpha+ciWidth)
  cat("Effect size = ", eSize, ", NH rejection rate over", S,  "simulations: \n")
  cat("95% CI lower = ", CI[1] , "\nempAlpha = ", empAlpha, "\n95% CI upper = ", CI[2], "\n")
  if ((abs(eSize) < 0.001) & (!((CI[1] < nomAlpha) & (CI[2] > nomAlpha)))) {
    cat("failed NH test.\n")
    failed <- failed + 1
  }
}
cat("total number of failed = ", failed, "\n")