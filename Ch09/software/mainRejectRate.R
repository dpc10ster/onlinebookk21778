# mainRejectRate.R
# Simulating testing of analysis method
rm( list = ls( all = TRUE ) )
library(RJafroc)
source( 'RoeMetzZSamples.R' )

nomAlpha <- 0.05;seed <- 1;set.seed( seed )
I <- 2;J <- 5;K1 <- 52;K2 <- 50;VarStrString <- "LH1"
VarStr <- RoeMetzVarStr( VarStrString )
cat("RMVarStr = ", unlist(VarStr), "\n")
mu2 <- sqrt(2)*qnorm(VarStr$auc);tau22 <- 0.0*mu2 # zero multipler sets NH condition
isBinned <- TRUE;NumBins <- 5 # data binning, R = 5
S <- 2000;reject <- 0
for (s in 1 : S) {
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
  ret <- StSignificanceTesting(dataset, FOM = "Wilcoxon", option = "RRRC")
  if (ret$pRRRC < nomAlpha) reject <- reject + 1
}
empAlpha <- reject/S
ciWidth <- qnorm(1-nomAlpha/2)*sqrt(nomAlpha*(1-nomAlpha)/S) # assuming true value is nomAlpha 
CI <- c(nomAlpha-ciWidth, nomAlpha+ciWidth)
cat("NH rejection rate over", S,  "simulations: \n")
cat("95% CI lower = ", CI[1] , "\nempAlpha = ", empAlpha, "\n95% CI upper = ", CI[2], "\n")
if (!((CI[1] < nomAlpha) & (CI[2] > nomAlpha))) cat("failed NH test.\n")
