#MainSampleSize.R
# source("MainAnalyzeData.R") # first we need to analyze the data both ways

# inferred  ROC sample size calculation
effectSize <- 0.05 # rounding down ResultsRoc$ciDiffTrtRRFC$Estimate = 0.0523
varCompOR <- ResultsRoc$varComp
varTR <- varCompOR$varCov[2]
cov1 <- varCompOR$varCov[3]
cov2 <- varCompOR$varCov[4]
cov3 <- varCompOR$varCov[5]
varEps <- varCompOR$varCov[6]  # diagonal element of covariance matrix
KStar <- dim(rocData$NL)[3]
ret <- SampleSizeGivenJ(J = 10, cov1 = cov1, cov2 = cov2, cov3 = cov3, 
                        varTR = varTR, varEps= varEps, KStar = KStar, effectSize =effectSize, 
                        option = "RRFC")
print(ret)

# wAFROC sample size calculation
effectSize <- 0.091283695 # TBA
varCompOR <- ResultsFroc$varComp
varTR <- varCompOR$varCov[2]
cov1 <- varCompOR$varCov[3]
cov2 <- varCompOR$varCov[4]
cov3 <- varCompOR$varCov[5]
varEps <- varCompOR$varCov[6] # diagonal element of covariance matrix
KStar <- dim(rocData$NL)[3]
ret <- SampleSizeGivenJ(J = 10, cov1 = cov1, cov2 = cov2, cov3 = cov3, 
                        varTR = varTR, varEps= varEps, KStar = KStar, effectSize =effectSize, 
                        option = "RRFC")
print(ret)

