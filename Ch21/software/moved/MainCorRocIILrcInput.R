source("initMainCorrocII.R")
options(digits=3)

ret <- ReadLrcFile("Franken1.lrc");
#ret <- ReadLrcFile("VanDyke.lrc")
zk1 <- ret$z1;zk2 <- ret$z2;
Jall <- length(zk1[1,,1]);K1 <- length(zk1[1,1,]);K2 <- length(zk2[1,1,])
jSelected <- 1
NumBins <- 5
retShell <- RunShell (NumBins,1,K1,K2,zk1,zk2)
if (length(retShell) == 1) stop("CorrocII failed") # skip degenerate readers

parametersOrig <- c(retShell$axOrig,retShell$bxOrig,retShell$ayOrig,retShell$byOrig,retShell$rhonOrig,retShell$rhosOrig)
zCovOrig <- retShell$CovOrig
cat("The 6 parameters are ", parametersOrig, "\n")
cat("The 2 sided pValue is ", retShell$pValue, "\n")
cat("The covariance matrix is:\n")
print(zCovOrig)
