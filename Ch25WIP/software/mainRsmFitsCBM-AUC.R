rm(list = ls()) # mainRsmFitsCBM.R
library(RJafroc)
library(bbmle)
library(stats)
library(ggplot2)
library(binom)
source("RsmFunctions.R")
source("Transforms.R")
source("optFunctions.R")
source("addArguments.R")
source("CBMFitR.R")
source("LLRoc.R")
source("FitRsmRocCurve.R")
source("PlotCBMRSM.R")
source("RJafrocIncludes.R")

stop("fix or delete me")

# included datasets
fileName <-  c("TONY", "VD", "FR", "FED", "JT", "MAG", "OPT", "PEN", "NICO", 
               "RUS", "DOB1", "DOB2", "DOB3", "FZR")
fileName <- fileName[fileName == "MAG"]
cat("fileName = ", fileName,"\n")
frocData <- loadDataFile(fileName)

lesionNum <- frocData$lesionNum
nLesDistr <- table(lesionNum)
if (length(nLesDistr) == 1) {
  nLesDistr <- c(lesionNum[1], 1)
  dim(nLesDistr) <- c(1, 2)
}else{
  nLesDistr <- t(rbind(as.numeric(unlist(attr(nLesDistr, "dimnames"))), as.vector(nLesDistr)))
}
rocData <- FROC2HrROC(frocData)
I <- length(rocData$modalityID)
J <- length(rocData$readerID)

retFileName <- paste0("ANALYZED/CBM-AUC/", "saveRetRoc", fileName)
if (!file.exists((retFileName))){
  s <- 1
  retSmRoc <- list()
  for (i in 1:I){
    for (j in 1:J){
      cat("i = ", i, ", j = ", j, "\n")
      CBM <- CBMFitR(rocData, i, j)
      tempRoc <- FitRsmRocCurve(rocData, i, j, AUCCbm = CBM$AUC, zetaCbm = CBM$cutoffs, nLesDistr = nLesDistr)
      retSmRoc[[s]] <- as.list(c(tempRoc,
                                 list(nLesDistr = nLesDistr, CBMAUC = CBM$AUC, CBMmu = CBM$mu, CBMalpha = CBM$alpha, CBMcutoffs = CBM$cutoffs, i = i, j = j)))
      s <- s + 1
    }
  }
  save(retSmRoc, file = retFileName)
}else{
  load(retFileName)
}

muRsm <- rep(NA, I*J)
lambdaP <- rep(NA, I*J)
nuP <- rep(NA, I*J)
muCbm <- rep(NA, I*J)
alpha <- rep(NA, I*J)
S <- rep(NA, I*J)
AucC <- rep(NA, I*J);C <- rep(NA, I*J)

s <- 1
for (i in 1:I){
  for (j in 1:J){
    empOp <- EmpiricalOpCharac(rocData, i, j, opChType = "ROC")$ROCPoints
    muCbm[s] <- retSmRoc[[s]]$CBMmu
    alpha[s] <- retSmRoc[[s]]$CBMalpha
    muRsm[s] <- as.numeric(retSmRoc[[s]]$mu$mu)
    lambdaP[s] <- as.numeric(retSmRoc[[s]]$lambdaP$lambdaP)
    nuP[s] <- as.numeric(retSmRoc[[s]]$nuP$nuP)
    nLesDistr <- retSmRoc[[s]]$nLesDistr
    fpf <- empOp$FPF; tpf <- empOp$TPF
    fpf <- fpf[-c(1, length(fpf))]
    tpf <- tpf[-c(1, length(tpf))]
    curveComp <- PlotCBMRSM(muCbm[s], alpha[s], muRsm[s], lambdaP[s], nuP[s], nLesDistr, fpf, tpf, i, j)
    print(curveComp)
    S[s] <- nuP[s] * exp(-lambdaP[s])
    retSmRoc[[s]]$S <- S[s]
    C[s] <- muRsm[s]
    AucC[s] <- pnorm(muRsm[s]/sqrt(2))
    retSmRoc[[s]]$C <- AucC[s]
    cat(i, j, ", mu = ", muRsm[s], ", lambdaP = ", lambdaP[s], ", nuP = ", nuP[s],
        ", S = ", retSmRoc[[s]]$S, ", AucC = ", retSmRoc[[s]]$C,
        ", muCbm = ", muCbm[s],  ", alpha = ", alpha[s],  ", AUC = ", as.numeric(retSmRoc[[s]]$AUC$AUC),
        "\n")
    s <- s + 1
  }
}

ndigits <- 5
df <- data.frame(muRsm = muRsm, muCbm = muCbm)
p <- ggplot(data = df, aes(x = muRsm, y = muCbm)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  geom_point()
print(p)
m <- lm(muCbm ~ muRsm, df)
cat("a = ", coef(m)[1],", b = ", coef(m)[2],", r2 = ", summary(m)$r.squared,"\n")

df <- data.frame(nuP = nuP, alpha = alpha)
p <- ggplot(data = df, aes(x = nuP, y = alpha)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  geom_point()
print(p)
m <- lm(alpha ~ nuP, df);
cat("a = ", coef(m)[1],", b = ", coef(m)[2],", r2 = ", summary(m)$r.squared,"\n")

df <- data.frame(S = S, AucC = AucC)
p <- ggplot(data = df, aes(x = S, y = AucC)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  geom_point()
print(p)
m <- lm(AucC ~ S, df);
cat("a = ", coef(m)[1],", b = ", coef(m)[2],", r2 = ", summary(m)$r.squared,"\n")
