# first fit URSM, find AUC, then impose this as constraint on RSM. two fits shown: URSM and RRSM
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
fileName <- fileName[fileName == "FZR"]
cat("fileName = ", fileName,"\n")
frocData <- loadDataFile(fileName)
I <- length(frocData$modalityID)
J <- length(frocData$readerID)
lesionNum <- frocData$lesionNum
nLesDistr <- table(lesionNum)
if (length(nLesDistr) == 1) {
  nLesDistr <- c(lesionNum[1], 1)
  dim(nLesDistr) <- c(1, 2)
}else{
  nLesDistr <- t(rbind(as.numeric(unlist(attr(nLesDistr, "dimnames"))), as.vector(nLesDistr)))
}

retFileName <- paste0("ANALYZED/URSM-AUC/", "saveRetRoc", fileName)
if (!file.exists((retFileName))){
  s <- 1
  retSmRoc <- list()
  for (i in 1:I){
    for (j in 1:J){
      cat("i = ", i, ", j = ", j, "\n")
      retURSM <- FitURsmRoc(frocData, i, j)
      tempRoc <- FitRsmRoc(frocData, i, j, constraint = "URSM")
      retSmRoc[[s]] <- as.list(c(tempRoc,
                                 list(nLesDistr = nLesDistr, URSMAUC = retURSM$AUC, URSMmu = retURSM$mu, URSMnuP = retURSM$nuP, URSMcutoffs = retURSM$cutoffs, i = i, j = j)))
      s <- s + 1
    }
  }
  save(retSmRoc, file = retFileName)
}else{
  load(retFileName)
}

muRsm <- rep(NA, I*J)
lambdaRsm <- rep(NA, I*J)
lambdaPRsm <- rep(NA, I*J)
nuRsm <- rep(NA, I*J)
nuPRsm <- rep(NA, I*J)
muURsm <- rep(NA, I*J)
lambdaURsm <- rep(NA, I*J)
lambdaPURsm <- rep(NA, I*J)
nuURsm <- rep(NA, I*J)
nuPURsm <- rep(NA, I*J)
S <- rep(NA, I*J)
AucC <- rep(NA, I*J);C <- rep(NA, I*J)

s <- 1
for (i in 1:I){
  for (j in 1:J){
    empOp <- EmpiricalOpCharac(frocData, i, j, opChType = "ROC")$ROCPoints
    muURsm[s] <- as.numeric(retSmRoc[[s]]$URSMmu$mu)
    lambdaPURsm[s] <- as.numeric(retSmRoc[[s]]$URSMlambdaP$lambdaP)
    nuPURsm[s] <- as.numeric(retSmRoc[[s]]$URSMnuP$nuP)
    lambdaURsm[s] <- as.numeric(retSmRoc[[s]]$URSMlambda$lambda)
    nuURsm[s] <- as.numeric(retSmRoc[[s]]$URSMnu$nu)
    
    muRsm[s] <- as.numeric(retSmRoc[[s]]$mu$mu)
    lambdaPRsm[s] <- as.numeric(retSmRoc[[s]]$lambdaP$lambdaP)
    nuPRsm[s] <- as.numeric(retSmRoc[[s]]$nuP$nuP)
    lambdaRsm[s] <- as.numeric(retSmRoc[[s]]$lambda$lambda)
    nuRsm[s] <- as.numeric(retSmRoc[[s]]$nu$nu)
    nLesDistr <- retSmRoc[[s]]$nLesDistr
    fpf <- empOp$FPF; tpf <- empOp$TPF
    fpf <- fpf[-c(1, length(fpf))]
    tpf <- tpf[-c(1, length(tpf))]
    
    curveComp <- PlotRsmOperatingCharacteristics(c(muRsm[s], muURsm[s]), c(lambdaRsm[s], lambdaURsm[s]), c(nuRsm[s], nuURsm[s]), nLesDistr, type = "ROC")$ROCPlot
    opPnts <- data.frame(FPF = fpf, TPF = tpf)
    curveComp <- curveComp + geom_point(data = opPnts, mapping = aes(x = FPF, y = TPF)) + labs(title = paste("i =", i, ", j =", j)) +
      scale_color_discrete(labels=c("RSM","URSM")) + theme(legend.title=element_blank())
    print(curveComp)
    S[s] <- nuPRsm[s] * exp(-lambdaPRsm[s])
    retSmRoc[[s]]$S <- S[s]
    C[s] <- muRsm[s]
    AucC[s] <- pnorm(muRsm[s]/sqrt(2))
    retSmRoc[[s]]$C <- AucC[s]
    cat(i, j, ", mu = ", muRsm[s], ", lambdaP = ", lambdaPRsm[s], ", nuP = ", nuPRsm[s],
        ", S = ", retSmRoc[[s]]$S, ", AucC = ", retSmRoc[[s]]$C,
        ", muURsm = ", muURsm[s],  ", lambdaPURsm = ", lambdaPURsm[s],  ", nuPURsm = ", nuPURsm[s], ", AUC = ", as.numeric(retSmRoc[[s]]$AUC$AUC),
        "\n")
    s <- s + 1
  }
}

ndigits <- 5
df <- data.frame(muRsm = muRsm, muURsm = muURsm)
p <- ggplot(data = df, aes(x = muRsm, y = muURsm)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  geom_point()
print(p)
m <- lm(muURsm ~ muRsm, df)
cat("a = ", coef(m)[1],", b = ", coef(m)[2],", r2 = ", summary(m)$r.squared,"\n")

df <- data.frame(nuPRsm = nuPRsm, nuPURsm = nuPURsm)
p <- ggplot(data = df, aes(x = nuPRsm, y = nuPURsm)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  geom_point()
print(p)
m <- lm(nuPURsm ~ nuPRsm, df);
cat("a = ", coef(m)[1],", b = ", coef(m)[2],", r2 = ", summary(m)$r.squared,"\n")

df <- data.frame(S = S, AucC = AucC)
p <- ggplot(data = df, aes(x = S, y = AucC)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  geom_point()
print(p)
m <- lm(AucC ~ S, df);
cat("a = ", coef(m)[1],", b = ", coef(m)[2],", r2 = ", summary(m)$r.squared,"\n")
