VarCov2 <- function (zjk1, zjk2, flag, FPFValue) # for difference FOM, radiologist minus CAD
{ 
  J <- length(zjk1[,1]) - 1 # number of human readers
  K1 <- length(zjk1[1,]);K2 <- length(zjk2[1,]);K <- K1 +  K2
  
  PsiJk = array(dim = c(J,K))  
  for (k in 1:K) {
    if (k <= K1) {
      zjk1_jk <- zjk1[,-k] 
      zjk2_jk <- zjk2
    } else {
      zjk1_jk <- zjk1 
      zjk2_jk <- zjk2[,-(k-K1)]
    }
    
    if (flag == "AUC") {
      aucCad_jk <- Wilcoxon(zjk1_jk[1, ], zjk2_jk[1, ])
      for (j in 1:J) PsiJk[j, k] <- Wilcoxon(zjk1_jk[j + 1, ], zjk2_jk[j + 1, ]) - aucCad_jk
    } else{
      temp <- array (dim = (J+1))
      if (flag == "PCL") {
        PCLCad_jk <- (LrocFoms(zjk1_jk[1,], zjk2_jk[1,], FPFValue))$PCL  # the array zjk2Cl is actually passed 
        for (j in 1:J)  PsiJk[j,k] <- (LrocFoms(zjk1_jk[j + 1,], zjk2_jk[j + 1,], FPFValue))$PCL - PCLCad_jk # do
      } else if (flag == "ALROC") {
        ALrocCad_jk <- (LrocFoms(zjk1_jk[1,], zjk2_jk[1,], FPFValue))$ALroc  # the array zjk2Cl is actually passed 
        for (j in 1:J)  PsiJk[j,k] <- (LrocFoms(zjk1_jk[j + 1,], zjk2_jk[j + 1,], FPFValue))$ALroc - ALrocCad_jk #do
      } else if (flag == "JAFROC") {
        JAFROCCad_jk <- Wilcoxon(zjk1_jk[1, ], zjk2_jk[1, ]) # the array zjk2Cl is actually passed
        for (j in 1:J) PsiJk[j, k] <- Wilcoxon(zjk1_jk[j + 1, ], zjk2_jk[j + 1, ]) - JAFROCCad_jk # do:
      } else stop("wrong flag value")  
    }
  } 
  
  Covariance <- array(dim = c(J, J))  
  for (j in 1:J){
    for (jp in 1:J){
      Covariance[j, jp] <- cov(PsiJk[j, ], PsiJk[jp, ])          
    }
  }  
  
  Var <- 0;count <- 0
  for (j in 1:J){    
    Var <- Var + Covariance[j, j] 
    count <- count + 1
  }
  Var <- Var / count 
  Var <- Var *(K-1)^2/K
  
  Cov2 <- 0;count <- 0
  for (j in 1:J){    
    for (jp in 1:J){
      if (jp != j){
        Cov2 <- Cov2 + Covariance[j, jp] 
        count <- count + 1
      }
    }
  }  
  Cov2 <- Cov2 / count 
  Cov2 <- Cov2 *(K-1)^2/K
  
  return (list (    
    Var = Var,
    Cov2 = Cov2
  ))  
}


DiffFomAnal2007Hillis53 <- function (zjk1, zjk2, flag, FPFValue, alpha = 0.05)
{  
  ret <- VarCov2(zjk1, zjk2, flag, FPFValue)
  Var <- ret$Var;  Cov2 <- ret$Cov2  
  J <- length(zjk1[,1]) - 1 
  
  thetajc <- array (dim = (J+1))
  if (flag == "AUC") {
    for (j in 1:(J+1)) thetajc[j] <- Wilcoxon(zjk1[j,],zjk2[j,])
  } else if (flag == "PCL") {
    for (j in 1:(J+1)) thetajc[j] <- (LrocFoms(zjk1[j,], zjk2[j,], FPFValue))$PCL # the array zjk2Cl is actually passed 
  } else if (flag == "ALROC") {
    for (j in 1:(J+1)) thetajc[j] <- (LrocFoms(zjk1[j,], zjk2[j,], FPFValue))$ALroc # do:
  } else if (flag == "JAFROC") {
    for (j in 1:(J+1)) thetajc[j] <- Wilcoxon(zjk1[j,], zjk2[j,]) # the array zjk2Cl is actually passed 
  } else stop("wrong flag value")  

  Psijc <- thetajc[2:(J+1)] - thetajc[1]
  
  MSR <- 0
  PsiMean <- mean(Psijc)
  for (j in 1:J){
    MSR <- MSR + (Psijc[j] - PsiMean)^2
  }
  MSR <- MSR / (J - 1)
  
  MSdenOR_single <- MSR + max(J * Cov2, 0)
  DdfhSingle <- MSdenOR_single^2 / (MSR^2 / (J - 1))
  TstatStar <- PsiMean / sqrt(MSdenOR_single/J)
  p_val <- 1 - pt(abs(TstatStar),DdfhSingle) + pt(-abs(TstatStar),DdfhSingle)
  if (p_val < alpha) reject_P <- 1 else reject_P <- 0
  
  CutPt <- qt(1 - alpha/2, DdfhSingle)  
  if (abs(TstatStar) > CutPt) reject_Tstat <- 1 else reject_Tstat <- 0 
  
  CI <- array(dim=2)
  CI[1] <- qt(alpha/2,df = DdfhSingle)
  CI[2] <- qt(1-alpha/2,df = DdfhSingle)
  CI <- CI * sqrt(MSdenOR_single/J)
  CI <- CI + PsiMean
  if (CI[1] * CI[2] > 0) reject_CI <- 1 else reject_CI <- 0
  
  if (reject_Tstat != reject_P) stop("Tstat ne P")
  if (reject_Tstat != reject_CI) stop("Tstat ne CI")
  
  return (list (
    Var = Var, 
    Cov2 = Cov2, 
    PsiMean = PsiMean,
    CI = CI,
    reject = reject_Tstat,
    ddfH = DdfhSingle,
    Tstat = TstatStar,
    p_val = p_val,
    thetajc = thetajc
  ))
}
