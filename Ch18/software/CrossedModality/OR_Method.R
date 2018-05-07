OR_Method <- function (I, J, K, cov)
{

  FOM <- cov$FOM
  Var <- cov$Var
  Cov1 <- cov$Cov1
  Cov2 <- cov$Cov2
  Cov3 <- cov$Cov3
  Cov2i <- cov$Cov2i
  Vari <- cov$Vari
  
  fom_mean <- mean(FOM)
  
  MST <- 0
  for (i in 1:I){
    MST <- MST + (mean(FOM[i, ]) - fom_mean)^2
  }
  MST <- J * MST / (I - 1)

  MSTR <- 0
  for (i in 1:I){
    for (j in 1:J){
      MSTR <- MSTR + (FOM[i, j] - mean(FOM[i, ]) - mean(FOM[, j]) + fom_mean)^2
    }
  }
  MSTR <- MSTR / ((J - 1) * (I - 1))
  
  # do random reader random case analysis first
  # do random reader random case analysis first
  # do random reader random case analysis first
  MS_DEN_DIFF_FOM_RRRC <- MSTR + max(J * (Cov2 - Cov3), 0)
  ndf <- (I-1)
  ddf_DIFF_FOM_RRRC <- MS_DEN_DIFF_FOM_RRRC^2 / (MSTR^2 / ((J - 1) * (I - 1)))
  F_DIFF_FOM_RRRC <- MST / MS_DEN_DIFF_FOM_RRRC
  
  pVal_DIFF_FOM_RRRC <- (1- pf(F_DIFF_FOM_RRRC, ndf, ddf_DIFF_FOM_RRRC)) 
  if (pVal_DIFF_FOM_RRRC < alpha) reject_P <- 1 else reject_P <- 0
  
  modality_IDs <- array("",dim = c(I))
  for (i in 1:I) {
    modality_IDs[i] <- paste(toString(i)," ")  
  }

  modality_means <- array(dim = I)
  for (i in 1:I) modality_means[i] <- mean(FOM[i,])
  
  modality_means_differences <- array(dim = c(I,I))
  modality_IDs_differences <- array("",dim = c(I,I))
  for (i1 in 1:(I-1)) {    
    for (i2 in (i1+1):I) {
      modality_IDs_differences[i1,i2] <- paste(toString( i2), "-", toString(i1), sep = "")
      modality_means_differences[i1,i2] <- modality_means[i2]- modality_means[i1]    
    }
  }
  
  temp <- modality_means_differences
  modality_means_differences <- modality_means_differences[!is.na(temp)]
  modality_IDs_differences <- modality_IDs_differences[!is.na(temp)]
  nDiffs <-  length(modality_means_differences)
  
  std_DIFF_FOM_RRRC <- sqrt(2*MS_DEN_DIFF_FOM_RRRC/J)
  CI_DIFF_FOM_RRRC <- array(dim = c(nDiffs, 3))
  for (i in 1 : nDiffs) {
    CI_DIFF_FOM_RRRC[i,1] <- modality_means_differences[i]
    CI_DIFF_FOM_RRRC[i,2] <- qt(alpha/2,df = ddf_DIFF_FOM_RRRC)*std_DIFF_FOM_RRRC + modality_means_differences[i]
    CI_DIFF_FOM_RRRC[i,3] <- qt(1-alpha/2,df = ddf_DIFF_FOM_RRRC)*std_DIFF_FOM_RRRC + modality_means_differences[i]
  }
  
  CutPt <- qf(1 - alpha, ndf, ddf_DIFF_FOM_RRRC)  
  if (F_DIFF_FOM_RRRC > CutPt) reject_F <- 1 else reject_F <- 0  
  if (reject_F != reject_P) stop("F ne P")
  
  # do fixed reader random case analysis next
  # do fixed reader random case analysis next
  # do fixed reader random case analysis next
  MS_DEN_DIFF_FOM_FRRC <- Var - Cov1 + (J-1)*(Cov2-Cov3)
  F_DIFF_FOM_FRRC <- MST / MS_DEN_DIFF_FOM_FRRC 
  ddf_DIFF_FOM_FRRC <- (I-1)*(K-1) # why not infinity?
  pVal_DIFF_FOM_FRRC <- (1- pf(F_DIFF_FOM_FRRC, ndf, ddf_DIFF_FOM_FRRC))
  #pVal_DIFF_FOM_FRRC <- (1- pchisq(F_DIFF_FOM_FRRC, ndf)) # why do they not use what is in the literature?
  
  std_DIFF_FOM_FRRC <-  sqrt(2*MS_DEN_DIFF_FOM_FRRC/J)
  CI_DIFF_FOM_FRRC <- array(dim = c(nDiffs, 3))
  for (i in 1 : nDiffs) {
    CI_DIFF_FOM_FRRC[i,1] <- modality_means_differences[i]
    CI_DIFF_FOM_FRRC[i,2] <- qt(alpha/2,df = ddf_DIFF_FOM_FRRC)*std_DIFF_FOM_FRRC + modality_means_differences[i]
    CI_DIFF_FOM_FRRC[i,3] <- qt(1-alpha/2,df = ddf_DIFF_FOM_FRRC)*std_DIFF_FOM_FRRC + modality_means_differences[i]
  }
  
  # do random reader fixed case analysis last
  # do random reader fixed case analysis last
  # do random reader fixed case analysis last
  F_DIFF_FOM_RRFC <- MST/MSTR
  ddf_DIFF_FOM_RRFC <- (I-1)*((J-1))
  pVal_DIFF_FOM_RRFC <- (1- pf(F_DIFF_FOM_RRFC, ndf, ddf_DIFF_FOM_RRFC))
  
  std_DIFF_FOM_RRFC <- sqrt(2*MSTR/J)
  CI_DIFF_FOM_RRFC <- array(dim = c(nDiffs, 3))
  for (i in 1 : nDiffs) {
    CI_DIFF_FOM_RRFC[i,1] <- modality_means_differences[i]
    CI_DIFF_FOM_RRFC[i,2] <- qt(alpha/2,df = ddf_DIFF_FOM_RRFC)*std_DIFF_FOM_RRFC + modality_means_differences[i]
    CI_DIFF_FOM_RRFC[i,3] <- qt(1-alpha/2,df = ddf_DIFF_FOM_RRFC)*std_DIFF_FOM_RRFC + modality_means_differences[i]
  }
  
  # confidence intervals for reader-averaged FOM in each modality
  # confidence intervals for reader-averaged FOM in each modality
  # confidence intervals for reader-averaged FOM in each modality
  MSR_i <- array(0, dim = I)
  for (i in 1:I) {
    for (j in 1:J) {
      MSR_i[i] <- MSR_i[i] + (FOM[i, j] - mean(FOM[i, ]))^2                              
    }
  }
  MSR_i <- MSR_i/(J-1)
  
  MS_DEN_DIFF_FOM_RRRCSingle_i <- array(dim = I)
  for (i in 1:I){
    MS_DEN_DIFF_FOM_RRRCSingle_i[i] <- MSR_i[i] + max(J*Cov2i[i],0)
  }
  
  df_AVG_FOM_RRRC <- (MS_DEN_DIFF_FOM_RRRCSingle_i/MSR_i)^2*(J-1)
  CI_AVG_FOM_RRRC <- array(dim = c(I,3))
  std_AVG_FOM_RRRC <- sqrt(MS_DEN_DIFF_FOM_RRRCSingle_i/J)
  for (i in 1:I){
    CI_AVG_FOM_RRRC[i,1] <- mean(FOM[i, ])
    CI_AVG_FOM_RRRC[i,2] <- qt(alpha/2,df = df_AVG_FOM_RRRC[i]) * std_AVG_FOM_RRRC[i] + mean(FOM[i, ]) #works
    CI_AVG_FOM_RRRC[i,3] <- qt(1-alpha/2,df = df_AVG_FOM_RRRC[i]) * std_AVG_FOM_RRRC[i] + mean(FOM[i, ]) #works
  }
  
  CI_AVG_FOM_FRRC <- array(dim = c(I,3))
  df_AVG_FOM_FRRC <- (K-1)
  std_AVG_FOM_FRRC <- sqrt((Vari + (J-1)*Cov2i)/J) # fixed by reverse engineering, intuition, from other results
  for (i in 1:I){
    CI_AVG_FOM_FRRC[i,1] <- mean(FOM[i, ])
    CI_AVG_FOM_FRRC[i,2] <- qt(alpha/2,df = df_AVG_FOM_FRRC) * std_AVG_FOM_FRRC[i] + mean(FOM[i, ]) #(almost) works for 1st modality
    CI_AVG_FOM_FRRC[i,3] <- qt(1-alpha/2,df = df_AVG_FOM_FRRC) * std_AVG_FOM_FRRC[i] + mean(FOM[i, ]) #off by about 1% for 2nd modality
  }
  # "PET-CT (10AR).xlsx" data
  # std dev for modality 2 should be 0.01835119
  # I am getting 0.01592989
  
  
  CI_AVG_FOM_RRFC <- array(dim = c(I,3))
  df_AVG_FOM_RRFC <- (J-1)
  std_AVG_FOM_RRFC <- sqrt(MSR_i/J)
  for (i in 1:I){
    CI_AVG_FOM_RRFC[i,1] <- mean(FOM[i, ])
     CI_AVG_FOM_RRFC[i,2] <- qt(alpha/2,df = df_AVG_FOM_RRFC) * std_AVG_FOM_RRFC[i] + mean(FOM[i, ]) #works
     CI_AVG_FOM_RRFC[i,3] <- qt(1-alpha/2,df = df_AVG_FOM_RRFC) * std_AVG_FOM_RRFC[i] + mean(FOM[i, ]) #works
  }
  
  rownames(FOM) <- modality_IDs
  colnames(CI_DIFF_FOM_RRRC) <- c("Mean", "Lower", "Upper") 
  rownames(CI_DIFF_FOM_RRRC) <- modality_IDs_differences 
  colnames(CI_DIFF_FOM_FRRC) <- c("Mean", "Lower", "Upper") 
  rownames(CI_DIFF_FOM_FRRC) <- modality_IDs_differences 
  colnames(CI_DIFF_FOM_RRFC) <- c("Mean", "Lower", "Upper") 
  rownames(CI_DIFF_FOM_RRFC) <- modality_IDs_differences 
  colnames(CI_AVG_FOM_RRRC) <- c("Mean", "Lower", "Upper") 
  rownames(CI_AVG_FOM_RRRC) <- modality_IDs 
  colnames(CI_AVG_FOM_FRRC) <- c("Mean", "Lower", "Upper") 
  rownames(CI_AVG_FOM_FRRC) <- modality_IDs 
  colnames(CI_AVG_FOM_RRFC) <- c("Mean", "Lower", "Upper") 
  rownames(CI_AVG_FOM_RRFC) <- modality_IDs 

  return (list(
    FOM = FOM,
    Var = Var,
    Cov1 = Cov1,
    Cov2 = Cov2,
    Cov3 = Cov3,
    Vari = Vari,
    Cov2i = Cov2i,
    ddf_DIFF_FOM_RRRC = ddf_DIFF_FOM_RRRC,
    ddf_DIFF_FOM_FRRC = ddf_DIFF_FOM_FRRC,    
    ddf_DIFF_FOM_RRFC = ddf_DIFF_FOM_RRFC,
    F_DIFF_FOM_RRRC = F_DIFF_FOM_RRRC,
    F_DIFF_FOM_FRRC = F_DIFF_FOM_FRRC,
    F_DIFF_FOM_RRFC = F_DIFF_FOM_RRFC,
    pVal_DIFF_FOM_RRRC = pVal_DIFF_FOM_RRRC,
    pVal_DIFF_FOM_FRRC = pVal_DIFF_FOM_FRRC,
    pVal_DIFF_FOM_RRFC = pVal_DIFF_FOM_RRFC,
    CI_DIFF_FOM_RRRC = CI_DIFF_FOM_RRRC,
    CI_DIFF_FOM_FRRC = CI_DIFF_FOM_FRRC,
    CI_DIFF_FOM_RRFC = CI_DIFF_FOM_RRFC,
    df_AVG_FOM_RRRC = df_AVG_FOM_RRRC,
    df_AVG_FOM_FRRC = df_AVG_FOM_FRRC,
    df_AVG_FOM_RRFC = df_AVG_FOM_RRFC,
    CI_AVG_FOM_RRRC = CI_AVG_FOM_RRRC,
    CI_AVG_FOM_FRRC = CI_AVG_FOM_FRRC,
    CI_AVG_FOM_RRFC = CI_AVG_FOM_RRFC
  ))
}