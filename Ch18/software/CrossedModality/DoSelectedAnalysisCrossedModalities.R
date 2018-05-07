DoSelectedAnalysisCrossedModalities <- function(filename, fom){
  
  cat("***Data file is ", filename, ', and time is ',  date(), '***\n')
  cat("fom = ", fom, '\n')
  cat("alpha = ", alpha, '\n')
  fom <- toupper(fom)
  
  ret <- Read_JAFROC_V2_CrossedModalities( filename )
  Nk <- ret$Nk
  NL <- ret$NL
  LL <- ret$LL
  Wk <- ret$Wk
  ModalityIDs1 <- ret$ModalityIDs1
  ModalityIDs2 <- ret$ModalityIDs2
  
  cat("***First modality levels = ", ModalityIDs1, "***\n")
  
  cat("***Second modality levels = ", ModalityIDs2, "***\n")
  
  I1 <- length(ModalityIDs1)
  I2 <- length(ModalityIDs2)
  J <- length(NL[1,1,,1,1])
  K <- length(NL[1,1,1,,1])
  
  if (fom == toupper("ROC"))  {
    Hratings <- array(dim = c(I1,I2,J,K))
    
    for (i1 in 1:I1) { 
      ret <- ConvertFROC2HRatings(NL[i1,,,,], LL[i1,,,,])
      Hratings[i1,,,] <- ret$Hratings
    }
    truth <- ret$truth
    
    cat("***RESULTS FOR HR INFERRED ROC ANALYSIS ***\n")
    
    JACK <- TRUE
    
    for (AvgIndx in 1:2) { #AvgIndx <- 1 # avg over i1 index  
      
      if (AvgIndx == 1) cat("Results when averaging over all values of first crossed modality\n") else
        cat("Results when averaging over all values of second crossed modality\n")
      
      cov <- OR_Covariances_ROC_CrossedModalities(JACK, Hratings, truth, AvgIndx)
      if (AvgIndx == 1) result_ROC <- OR_Method(I2,J,K,cov) else result_ROC <- OR_Method(I1,J,K,cov)
      
      print(result_ROC)
      FOM <- result_ROC$FOM
      cat("Showing reader-averaged FOMs\n")
      for (i2 in 1:length(FOM[,2])) {
        cat(mean(FOM[i2,]),"\n")  
      }
      cat("\n")
    }
    
  } else {
    cat("***RESULTS FOR wAFROC ANALYSIS ***\n")
    
    JACK <- TRUE
    
    for (AvgIndx in 1:I1) { #AvgIndx <- 1 # avg over i1 index; AvgIndx <- 2 # avg over i2 index 
      
      if (AvgIndx == 1) cat("Results when averaging over all values of first crossed modality\n") else
        cat("Results when averaging over all values of second crossed modality\n")
      
      cov <- OR_Covariances_wAfroc_CrossedModalities(JACK, fom, NL,LL,Nk,Wk,AvgIndx)
      if (AvgIndx == 1) result_JAFROC <- OR_Method(I2,J,K,cov) else result_JAFROC <- OR_Method(I1,J,K,cov)
      
      print(result_JAFROC)
      FOM <- result_JAFROC$FOM
      cat("Showing reader-averaged FOMs\n")
      for (i2 in 1:length(FOM[,2])) {
        cat(mean(FOM[i2,]),"\n")  
      }
      cat("\n")
    }
  }
  
  cat("AllDone!")
  
  
}