GoodnessOfFitCBM <- function(mu, alpha, zetas, FP, TP){
  while (any(FP < 5)){
    for (b in 1:length(FP)){
      if (FP[b] < 5){
        if (b == 1){
          FP[b + 1] <- FP[b + 1] + FP[b]
          FP <- FP[-b]
        }else{
          FP[b - 1] <- FP[b - 1] + FP[b]
          FP <- FP[-b]
        }
      }
    }
  }
  
  col
}