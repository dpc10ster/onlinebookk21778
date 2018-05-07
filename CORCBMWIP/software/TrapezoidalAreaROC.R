# extended it to handle most general case
# nice thing about this code is that the normal and abnormal cases can be mixed; the truth vector tells
# which is which
TrapezoidalAreaROC <- function (RocRatings, truth)
{
  MultipleModality <- FALSE
  NestedModality <- FALSE
  data <- dim(RocRatings)
  if (length(data) == 4) {
    NestedModality <- TRUE
    I1 <- data[1]
    I2 <- data[2]
    J <- data[3]
  } else if (length(data) == 3) {
    MultipleModality <- TRUE
    I <- data[1]
    J <- data[2]
  } else if (length(data) == 2) {
    J <- data[1]
  } else stop("Error in TrapezoidalAreaROC function")
  
  K2 <- sum(truth)
  K1 <- length(truth) - K2
  
  K1Indx <- (1:length(truth))[!as.logical(truth)] #index of incorrect cases; or normal cases
  K2Indx <- (1:length(truth))[as.logical(truth)] #index of correct cases; or abnormal cases
  
  if (NestedModality && !MultipleModality){
    AUC <- array( dim = c(I1,I2,J))
    for (i1 in 1 : I1){
      for (i2 in 1 : I2){
        for (j in 1 : J) {
          S <- 0
          for (k in K1Indx) {    
            S <- S + sum(RocRatings[i1,i2,j,k] < RocRatings[i1,i2,j,K2Indx]) + 0.5 * sum(RocRatings[i1,i2,j,k] == RocRatings[i1,i2,j,K2Indx])
          }
          AUC[i1,i2,j] <- S/(K2 * K1)
        }
      }
    }
  }
  
  if (!NestedModality && MultipleModality){
    AUC <- array( dim = c(I, J))
    for (i in 1 : I){
      for (j in 1 : J) {
        S <- 0
        for (k in K1Indx) {    
          S <- S + sum(RocRatings[i, j, k] < RocRatings[i, j, K2Indx]) + 0.5 * sum(RocRatings[i, j, k] == RocRatings[i, j, K2Indx])
        }
        AUC[ i, j ] <- S/(K2 * K1)
      }
    }    
  }
  
  if (!NestedModality && !MultipleModality){
    AUC <- array( dim = J)
    for (j in 1 : J) {
      S <- 0
      for (k in K1Indx) {    
        S <- S + sum(RocRatings[j,k] < RocRatings[j,K2Indx]) + 0.5 * sum(RocRatings[j, k] == RocRatings[j, K2Indx])
      }
      AUC[ j ] <- S/(K2 * K1)   
    }
  }
  
  return (AUC)
}