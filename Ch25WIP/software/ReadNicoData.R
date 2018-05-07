ReadNicoData <- function (RADIOLOGISTS, residents) 
{
  truth <- read.table("jaf_truth.txt", sep = ",")
  CaseID <- truth$V1
  LesionID <- truth$V2
  K1 <- length(LesionID[LesionID == 0])
  K2 <- length(LesionID[LesionID == 1])
  K <- c(K1, K2)
  CaseIDNor <- CaseID[LesionID == 0]
  CaseIDAbn <- CaseID[LesionID == 1]
  
  data <- read.table("findings.txt", sep = "")
  CaseID <- data$V1
  ReaderID <- data$V2
  # increment each rating by one, so it extends from 1 to 101; dpc 12/14/16
  Rating <- data$V13 + 1 #  sum(is.na(Rating)) should give 0
  CL <- as.logical(data$V14)
  
  ReaderIDArray <- unique(ReaderID)
  temp <- match(residents, ReaderIDArray)
  Radiologists <- ReaderIDArray[-temp]
  if (RADIOLOGISTS) ReaderIDArray <- Radiologists 
  J <- length(ReaderIDArray)
  
  zjkt <- array(dim = c(J, max(c(K1,K2)), 2))
  CLArray <- array(FALSE, dim = c(J, K2))
  for (i in 1:length(CaseID)){
    if ( CaseID[i] %in% CaseIDNor ) t <- 1 else t <- 2
    j <- match(ReaderID[i], ReaderIDArray)
    if (is.na(j)) next
    if (t == 1) k <- match(CaseID[i], CaseIDNor) else k <- match(CaseID[i], CaseIDAbn)
    #if (t == 2) cat("case, j,k,t, =  ", CaseID[i], j, k, t, "\n")
    if (t == 1) {      
      if (is.na(zjkt[j,k,t])) zjkt[j,k,t] <- Rating[i] else zjkt[j,k,t] <- max(c(zjkt[j,k,t],Rating[i]))
    } else {      
      if (is.na(zjkt[j,k,t])) {
        zjkt[j,k,t] <- Rating[i]
        CLArray[j,k] <- CL[i]
      } else {
        if (Rating[i] > zjkt[j,k,t]) {
          zjkt[j,k,t] <- Rating[i]
          CLArray[j,k] <- CL[i]
        }
      }      
    }
  }
  UnMarked <- 0
  zjkt[is.na(zjkt)] <- UnMarked  # sum(is.na(zjkt)) should give zero
  zjk1 <- zjkt[,1:K1,1]
  zjk2Cl <- array(dim = c(J, K2))
  temp1 <- (CLArray == TRUE)
  zjk2 <- zjkt[,1:K2,2]
  zjk2Cl[temp1] <- zjk2[temp1]
  
  zjk2Il <- array(dim = c(J, K2))
  temp2 <- (CLArray == FALSE)
  zjk2Il[temp2] <- zjk2[temp2]
  
  # sum(is.na(zjk1)) should give 0
  zjk2Cl[is.na(zjk2Cl)] <- UnMarked
  zjk2Il[is.na(zjk2Il)] <- UnMarked
  
  # make CAD the first reader
  zjk1 <- zjk1[c(J,1:(J-1)),]
  zjk2Cl <- zjk2Cl[c(J,1:(J-1)),]
  zjk2Il <- zjk2Il[c(J,1:(J-1)),]
  
  return (list (
    zjk1 = zjk1,
    zjk2Cl = zjk2Cl,
    zjk2Il = zjk2Il
  ))
}  
