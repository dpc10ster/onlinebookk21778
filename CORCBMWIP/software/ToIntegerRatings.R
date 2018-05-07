# Included functions are:
#
# CalculateZetas <- function (Z1,Z2)
# ToIntegerRatings.R <- function (f,t,DesiredNumBins)
#
#
# categorized integer or floating point z-samples
# to 5 bins for compatibility with CorROCII
# considered using CLabRoc, but extra dimensions of Covariancd matrix
# could make positive definiteness problem worse
CalculateZetas <- function (Z1,Z2)
{
  
  K1 <- length(Z1)
  K2 <- length(Z2)
  
  #   Data <- array(dim=c(K1+K2,2))
  #   Data[1:K1,1] <- Z1;Data[1:K1,2] <- -1
  #   Data[(K1+1):(K1+K2),1] <- Z2;Data[(K1+1):(K1+K2),2] <- +1
  #   ret <-  sort(Data[,1],index.return=TRUE)
  #   Data1 <- Data[ret$ix,]
  #   
  #   fp1 <- Data1[which(Data1[,2]==-1),1]
  #   tp1 <- Data1[which(Data1[,2]==+1),1]  
  fp1 <- Z1
  tp1 <- Z2
  fp <- Z1
  tp <- Z2
  maxCtffs <- 200
  zetas <- array(dim=maxCtffs)
  z1 <- array(dim=maxCtffs)
  z2 <- array(dim=maxCtffs)
  r <- 1
  while(1) {
    if (r > maxCtffs) stop("Need more cutoffs")
    zetas[r] <- min( c( fp1, tp1 ) )
    fp1 <- fp[fp > zetas[r]]
    tp1 <- tp[tp > zetas[r]]
    z1[r] <- length( fp )-length( fp1 )
    z2[r] <- length( tp )-length( tp1 )
    r <- r+1
    fp <- fp1
    tp <- tp1
    if( sum(z1[!is.na(z1)]) == K1 && sum(z2[!is.na(z2)]) == K2) break
  }
  zetas <- zetas[!is.na(zetas)]
  return(zetas)
  
}


# categorized integer or floating point z-samples
# to DesiredNumBins bins for compatibility with CorROCII
# considered using CLabRoc, but extra dimensions of Covariance matrix
# could make positive definiteness problem worse
ToIntegerRatings <- function (f,t,DesiredNumBins)
{
  zetas <- CalculateZetas(f,t)    
    
  fb <- as.vector(table(cut(f,c(-Inf, zetas, Inf)))) 
  tb <- as.vector(table(cut(t,c(-Inf, zetas, Inf)))) 
  if (length(fb) != length(tb)) stop
  RocDataTable <- array(dim = c(2,length(fb)))
  RocDataTable[1,1:length(fb)] <- fb;
  RocDataTable[2,1:length(tb)] <- tb;
    
  minValue <- 0
  while(1) {
    R1 <- length(RocDataTable[1,])
    if (R1 == DesiredNumBins) break    
    if ((length(which(RocDataTable[1,] == minValue)) == 0) && 
      (length(which(RocDataTable[2,] == minValue)) == 0))  minValue <- minValue+1  
    for (t1 in 1:2) {
      R1 <- length(RocDataTable[1,])
      if (R1 == DesiredNumBins) break    
      minValueColumn <- which(RocDataTable[t1,] == minValue)[1]
      if (is.na(minValueColumn)) next
      if (minValueColumn == R1) {
        RocDataTable[1,R1-1]  <- sum(RocDataTable[1,(R1-1):R1]) 
        RocDataTable[2,R1-1]  <- sum(RocDataTable[2,(R1-1):R1])
        RocDataTable  <- RocDataTable[,-minValueColumn]
        zetas <- zetas[-(minValueColumn-1)]        
      } else if (minValueColumn == 1) {
        RocDataTable[1,2]  <- sum(RocDataTable[1,1:2]) 
        RocDataTable[2,2]  <- sum(RocDataTable[2,1:2])
        RocDataTable  <- RocDataTable[,-minValueColumn]
        zetas <- zetas[-1]        
      } else {
        RocDataTable[1,minValueColumn-1]  <- sum(RocDataTable[1,(minValueColumn-1):minValueColumn]) 
        RocDataTable[2,minValueColumn-1]  <- sum(RocDataTable[2,(minValueColumn-1):minValueColumn])
        RocDataTable  <- RocDataTable[,-minValueColumn]
        zetas <- zetas[-(minValueColumn-1)]
      }
    }
    if (R1 == (DesiredNumBins)) break      
  }
  
  fb <- as.vector(table(cut(f,c(-Inf, zetas, Inf)))) 
  tb <- as.vector(table(cut(t,c(-Inf, zetas, Inf)))) 
  if (length(fb) != length(tb)) stop("lengths unequal in ToIntegerRatings")
  
  f <- cut(f,c(-Inf, zetas, Inf),labels=FALSE) 
  t <- cut(t,c(-Inf, zetas, Inf),labels=FALSE) 
  
  return(list(fb=fb, tb=tb, f=f,t=t,zetas=zetas))
}
