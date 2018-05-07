lesionDistribution <- function(frocData)
{  
  lesionNum <- frocData$lesionNum
  nLesDistr <- table(lesionNum)
  if (length(nLesDistr) == 1) {
    nLesDistr <- c(lesionNum[1], 1)
    dim(nLesDistr) <- c(1, 2)
  }else{
    nLesDistr <- t(rbind(as.numeric(unlist(attr(nLesDistr, "dimnames"))), as.vector(nLesDistr)))
  }
  return(nLesDistr)
}