WilcoxonDataTable <- function (RocDataTable)
{
  
  zk1  <- rep(1:length(RocDataTable[1,]),RocDataTable[1,])#convert frequency table to array
  zk2  <- rep(1:length(RocDataTable[2,]),RocDataTable[2,])#do:
  
  K1 = length(zk1)
  K2 = length(zk2)
  
  W <- 0
  for (k1 in 1:K1) {
    W <- W + sum(zk1[k1] < zk2)
    W <- W + 0.5 * sum(zk1[k1] == zk2)
  }
  W <- W/K1/K2
  
  return (list(
    Az = W)
  )
  
}
