VarPopSampling <- function ( K, mu, sigma, zetas, P)
{
  AUC <- array(dim = P) 
  for (p in 1:P)
  {
    z <- GenerateCaseSamples(K, mu, sigma, zetas)
    z1 <- z$z1
    z2 <- z$z2
    AUC[p] <- Wilcoxon(z1, z2)
    # RocCountsTableSimPop <- SimulateRocCountsTable(K, mu, sigma, zetas)
    # RocCountsTableSimPop[is.na(RocCountsTableSimPop )] <- 0#replace NAs with zeroes
    # AUC[p] <- (WilcoxonCountsTable(RocCountsTableSimPop))$Az
  }
  
  return (var(AUC))
}
# if (FOM == "Az") { # AUC for observed data
#   RocCountsTable  <- array(dim = c(2,length(RocCountsTable[1,])))
#   bsTable <- table(z1)
#   RocCountsTable[1, as.numeric(names(bsTable))] <- bsTable#convert array to frequency table
#   bsTable <- table(z2)
#   RocCountsTable[2, as.numeric(names(bsTable))] <- bsTable #do:
#   RocCountsTable[is.na(RocCountsTable )] <- 0 #replace NAs with zeroes
#   ret <- RocfitR(RocCountsTable)                 
# } else {
