rm( list = ls()) # mainBootstrapSd.R

source("Transforms.R")
source("LL.R")
source("RocfitR.R")
source("RocOperatingPoints.R")
source("FixRocCountsTable.R")
source("WilcoxonCountsTable.R")

FOM  <- "Az"
#FOM  <- "Wilcoxon"
B <- 200;seed <- 1;set.seed(seed)
cat("FOM = ", FOM, ", seed = ", seed, ", B = ", B, "\n")
RocCountsTable = array(dim = c(2,5))
RocCountsTable[1,]  <- c(30,19,8,2,1)
RocCountsTable[2,]  <- c(5,6,5,12,22)

K <- c(sum(RocCountsTable[1,]), sum(RocCountsTable[2,])) # this is the K vector 

if (FOM == "Az") { # AUC for observed data
  ret <- RocfitR(RocCountsTable)
  if (ret$Az == -1) stop("RocfitR did not converge on original data")
} else {
  ret <- WilcoxonCountsTable(RocCountsTable) 
}
Az  <- ret$Az

# ready to bootstrap; first put the counts data into a linear form
z1  <- rep(1:length(RocCountsTable[1,]),RocCountsTable[1,])#convert counts table to array
z2  <- rep(1:length(RocCountsTable[2,]),RocCountsTable[2,])#do:
AUC <- array(dim = B)#to save the bs AUC values
for ( b in 1 : B){
  while (1) {
    RocCountsTable_bs  <- array(dim = c(2,length(RocCountsTable[1,])))
    k1_b <- ceiling( runif( K[ 1 ] ) * K[ 1 ] ) # bs indices for non-diseased    
    k2_b <- ceiling( runif( K[ 2 ] ) * K[ 2 ] ) # bs indices for diseased
    bsTable <- table(z1[k1_b])
    RocCountsTable_bs[1, as.numeric(names(bsTable))] <- bsTable#convert array to frequency table
    bsTable <- table(z2[k2_b])
    RocCountsTable_bs[2, as.numeric(names(bsTable))] <- bsTable #do:
    RocCountsTable_bs[is.na(RocCountsTable_bs )] <- 0 #replace NAs with zeroes
    if (FOM == "Az") {
      temp <- RocfitR(RocCountsTable_bs)                # AUC for observed data 
    } else {
      temp <- WilcoxonCountsTable(RocCountsTable_bs) 
    }    
    AUC[b]  <- temp$Az
    if (AUC[b] != -1) break # a return of -1 means Az did not converge
  }
}
Var <- var(AUC)
stdAUC  <- sqrt(Var)
cat("OrigAUC = ", Az, ", meanAUCbs = ", mean(AUC), ", stdAUC = ", stdAUC, "\n")
