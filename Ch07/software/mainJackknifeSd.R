rm( list = ls()) # MainJackknifeSd.R

source("Transforms.R")
source("LL.R")
source("RocfitR.R")
source("RocOperatingPoints.R")
source("FixRocCountsTable.R")
source("WilcoxonCountsTable.R")

FOM  <- "Az"
#FOM  <- "Wilcoxon"
cat("FOM = ", FOM, "\n")
RocCountsTable = array(dim = c(2,5))
RocCountsTable[1,]  <- c(30,19,8,2,1)
RocCountsTable[2,]  <- c(5,6,5,12,22)

K <- c(sum(RocCountsTable[1,]), sum(RocCountsTable[2,])) # this is the K vector 

if (FOM == "Az") { # AUC for observed data
  Az <- RocfitR(RocCountsTable)     
  if (Az$Az == -1) stop("RocfitR did not converge on original data")
} else {
  Az <- WilcoxonCountsTable(RocCountsTable) 
}
Az = Az$Az # AUC for observed data

z1  <- rep(1:length(RocCountsTable[1,]),RocCountsTable[1,])#convert frequency table to array
z2  <- rep(1:length(RocCountsTable[1,]),RocCountsTable[2,])#do:

AUC  <- array(dim = sum(K)); Y <- array(dim = sum(K))
z_jk  <- array(dim = sum(K))
for ( k in 1 : sum(K)){
  RocCountsTable_jk  <- array(dim = c(2,length(RocCountsTable[1,])))
  if ( k <= K[ 1 ]){
    z1_jk <- z1[ -k ]
    z2_jk <- z2
  }else{
    z1_jk <- z1
    z2_jk <- z2[ -(k - K[ 1 ]) ] 
  }  
  RocCountsTable_jk[1,1:length(table(z1_jk))]  <- table(z1_jk)#convert array to frequency table
  RocCountsTable_jk[2,1:length(table(z2_jk))]  <- table(z2_jk)#do:
  RocCountsTable_jk[is.na(RocCountsTable_jk)] <- 0#replace NAs with zeroes
  if (FOM == "Az") {
    temp <- RocfitR(RocCountsTable_jk)                # AUC for observed data 
  } else {
    temp <- WilcoxonCountsTable(RocCountsTable_jk) 
  }      
  AUC[k]  <- temp$Az
  Y[k] <- sum(K)*Az - (sum(K)-1)*AUC[k]
  if (AUC[k] == -1) stop("RocfitR did not converge in jackknife loop") 
}

Var <- var(AUC) * ( sum(K) - 1)^2 / sum(K) #Efron and Stein's paper
stdAUC  <- sqrt(Var)
cat("OrigAUC = ", Az, "\njackknifeMeanAuc = ", mean(AUC), "\nstdAUC = ", stdAUC, "\n")
