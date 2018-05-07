rm( list = ls()) # mainCalSimulator.R

source("Transforms.R")
source("LL.R")
source("RocfitR.R")
source("RocOperatingPoints.R")
source("FixRocCountsTable.R")
source("SimulateRocCountsTable.R")
source("WilcoxonCountsTable.R")
options(digits = 4)
FOM  <- "Az"
#FOM  <- "Wilcoxon"

seed <- 2;set.seed(seed);P <- 2000#number of pop. samples
cat("seed = ", seed, ", FOM = ", FOM, ", P = ", P, "\n")

RocCountsTable = array(dim = c(2,5))
RocCountsTable[1,]  <- c(30,19,8,2,1)
RocCountsTable[2,]  <- c(5,6,5,12,22)
K <- c(sum(RocCountsTable[1,]), sum(RocCountsTable[2,])) # this is the K vector 

# to build the model we have to do a parametric fit first
ret <- RocfitR(RocCountsTable) # AUC for observed data
AUC_org  <- ret$Az;a <- ret$a;b  <- ret$b;zeta <- ret$zeta;mu  <- a/b; sigma <- 1/b;
zeta <- zeta/b # need to also scale zetas

if (FOM == "Az") { # AUC for observed data
  AUC_org <- RocfitR(RocCountsTable)                 
} else {
  AUC_org <- WilcoxonCountsTable(RocCountsTable) 
}
AUC_org  <- AUC_org$Az;
cat("Calibrated simulator values: a, b, zetas:\n", ret$a, ret$b, ret$zeta, "\n")
AUC <- array(dim = P)#to save the pop sample AUC values
a <- array(dim = P);b <- array(dim = P)
for ( p in 1 : P){
  while (1) {
    RocCountsTableSimPop <- SimulateRocCountsTable(K, mu, sigma, zeta)
    RocCountsTableSimPop[is.na(RocCountsTableSimPop )] <- 0#replace NAs with zeroes
    if (FOM == "Az") {
      temp <- RocfitR(RocCountsTableSimPop)                # AUC for observed data 
      if (temp[1] != -1) {# a return of -1 means RocFitR did not converge
        AUC[p]  <- temp$Az;a[p] <- temp$a;b[p] <- temp$b
        break 
      }
    } else {
      AUC[p] <- (WilcoxonCountsTable(RocCountsTableSimPop))$Az
      break 
    }    
  }
}
Var <- var(AUC)
stdAUC  <- sqrt(Var)

cat("seed = ", seed, 
    "\nOrigAUC = ", AUC_org, 
    "\nmeanAUC = ", mean(AUC), 
    "\nstdAUC = ", stdAUC, "\n")
