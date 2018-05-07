rm(list = ls()) #mainSsOrhCheck2.R
library(RJafroc)
alpha <- 0.05
effectSize <-0.05
KStar <- 114 # VanDyke Data
# following are proproc generated values from Hillis 2011 for Van Dyke Data
Cov1 <- 0.000351859
Cov2 <- 0.000346505
Cov3 <- 0.000221453
Var <- 0.001393652
msTR <- 0.000623
VarTR <- msTR - Var + Cov1 + max(Cov2-Cov3,0)
VarTR <- max(VarTR,0)
cat("VarTR = ", VarTR, "Var = ", Var, 
    ", Cov1 = ", Cov1, ", Cov2 = ", Cov2, ", Cov3  = ", Cov3, 
    ", effectSize = ", effectSize, "\n")
powTab <- SsPowerTable(alpha = alpha, effectSize = effectSize, desiredPower = 0.8,  
                       method = "ORH", option = "RRRC", Cov1, Cov2, Cov3, VarTR, Var, KStar)
print(powTab) 
# above compare to Table 5 left quarter in Hillis 2011; exact match
# 