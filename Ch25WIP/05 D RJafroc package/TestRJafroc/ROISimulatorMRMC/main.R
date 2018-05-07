rm( list = ls() )
require("xlsx")
source("ROISimulator.R")
source("SaveJafrocFile.R")

seed.value <- 3
set.seed(seed.value)

I <- 2
J <- 5
mu <- 2
K = c(50, 40)
Q <- 4
Deltamu_i <- c(0, 0.4)#, 0.2 ) # ! the number of elements must match the number of treatments

# origin of next 6 numbers are from Roe and Metz validation paper, which needs validation!!
Var_R <- 0.2
Var_mR <- 0.005

# next four numbers have to add up to one
Var_C <- 0.7
Var_mC <- 0.05
Var_RC <- 0.2
Var_e <- 0.05

if ( abs(Var_C + Var_mC + Var_RC + Var_e - 1) > 1e-4 ) stop("check variance components")

rho_C <- 0.1
rho_RC <- 0.1
rho_mC <- 0.9
rho_e <- 0.9

ret <- ROISimulator (I, J, K, Q, mu, Deltamu_i, Var_R, Var_C, Var_mR, 
                                 Var_mC, Var_RC, Var_e, rho_C, rho_RC, rho_mC, rho_e) 
Rijkttlss <- ret$Rijkttlss
isDiseased <- ret$isDiseased

SaveJafrocFile ( K, Rijkttlss , isDiseased, "RoiData.xlsx") 




