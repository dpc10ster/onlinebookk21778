# mainNpvPpv.R
rm(list = ls())
# disease prevalence in US screening mammography
prevalence <- 0.005
FPF <- 0 # typical operating point
TPF <- 0 # do:
specificity <- 1-FPF
sensitivity <- TPF
NPV <- (1-prevalence)*(specificity)/
  ((1-prevalence)*(specificity) + 
     prevalence*(1-sensitivity))
PPV <- prevalence*sensitivity/
  (prevalence*sensitivity + 
     (1-prevalence)*(1-specificity))
cat("NPV = ", NPV, "\nPPV = ", PPV, "\n")
accuracy <-(1-prevalence)*
  (specificity)+(prevalence)*(sensitivity)
cat("accuracy = ", accuracy, "\n")
