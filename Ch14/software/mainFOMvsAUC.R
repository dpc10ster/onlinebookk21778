rm(list = ls()) # MainFOMvsAUC.R; demonstration of equivalences between FOMs and AUCs
library(RJafroc)
library(caTools)

seed <- 1;set.seed(seed)
Lmax <- 2
# following generates 1 to Lmax lesions per dis. case
K1 <- 7;K2 <- 8;Lk2 <- ceiling(runif(K2, 0, Lmax)) 
cat("K1 = ", K1, ", K2 = ", K2, "\n")
mu <- 1.5;lambda <- 0.8;nu <- 0.8 ;zeta1 <- -1
frocDataRaw <- SimulateFrocDataset(
  mu = mu, 
  lambda = lambda, 
  nu = nu, 
  I = 1, 
  J = 1,
  K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1)

# compare afrocAUC vs. afrocFOM
# compare afrocAUC vs. afrocFOM
afrocPlot <- PlotEmpiricalOperatingCharacteristics(
  frocDataRaw,1,1,opChType = "AFROC")
afrocFOM <- signif(as.numeric(UtilFigureOfMerit(frocDataRaw, FOM = "AFROC")), digits = 8)
# trapz(x, y) function, below, returns the trapezoid integral defined by x and y, 
# which can be used to calculate the trapezoidal area in our example. Use signif() function 
# to round the results to 8 significant digits
afrocAUC <- signif(trapz(afrocPlot$Points$genAbscissa, afrocPlot$Points$genOrdinate), digits = 8) 
print(afrocAUC == afrocFOM)
cat("afrocFOM, afrocAUC = ", afrocAUC, "\n")

# compare wafrocAUC vs. wafrocFOM
# compare wafrocAUC vs. wafrocFOM
wafrocPlot <- PlotEmpiricalOperatingCharacteristics(frocDataRaw,1,1,opChType = "wAFROC")
wafrocFOM <- signif(as.numeric(UtilFigureOfMerit(frocDataRaw, FOM = "wAFROC")), digits = 8)
wafrocAUC <- signif(trapz(wafrocPlot$Points$genAbscissa, wafrocPlot$Points$genOrdinate), digits = 8)
print(wafrocAUC == wafrocFOM) # shows TRUE if equality is satisfied
cat("wafrocFOM, wafrocAUC = ", wafrocAUC, "\n")

# compare afroc1AUC vs. afroc1FOM
# compare afroc1AUC vs. afroc1FOM
afroc1Plot <- PlotEmpiricalOperatingCharacteristics(frocDataRaw,1,1,opChType = "AFROC1")
afroc1FOM <- signif(as.numeric(UtilFigureOfMerit(frocDataRaw, FOM = "AFROC1")), digits = 8)
afroc1AUC <- signif(trapz(afroc1Plot$Points$genAbscissa, afroc1Plot$Points$genOrdinate), digits = 8)
print(afroc1AUC == afroc1FOM) # shows TRUE if equality is satisfied
cat("afroc1FOM, afrocAUC1 = ", afroc1AUC, "\n")

# compare wafroc1AUC vs. wafroc1FOM
# compare wafroc1AUC vs. wafroc1FOM
wafroc1Plot <- PlotEmpiricalOperatingCharacteristics(frocDataRaw,1,1,opChType = "wAFROC1")
wafroc1FOM <- signif(as.numeric(UtilFigureOfMerit(frocDataRaw, FOM = "wAFROC1")), digits = 8)
wafroc1AUC <- signif(trapz(wafroc1Plot$Points$genAbscissa, wafroc1Plot$Points$genOrdinate), digits = 8)
print(wafroc1AUC == wafroc1FOM) # shows TRUE if equality is satisfied
cat("wafroc1FOM, wafroc1AUC = ", wafroc1AUC, "\n")

# compare frocAUC vs. frocFOM
# compare frocAUC vs. frocFOM
frocPlot <- PlotEmpiricalOperatingCharacteristics(frocDataRaw,1,1,opChType = "FROC")
frocFOM <- signif(as.numeric(UtilFigureOfMerit(frocDataRaw, FOM = "FROC")), digits = 8)
frocAUC <- signif(trapz(frocPlot$Points$genAbscissa, frocPlot$Points$genOrdinate), digits = 8)
print(frocAUC == frocFOM) # shows TRUE if equality is satisfied
cat("frocFOM, frocAUC = ", frocAUC, "\n")
