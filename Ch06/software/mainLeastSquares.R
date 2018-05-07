rm(list = ls()) # mainLeastSquares.R # freeze line numbers
library(ggplot2)
# ML estimates of a and b 
# (from Eng JAVA program)
# a <- 1.3204; 
# b <- 0.6075 
# # these are not used in program; 
# just there for comparison

# these are from Table 6.11, last two rows
FPF <- c(0.017, 0.050, 0.183, 0.5)  
TPF <- c(0.440, 0.680, 0.780, 0.900)

# apply the PHI_INV function
phiInvFPF <- qnorm(FPF)
phiInvTPF <- qnorm(TPF)

fit <- lm(phiInvTPF~phiInvFPF)
print(fit)
pointsData <- data.frame(
  phiInvFPF = phiInvFPF, 
  phiInvTPF = phiInvTPF)
pointsPlot <- ggplot(
  data = pointsData, 
  mapping = 
    aes(x = phiInvFPF, 
        y = phiInvTPF)) + 
  geom_point(size = 5) + 
  theme(
    axis.title.y = element_text(size = 25,face="bold"),
    axis.title.x = element_text(size = 30,face="bold")) +
  geom_abline(
    slope = fit$coefficients[2], 
    intercept = fit$coefficients[1], size = 2)
print(pointsPlot)