install.packages(pkgs = "~/Desktop/RJafroc_0.0-1.tar.gz", repo = NULL, type = "source")
rm(list = ls())
library(RJafroc)

pmfLesionDistribution <- rbind(c(1, 1), c(2, 0))
OperatingCharacteristics(mu = 0, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 1, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 3, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 4, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 5, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")


OperatingCharacteristics(mu = 2, lambda = 1, beta = 0, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 1, beta = 2, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 1, beta = 3, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 1, beta = 4, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 1, beta = 5, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")

OperatingCharacteristics(mu = 2, lambda = 0, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 2, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 3, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 4, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
OperatingCharacteristics(mu = 2, lambda = 5, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")

pmfLesionDistribution <- rbind(c(1, 1), c(2, 0))
OperatingCharacteristics(mu = 2, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
pmfLesionDistribution <- rbind(c(1, .75), c(2, .25))
OperatingCharacteristics(mu = 2, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
pmfLesionDistribution <- rbind(c(1, .5), c(2, .5))
OperatingCharacteristics(mu = 2, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
pmfLesionDistribution <- rbind(c(1, .25), c(2, .75))
OperatingCharacteristics(mu = 2, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
pmfLesionDistribution <- rbind(c(1, 0), c(2, 1))
OperatingCharacteristics(mu = 2, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")
pmfLesionDistribution <- rbind(c(1, 0), c(2, 0.5), c(3, 0.5))
OperatingCharacteristics(mu = 2, lambda = 1, beta = 1, pmfLesionDistribution = pmfLesionDistribution, type = "ROC")











for(mu in seq(0,5,1)) {
  OperatingCharacteristics(mu = mu, lambda = 1, beta = .5, pmfLesionDistribution = pmfLesionDistribution)
}


pmfLesionDistribution <- rbind(c(1, 35), c(2, 25), c(3, 10))
OperatingCharacteristics(mu = c(2, 3), lambda = c(1, 1.5), beta = c(0.6, 0.8),
                         pmfLesionDistribution = pmfLesionDistribution)

pmfLesionDistribution <- rbind(c(1, 35), c(2, 25), c(3, 10))
OperatingCharacteristics(mu = c(2, 3), lambda = c(1, 1), beta = c(0.8, 0.8),
                         pmfLesionDistribution = pmfLesionDistribution)

ret  <- OperatingCharacteristics(mu = 2.5, lambda = 1, beta = 0.8,
                                 pmfLesionDistribution = pmfLesionDistribution)

ret  <- OperatingCharacteristics(mu = 2.5, lambda = 1, beta = 0.8,
                                 pmfLesionDistribution = pmfLesionDistribution)



