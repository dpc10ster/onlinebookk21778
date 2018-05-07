#' Operating characteristics and highest rating pdfs for specified values of search model parameters for up to two treatments
#' 
#' Visualize the predicted ROC, AFROC, FROC and pdfs (probability density functions of highest rating, 
#' for non-diseased and diseased cases) for up to 2 specified sets of search model parameters 
#' corresponding to up to two treatments.
#' 
#' @param mu The mean(s) of the Gasussian distribution(s) for the signal sites (continuous ratings of lesions that 
#' are found by the search mechanism) for up to two treatments.
#' @param lambda The Poisson distribution parameter(s), which describes the random number of noise sites 
#' (suspicious regions that do not correspond to actual lesions) per case for up to two treatments.
#' @param nu 1 - exp(nu*mu) is the success probability of the binomial distribution(s) describing the random number of signal sites 
#' (suspicious regions that correspond to actual lesions) per case, for up to two treatments.
#' @param pmfLesionDistribution The probability mass function (pmf) of the lesion distribution. 
#' A 2D array, where the first dimension contains the integers 1 through the maximum number 
#' of lesions per case over the whole dataset, and the second dimension contains the fraction of diseased cases 
#' with the number of lesions specified in the first dimension. The sum over the
#' second dimension must equal one.
#' @param type The type of operating characteristic desired: can be "\code{ROC}", "\code{AFROC}", 
#' "\code{FROC}" or "\code{PDFs}" or "\code{ALL}".
#' @param legendPosition The positioning of the legend: "\code{right}", "\code{left}", "\code{top}" or "\code{bottom}".
#' 
#' @return It returns a list of 6 elements containing four \pkg{ggplot2} objects 
#' (ROCPlot, FROCPlot, AFROCPlot and PDFPlot) and two area measures (each of which can have up to two elements), 
#' the area under the search model predicted ROC curves
#' in up to two treatments, and the area under the search model predicted AFROC curves
#' in up to two treatments. \bold{The area under the FROC is a poor measure of performance}.
#' 
#' 
#' @import ggplot2
#' 
#' 
#' 
#' @references
#' Chakraborty, D. P. (2006). ROC curves predicted by a model of visual search. Physics in Medicine and Biology, 51(14), 3463-82.
#' 
#' Chakraborty, D. P., & Yoon, H.-J. (2008). Operating characteristics predicted by models for diagnostic tasks involving lesion localization. Medical Physics, 35(2), 435.
#' 
#' Chakraborty, D. P. (2006). A search model and figure of merit for observer data acquired according to the free-response paradigm. Physics in Medicine and Biology, 51(14), 3449-62.
#' 
#' @examples
#' pmfLesionDistribution <- rbind(c(1, 0.2), c(2, 0.4), c(3, 0.1), c(2, 0.3))
#' lesionWeights <- rbind(c(1.0, -Inf, -Inf), c(0.4,  0.6, -Inf), c(0.2,  0.3,  0.5), c(0.7,  0.3, -Inf))
#' OperatingCharacteristics(mu = c(2, 3), lambda = c(1, 1.5), nu = c(0.6, 0.8),
#'    pmfLesionDistribution = pmfLesionDistribution, lesionWeights = lesionWeights, legendPosition = "bottom")
#' ## For mu = 2, lambda = 1, nu = 0.6, in one treatment and mu = 3, lambda = 1.5, nu = 0.8, 
#' ## in the other treatment. 70% of the diseased cases have a single lesion, 20% have have two lesions and 10% have 3 lesions.
#' 
OperatingCharacteristics <- function(mu, lambda, nu, pmfLesionDistribution, lesionWeights, type = "ALL", legendPosition = "bottom"){
  if (!all(c(length(mu) == length(lambda), length(mu) == length(nu))))
    stop("Parameters have different lengths.")
  
  if (nrow(pmfLesionDistribution) != nrow(lesionWeights)){
    stop("lesionWeights and pmfLesionDistribution must have same number of rows.")
  }
  
  if ((length(pmfLesionDistribution) == 1) && is.integer(pmfLesionDistribution)){
    pmfLesionDistribution <- c(pmfLesionDistribution, 1)
  }else{
    if (ncol(pmfLesionDistribution) != 2){
      stop("pmfLesionDistribution must have two columns")
    }
  }
  
  for (r in 1:nrow(lesionWeights)){
    rowWeight <- lesionWeights[r, ]
    nWeight <- sum(rowWeight != -Inf)
    if (sum(rowWeight[rowWeight != -Inf]) != 1){
      errMsg <- sprintf("Line %d of lesion weights matrix should be summed up to 1.", r)
      stop(errMsg)
    }
    if (nWeight != pmfLesionDistribution[r , 1]){
      errMsg <- sprintf("The number of elements in the line %d of lesion weights matrix is different the number of lesion in lesion distribution matrix.", r)
      stop(errMsg)
    }
  }
  
  plotStep <- 0.005
  zeta <- seq(from = -20, to = 20, by = plotStep)
  
  ROCPlot <- NA
  LROCPlot <- NA
  FROCPlot <- NA
  AFROCPlot <- NA
  PDFPlot <- NA
  
  ROCPoints <- data.frame(FPF = NULL, TPF = NULL, Treatment = NULL)
  ROCDashes <- data.frame(FPF = NULL, TPF = NULL, Treatment = NULL)
  LROCPoints <- data.frame(FPF = NULL, TPF = NULL, Treatment = NULL)
  LROCDashes <- data.frame(FPF = NULL, TPF = NULL, Treatment = NULL)
  FROCPoints <- data.frame(NLF = NULL, LLF = NULL, Treatment = NULL)
  FROCDashes <- data.frame(NLF = NULL, LLF= NULL, Treatment = NULL)
  AFROCPoints <- data.frame(FPF = NULL, LLF= NULL, Treatment = NULL)
  AFROCDashes <- data.frame(FPF = NULL, LLF= NULL, Treatment = NULL)
  wAFROCPoints <- data.frame(FPF = NULL, wLLF= NULL, Treatment = NULL)
  wAFROCDashes <- data.frame(FPF = NULL, wLLF= NULL, Treatment = NULL)
  abnPDFPoints <- data.frame(pdf = NULL, highestZSample = NULL, Treatment = NULL)
  norPDFPoints <- data.frame(pdf = NULL, highestZSample = NULL, Treatment = NULL)
  aucROC <- rep(NA, length(mu))
  aucLROC <- aucROC
  aucAFROC <- aucROC
  aucWAFROC <- aucROC
  lambdaP <- lambda
  nuP <- nu
  for (i in 1:length(mu)){
    if (nu[i] < 0 ) stop("nu must be non-negative")
    
    lambdaP[i] <- lambda[i] / mu[i]
    if (abs(nu[i] * mu[i]) <= 1e-6 ) nuP[i] <- 1e-6 else nuP[i] <- (1-exp(-nu[i] * mu[i]))
    FPF <- sapply(zeta, xROC, lambdaP = lambdaP[i])
    TPF <- sapply(zeta, yROC, mu = mu[i], lambdaP = lambdaP[i], nuP = nuP[i], pmfLesionDistribution = pmfLesionDistribution)
    PCL <- sapply(zeta, yLROC, mu = mu[i], lambdaP = lambdaP[i], nuP = nuP[i])
    ROCPoints <- rbind(ROCPoints, data.frame(FPF = FPF, TPF = TPF, Treatment = as.character(i)))
    ROCDashes <- rbind(ROCDashes, data.frame(FPF = c(FPF[1], 1), TPF = c(TPF[1], 1), Treatment = as.character(i)))
    
    LROCPoints <- rbind(LROCPoints, data.frame(FPF = FPF, PCL = PCL, Treatment = as.character(i)))
    LROCDashes <- rbind(LROCDashes, data.frame(FPF = c(FPF[1], 1), PCL = c(PCL[1], 1), Treatment = as.character(i)))
    
    NLF <- sapply(zeta, xFROC, lambdaP = lambdaP[i])
    LLF <- sapply(zeta, yFROC, mu = mu[i], nuP = nuP[i])
    FROCPoints <- rbind(FROCPoints, data.frame(NLF = NLF, LLF = LLF, Treatment = as.character(i)))
    
    AFROCPoints <- rbind(AFROCPoints, data.frame(FPF = FPF, LLF = LLF, Treatment = as.character(i)))
    AFROCDashes <- rbind(AFROCDashes, data.frame(FPF = c(FPF[1], 1), LLF = c(LLF[1], 1), Treatment = as.character(i)))
    
    wLLF <- sapply(zeta, yWAFROC, mu = mu[i], nuP = nuP[i], pmfLesionDistribution = pmfLesionDistribution, lesionWeights = lesionWeights)
    wAFROCPoints <- rbind(wAFROCPoints, data.frame(FPF = FPF, wLLF = wLLF, Treatment = as.character(i)))
    wAFROCDashes <- rbind(wAFROCDashes, data.frame(FPF = c(FPF[1], 1), wLLF = c(wLLF[1], 1), Treatment = as.character(i)))
    
    deltaFPF <- FPF[1:(length(FPF) - 1)] - FPF[2:length(FPF)]  
    if(type == "ALL" || type == "PDFs"){     
      pdfNor <- deltaFPF / plotStep
      norPDFPoints <- rbind(norPDFPoints, 
                            data.frame(pdf = pdfNor[pdfNor > 1e-6], highestZSample = zeta[-1][pdfNor > 1e-6], 
                                       Treatment = as.character(i), class = "Nor"))
      deltaTPF <- TPF[1:(length(TPF) - 1)] - TPF[2:length(TPF)]
      pdfAbn <- deltaTPF / plotStep
      abnPDFPoints <- rbind(abnPDFPoints, 
                            data.frame(pdf = pdfAbn[pdfAbn > 1e-6], highestZSample = zeta[-1][pdfAbn > 1e-6], 
                                       Treatment = as.character(i), class = "Abn"))
    }
    
    #aucROC[i] <- sum((TPF[1:(length(TPF) - 1)] + TPF[2:length(TPF)]) * deltaFPF / 2) + (TPF[1] + 1) * (1 - FPF[1]) / 2
    maxFPF <- xROC(-20, lambdaP[i])
    maxTPF <- yROC(-20, mu[i], lambdaP[i], nuP[i], pmfLesionDistribution)
    AUC <- integrate(intROC, 0, maxFPF, mu = mu[i], lambdaP = lambdaP[i], nuP = nuP[i], pmfLesionDistribution = pmfLesionDistribution)$value
    aucROC[i] <- AUC + (1 + maxTPF) * (1 - maxFPF) / 2
    
    #aucAFROC[i] <- sum((LLF[1:(length(LLF) - 1)] + LLF[2:length(LLF)]) * deltaFPF / 2) + (LLF[1] + 1) * (1 - FPF[1]) / 2
    maxLLF <- yFROC(-20, mu[i], nuP[i])
    AUC <- integrate(intAFROC, 0, maxFPF, mu = mu[i], lambdaP = lambdaP[i], nuP = nuP[i])$value
    aucAFROC[i] <- AUC + (1 + maxLLF) * (1 - maxFPF) / 2
    
    maxWLLF <- yWAFROC(-20, mu[i], nuP[i], pmfLesionDistribution, lesionWeights)
    AUC <- integrate(intWAFROC, 0, maxFPF, mu = mu[i], lambdaP = lambdaP[i], nuP = nuP[i], pmfLesionDistribution, lesionWeights)$value
    aucWAFROC[i] <- AUC + (1 + maxWLLF) * (1 - maxFPF) / 2
    
    maxPCL <- yLROC(-20, mu[i], lambdaP[i], nuP[i])
    AUC <- integrate(intLROC, 0, maxFPF, mu = mu[i], lambdaP = lambdaP[i], nuP = nuP[i])$value
    aucLROC[i] <- AUC + (1 + maxPCL) * (1 - maxFPF) / 2 
  }
  
  if(type == "ALL" || type == "ROC"){
    ROCPlot <- with(ROCPoints, {
      ggplot(data = ROCPoints) + geom_line(aes(x = FPF, y = TPF, color = Treatment)) + 
        geom_line(data = ROCDashes, aes(x = FPF, y = TPF, color = Treatment), linetype = 2) +       
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + theme(legend.position = legendPosition)  
    })
  }
  
  if(type == "ALL" || type == "FROC"){
    FROCPlot <- with(FROCPoints, {
      ggplot(data = FROCPoints) + geom_line(aes(x = NLF, y = LLF, color = Treatment))  +       
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + theme(legend.position = legendPosition)    
    })
  }
  
  if(type == "ALL" || type == "AFROC"){
    AFROCPlot <- with(AFROCPoints, {
      ggplot(data = AFROCPoints) + geom_line(aes(x = FPF, y = LLF , color = Treatment)) + 
        geom_line(data = AFROCDashes, aes(x = FPF, y = LLF, color = Treatment), linetype = 2) +       
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + theme(legend.position = legendPosition)   
    })
  }
  
  if(type == "ALL" || type == "wAFROC"){
    wAFROCPlot <- with(wAFROCPoints, {
      ggplot(data = wAFROCPoints) + geom_line(aes(x = FPF, y = wLLF , color = Treatment)) + 
        geom_line(data = wAFROCDashes, aes(x = FPF, y = wLLF, color = Treatment), linetype = 2) +       
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + theme(legend.position = legendPosition)   
    })
  }
  
  if(type == "ALL" || type == "LROC"){
    LROCPlot <- with(LROCPoints, {
      ggplot(data = LROCPoints) + geom_line(aes(x = FPF, y = PCL , color = Treatment)) + 
        geom_line(data = LROCDashes, aes(x = FPF, y = PCL, color = Treatment), linetype = 2) +       
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + theme(legend.position = legendPosition)   
    })
  }
  
  if(type == "ALL" || type == "PDFs"){
    if (legendPosition == "top" || legendPosition == "bottom"){
      legendArrange = "horizontal"
    }else{
      legendArrange = "vertical"
    }
    PDFPoints <- rbind(norPDFPoints, abnPDFPoints)
    PDFPlot <- with(PDFPoints, {
      ggplot(data = PDFPoints, aes(x = highestZSample, y = pdf, color = Treatment, linetype = class)) + 
        geom_line()  + theme(legend.position = legendPosition, legend.box = legendArrange) + labs(x = "Highest Z Sample") 
    })
  }
  
  
  return(list(
    ROCPlot = ROCPlot,
    AFROCPlot = AFROCPlot,
    wAFROCPlot = wAFROCPlot,
    LROCPlot = LROCPlot,
    FROCPlot = FROCPlot,
    PDFPlot = PDFPlot,
    aucROC = aucROC,
    aucAFROC = aucAFROC,
    aucWAFROC = aucWAFROC,
    aucLROC = aucLROC
  ))
}

erf <- function(x){
  return (2 * pnorm(sqrt(2) * x) - 1)
}

xROC <- function(zeta, lambdaP){
  # returns FPF, the abscissa of ROC curve
  return( 1 - exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2))))
}

yROC <- function(zeta, mu, lambdaP, nuP, pmfLesionDistribution){
  # returns TPF, the ordinate of ROC curve
  fl <- pmfLesionDistribution[, 2] / sum(pmfLesionDistribution[, 2])
  TPF <- 0
  for (i in 1:nrow(pmfLesionDistribution)){
    TPF <- TPF + fl[i] * (1 - (1 - nuP/2 + nuP/2  *erf( (zeta - mu) / sqrt(2) ))^pmfLesionDistribution[i, 1] * exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2))))
  }
  return (TPF)
}

intROC <- function(FPF, mu, lambdaP, nuP, pmfLesionDistribution){
  # returns TPF, the ordinate of ROC curve; takes FPF as the variable. 
  # AUC is calculated by integrating this function in terms of FPF
  tmp <- 1 / lambdaP * log(1 - FPF) + 1
  tmp[tmp < 0] <- pnorm(-20)
  zeta <- qnorm(tmp)
  TPF <- yROC(zeta, mu, lambdaP, nuP, pmfLesionDistribution)
  return (TPF)
}

xFROC <- function(zeta, lambdaP){
  # returns NLF, the abscissa of FROC curve
  NLF <- lambdaP * (1 - pnorm(zeta))
  return(NLF)
}

yFROC <- function(zeta, mu, nuP){
  # returns LLF, the ordinate of FROC, AFROC curve
  LLF <- nuP * (1 - pnorm(zeta - mu))
  return(LLF)
}

intAFROC <- function(FPF, mu, lambdaP, nuP){
  # returns LLF, the ordinate of AFROC curve; takes FPF as the variable. 
  # AUC is calculated by integrating this function in terms of FPF
  tmp <- 1 / lambdaP * log(1 - FPF) + 1
  tmp[tmp < 0] <- pnorm(-20)
  zeta <- qnorm(tmp)
  LLF <- yFROC(zeta, mu, nuP)
  return(LLF)
}

yWAFROC <- function(zeta, mu, nuP, pmfLesionDistribution, lesionWeights){
  # returns wLLFL, the ordinate of wAFROC curve
  fl <- pmfLesionDistribution[, 2] / sum(pmfLesionDistribution[, 2])
  wLLF <- 0
  for (L in 1:nrow(pmfLesionDistribution)){
    nLesion <- pmfLesionDistribution[L, 1] 
    # nLesion is the first element in the row L of pmfLesionDistribution, 
    # which is the number of lesions for this lesion weights distributions condition
    wLLFTmp <- 0
    for (l in 1:nLesion){
      # l is the number of sucesses with number of lesions nLesion
      wLLFTmp <- wLLFTmp + sum(lesionWeights[L, 1:l]) * dbinom(l, nLesion, nuP) * (1 - pnorm(zeta - mu))
      
    }
    wLLF <- wLLF + fl[L] * wLLFTmp
  }
  return(wLLF)
}

intWAFROC <- function(FPF, mu, lambdaP, nuP, pmfLesionDistribution, lesionWeights){
  # returns wLLF, the ordinate of AFROC curve; takes FPF as the variable. 
  # AUC is calculated by integrating this function in terms of FPF
  tmp <- 1 / lambdaP * log(1 - FPF) + 1
  tmp[tmp < 0] <- pnorm(-20)
  zeta <- qnorm(tmp)
  wLLF <- sapply(zeta, yWAFROC, mu = mu, nuP = nuP, pmfLesionDistribution, lesionWeights)
  return(wLLF)
}

yLROC <- function(zeta, mu, lambdaP, nuP){
  # returns PCL, the ordinate of LROC curve
  total <- integrate(intPCL, -Inf, Inf, mu = mu, lambdaP = lambdaP)$value
  return(nuP * (total - integrate(intPCL, -Inf, zeta, mu = mu, lambdaP = lambdaP)$value))
}

intPCL <- function(y, mu, lambdaP){
  # the integral to calculate PCL
  return( dnorm(mu - y) * exp(-lambdaP/2 + lambdaP/2 * erf(y/sqrt(2))) )
}

intLROC <- function(FPF, mu, lambdaP, nuP){
  # returns PCL, the ordinate of LROC curve; takes FPF as the variable. 
  # AUC is calculated by integrating this function in terms of FPF
  tmp <- 1 / lambdaP * log(1 - FPF) + 1
  tmp[tmp < 0] <- pnorm(-20)
  zeta <- qnorm(tmp)
  PCL <- sapply(zeta, yLROC, mu = mu, lambdaP = lambdaP, nuP = nuP)
  return(PCL)
}