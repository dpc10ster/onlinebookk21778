PlotsValues <- function (I, saveSim, CI, bsMean, bsVar, nhCondition) 
{
  # test simulated cov1, cov2 etc. and fom wrt 95% CI's: expect a few failures
  if (I == 1) {
    Failed <- array(NA, dim = c(3, S)) # cov2, var, fom
    compStr <- c("Cov2","Var","fom")
  } else {
    Failed <- array(NA, dim = c(4 + I, S)) # cov1, cov2, cov3, var, fom1, fom2
    compStr <- c("Cov1","Cov2","Cov3","Var",paste0("fom", 1:I))
  }
  # generate histogram plots with 95% confidence intervals 
  if (!nhCondition) {
    for (n in 1:length(compStr)){
      cat(compStr[n], "\n")
      cat("FailureRate:", "\n")
      cat("MeanSIM:", "\n")
      cat("ORIG:", "\n")
      cat("MeanBOOT:", "\n")
      cat("95%CI:", "\n")
      cat("95%CISIM:", "\n")
      cat("StdBOOT:", "\n")
      cat("StdSIM:", "\n")
      cat("Ignore:", "\n")
      cat("Ignore:", "\n")
      cat("Ignore:", "\n")
      cat("Ignore:", "\n")
    }
    
    for (n in 1:length(compStr)){
      temp <- (sapply(saveSim, "[[", n)) # nth element in compStr
      meanTemp <- mean(temp)
      Failed[n, ] <- !((CI[n, 1] < temp) & (CI[n, 2] > temp))
      simuData <- data.frame(simu = temp)
      names(simuData) <- compStr[n]
      histogram <- ggplot(data = simuData, mapping = aes_string(compStr[n])) + geom_histogram(fill = "grey") +
        geom_segment(aes(x = CI[n, 1], y = 0, xend = CI[n, 1], yend = 0), linetype = 2) + 
        geom_segment(aes(x = CI[n, 2], y = 0, xend = CI[n, 2], yend = 0), linetype = 2)
      maxCounts <-suppressMessages(max(ggplot_build(histogram)$data[[1]]$count))
      histogram <- histogram + 
        geom_segment(aes(x = CI[n, 1], y = 0, xend = CI[n, 1], yend = maxCounts + 1), linetype = 2) + 
        geom_segment(aes(x = CI[n, 2], y = 0, xend = CI[n, 2], yend = maxCounts + 1), linetype = 2) + 
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                           axis.text = element_text(size = 15), axis.title = element_text(size = 15)) + 
        scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
      suppressMessages(print(histogram))
      cat(sum(Failed[n, ])/S, "\n")
      cat(meanTemp, "\n")
      cat(ret$orgVarCompFom[n], "\n")
      cat(bsMean[n], "\n")
      cat(CI[n, 2] - CI[n, 1], "\n")
      simuCI <- quantile(temp, c(0.025, 0.975))
      cat(simuCI[2] - simuCI[1], "\n")
      bsSd <- sqrt(bsVar[n])
      cat(bsSd, "\n")
      sdTemp <- sqrt(var(temp))
      cat(sdTemp, "\n")
      meanRatio <- meanTemp / ret$orgVarCompFom[n]
      cat(meanRatio, "\n")
      meanRatio <- meanTemp / bsMean[n]
      cat(meanRatio, "\n")
      ciRatio <- (simuCI[2] - simuCI[1]) / (CI[n, 2] - CI[n, 1])
      cat(ciRatio, "\n")
      sdRatio <- sdTemp / bsSd
      cat(sdRatio, "\n\n")
    }
  }
  if (nhCondition == TRUE){
    cat("Mean (FOM2 - FOM1):", "\n")
    cat("Rej rate:", "\n")
    fom1 <- (sapply(saveSim, "[[", 5))
    fom2 <- (sapply(saveSim, "[[", 6))
    fomDiff <- data.frame(fomDiff = fom2 - fom1)
    #     histogram <- ggplot(data = fomDiff, mapping = aes_string(names(fomDiff))) + geom_histogram(fill = "grey") + 
    #       theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    #                          axis.text = element_text(size = 15), axis.title = element_text(size = 15)) + 
    #       scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    #     suppressMessages(print(histogram))
    cat(mean(fomDiff$fomDiff), "\n")
    pValue <- data.frame(pValue = sapply(saveSim, "[[", 4 + I + 3))
    names(pValue) <- "pValue"
    #     histogram <- ggplot(data = pValue, mapping = aes_string(names(pValue))) + geom_histogram(fill = "grey")
    #     maxCounts <-suppressMessages(max(ggplot_build(histogram)$data[[1]]$count))
    #     histogram <- histogram + 
    #       geom_segment(aes(x = 0.05, y = 0, xend = 0.05, yend = maxCounts + 1), linetype = 2) + 
    #       theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    #                          axis.text = element_text(size = 15), axis.title = element_text(size = 15)) + 
    #       scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    #     suppressMessages(print(histogram))
    cat(sum(pValue <= 0.05)/S, "\n")
  }
  
}