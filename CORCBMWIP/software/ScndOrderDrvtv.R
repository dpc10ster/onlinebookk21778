ScndOrder <- function(parameters, nBins, Kdd1, Kdd2){
  allParameters <- names(parameters)
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(parameters[allParameters[which(zetaXPos == 1)]])
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(parameters[allParameters[which(zetaYPos == 1)]])
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  muX <- parameters$muX
  muY <- parameters$muY
  alpha <- parameters$alpha
  rho1 <- parameters$rhoNor
  rho2 <- parameters$rhoAbn
  
  scndOrderDrvtv <- array(0, dim = c(length(parameters), length(parameters)))
  scndOrderDrvtv <- scndOrderDrvtv + ScndMuX(parameters, nBins, Kdd1, Kdd2)
  return(scndOrderDrvtv)
}

ScndMuX <- function(parameters, nBins, Kdd1, Kdd2){
  allParameters <- names(parameters)
  nParam <- which(regexpr("muX",  allParameters) == 1)
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(parameters[allParameters[which(zetaXPos == 1)]])
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(parameters[allParameters[which(zetaYPos == 1)]])
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  muX <- parameters$muX
  muY <- parameters$muY
  alpha <- parameters$alpha
  rho1 <- parameters$rhoNor
  rho2 <- parameters$rhoAbn
  
  scndOrderDrvtv <- array(0, dim = c(length(parameters), length(parameters)))
  for (s in 1:nBins){
    for (r in 1:nBins){
      den <- F2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
      - F2(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2(zetaX[r], zetaY[s], muX, muY, alpha, rho2)
      for (par in 1:length(parameters)){
        if (par %in% c(1:5)){
          scndOrderDrvtv[1, par] <- scndOrderDrvtv[1, par] +
            switch (par,
                    "1" = Kdd2 * ((F2dmuX(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuX(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                   - F2dmuX(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuX(zetaX[r], zetaY[s], muX, muY, alpha, rho2))^2 / den),
                    "2" = Kdd2 * ((F2dmuX(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuX(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                   - F2dmuX(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuX(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                    (F2dmuY(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuY(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                     - F2dmuY(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuY(zetaX[r], zetaY[s], muX, muY, alpha, rho2))
                                  / den),
                    "3" = Kdd2 * ((F2dmuX(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuX(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                   - F2dmuX(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuX(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                    (F2dalpha(zetaX[r + 1], zetaY[s + 1], muX, muY, rho2) - F2dalpha(zetaX[r], zetaY[s + 1], muX, muY, rho2) 
                                     - F2dalpha(zetaX[r + 1], zetaY[s], muX, muY, rho2) + F2dalpha(zetaX[r], zetaY[s], muX, muY, rho2))
                                  / den),
                    "4" = 0,
                    "5" = Kdd2 * ((F2dmuX(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuX(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                   - F2dmuX(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuX(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                    (F2drho2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2drho2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                     - F2drho2(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2drho2(zetaX[r], zetaY[s], muX, muY, alpha, rho2))
                                  / den)
            )
        }
        if (par == 5 + r){
          scndOrderDrvtv[1, par] <- scndOrderDrvtv[1, par] + Kdd2 * ((F2dmuX(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuX(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                                                      - F2dmuX(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuX(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                       (F2dx(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dx(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2))
                                                                     / den)
        }
        if (par == 5 + length(zetaX) - 2 + s){
          scndOrderDrvtv[1, par] <- scndOrderDrvtv[1, par] + Kdd2 * ((F2dmuX(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuX(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                                                      - F2dmuX(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuX(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                       (F2dy(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dy(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2))
                                                                     / den)
        }
        if (any(is.infinite(scndOrderDrvtv))){
          dummyStop <- 1
        }
      }
    }
  }
  scndOrderDrvtv[, 1] <- scndOrderDrvtv[1, ]
  return(-scndOrderDrvtv)
}


ScndMuY <- function(parameters, nBins, Kdd1, Kdd2){
  allParameters <- names(parameters)
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(parameters[allParameters[which(zetaXPos == 1)]])
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(parameters[allParameters[which(zetaYPos == 1)]])
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  muX <- parameters$muX
  muY <- parameters$muY
  alpha <- parameters$alpha
  rho1 <- parameters$rhoNor
  rho2 <- parameters$rhoAbn
  
  scndOrderDrvtv <- array(0, dim = c(length(parameters), length(parameters)))
  for (s in 1:nBins){
    for (r in 1:nBins){
      den <- F2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
      - F2(zetaX[r], zetaY[s], muX, muY, alpha, rho2) + F2(zetaX[r], zetaY[s], muX, muY, alpha, rho2)
      for (par in 1:length(parameters)){
        if (par %in% c(1:5)){
          scndOrderDrvtv[2, par] <- scndOrderDrvtv[2, par] +
            switch (par,
                    "1" = 0,
                    "2" = Kdd2 * ((F2dmuY(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuY(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                   - F2dmuY(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuY(zetaX[r], zetaY[s], muX, muY, alpha, rho2))^2
                                  / den),
                    "3" = Kdd2 * ((F2dmuY(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuY(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                   - F2dmuY(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuY(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                    (F2dalpha(zetaX[r + 1], zetaY[s + 1], muX, muY, rho2) - F2dalpha(zetaX[r], zetaY[s + 1], muX, muY, rho2) 
                                     - F2dalpha(zetaX[r + 1], zetaY[s], muX, muY, rho2) + F2dalpha(zetaX[r], zetaY[s], muX, muY, rho2))
                                  / den),
                    "4" = 0,
                    "5" = Kdd2 * ((F2dmuY(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuY(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                   - F2dmuY(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuY(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                    (F2drho2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2drho2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                     - F2drho2(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2drho2(zetaX[r], zetaY[s], muX, muY, alpha, rho2))
                                  / den)
            )
        }
        if (par == 5 + r){
          scndOrderDrvtv[2, par] <- scndOrderDrvtv[2, par] + Kdd2 * ((F2dmuY(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuY(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                                                      - F2dmuY(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuY(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                       (F2dx(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dx(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2))
                                                                     / den)
        }
        if (par == 5 + length(zetaX) - 2 + s){
          scndOrderDrvtv[2, par] <- scndOrderDrvtv[2, par] + Kdd2 * ((F2dmuY(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dmuY(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                                                      - F2dmuY(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dmuY(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                       (F2dy(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dy(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2))
                                                                     / den)
        }
        if (any(is.infinite(scndOrderDrvtv))){
          dummyStop <- 1
        }
      }
    }
  }
  scndOrderDrvtv[, 2] <- scndOrderDrvtv[2, ]
  return(-scndOrderDrvtv)
}


ScndAlpha <- function(parameters, nBins, Kdd1, Kdd2){
  allParameters <- names(parameters)
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(parameters[allParameters[which(zetaXPos == 1)]])
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(parameters[allParameters[which(zetaYPos == 1)]])
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  muX <- parameters$muX
  muY <- parameters$muY
  alpha <- parameters$alpha
  rho1 <- parameters$rhoNor
  rho2 <- parameters$rhoAbn
  
  scndOrderDrvtv <- array(0, dim = c(length(parameters), length(parameters)))
  for (s in 1:nBins){
    for (r in 1:nBins){
      den <- F2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
      - F2(zetaX[r], zetaY[s], muX, muY, alpha, rho2) + F2(zetaX[r], zetaY[s], muX, muY, alpha, rho2)
      for (par in 1:length(parameters)){
        if (par %in% c(1:5)){
          scndOrderDrvtv[3, par] <- scndOrderDrvtv[3, par] +
            switch (par,
                    "1" = 0,
                    "2" = 0,
                    "3" = Kdd2 * ((F2dalpha(zetaX[r + 1], zetaY[s + 1], muX, muY, rho2) - F2dalpha(zetaX[r], zetaY[s + 1], muX, muY, rho2) 
                                   - F2dalpha(zetaX[r + 1], zetaY[s], muX, muY, rho2) + F2dalpha(zetaX[r], zetaY[s], muX, muY, rho2))^2
                                  / den),
                    "4" = 0,
                    "5" = Kdd2 * ((F2dalpha(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dalpha(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                   - F2dalpha(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dalpha(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                    (F2drho2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2drho2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                     - F2drho2(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2drho2(zetaX[r], zetaY[s], muX, muY, alpha, rho2))
                                  / den)
            )
        }
        if (par == 5 + r){
          scndOrderDrvtv[3, par] <- scndOrderDrvtv[3, par] + Kdd2 * ((F2dalpha(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dalpha(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                                                      - F2dalpha(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dalpha(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                       (F2dx(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dx(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2))
                                                                     / den)
        }
        if (par == 5 + length(zetaX) - 2 + s){
          scndOrderDrvtv[3, par] <- scndOrderDrvtv[3, par] + Kdd2 * ((F2dalpha(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dalpha(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                                                      - F2dalpha(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dalpha(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                       (F2dy(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dy(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2))
                                                                     / den)
        }
        if (any(is.infinite(scndOrderDrvtv))){
          dummyStop <- 1
        }
      }
    }
  }
  scndOrderDrvtv[, 3] <- scndOrderDrvtv[3, ]
  return(-scndOrderDrvtv)
}


ScndRho1 <- function(parameters, nBins, Kdd1, Kdd2){
  allParameters <- names(parameters)
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(parameters[allParameters[which(zetaXPos == 1)]])
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(parameters[allParameters[which(zetaYPos == 1)]])
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  muX <- parameters$muX
  muY <- parameters$muY
  alpha <- parameters$alpha
  rho1 <- parameters$rhoNor
  rho2 <- parameters$rhoAbn
  
  scndOrderDrvtv <- array(0, dim = c(length(parameters), length(parameters)))
  for (s in 1:nBins){
    for (r in 1:nBins){
      den <- F1(zetaX[r + 1], zetaY[s + 1], rho1) - F1(zetaX[r], zetaY[s + 1], rho1) 
      - F1(zetaX[r], zetaY[s], rho1) + F1(zetaX[r], zetaY[s], rho1)
      for (par in 1:length(parameters)){
        if (par %in% c(1:5)){
          scndOrderDrvtv[4, par] <- scndOrderDrvtv[4, par] +
            switch (par,
                    "1" = 0,
                    "2" = 0,
                    "3" = 0,
                    "4" = Kdd2 * ((F1drho1(zetaX[r + 1], zetaY[s + 1], rho1) - F1drho1(zetaX[r], zetaY[s + 1], rho1) 
                                   - F1drho1(zetaX[r + 1], zetaY[s], rho1) + F1drho1(zetaX[r], zetaY[s], rho1))^2
                                  / den),
                    "5" = 0
            )
        }
        if (par == 5 + r){
          scndOrderDrvtv[4, par] <- scndOrderDrvtv[4, par] + Kdd1 * ((F1drho1(zetaX[r + 1], zetaY[s + 1], rho1) - F1drho1(zetaX[r], zetaY[s + 1], rho1) 
                                                                      - F1drho1(zetaX[r + 1], zetaY[s], rho1) + F1drho1(zetaX[r], zetaY[s], rho1)) * 
                                                                       (F1dx(zetaX[r + 1], zetaY[s + 1], rho1) - F1dx(zetaX[r + 1], zetaY[s], rho1))
                                                                     / den)
        }
        if (par == 5 + length(zetaX) - 2 + s){
          scndOrderDrvtv[4, par] <- scndOrderDrvtv[4, par] + Kdd1 * ((F1drho1(zetaX[r + 1], zetaY[s + 1], rho1) - F1drho1(zetaX[r], zetaY[s + 1], rho1) 
                                                                      - F1drho1(zetaX[r + 1], zetaY[s], rho1) + F1drho1(zetaX[r], zetaY[s], rho1)) * 
                                                                       (F1dy(zetaX[r + 1], zetaY[s + 1], rho1) - F1dy(zetaX[r], zetaY[s + 1], rho1))
                                                                     / den)
        }
        if (any(is.infinite(scndOrderDrvtv))){
          dummyStop <- 1
        }
      }
    }
  }
  scndOrderDrvtv[, 4] <- scndOrderDrvtv[4, ]
  return(-scndOrderDrvtv)
}


ScndRho2 <- function(parameters, nBins, Kdd1, Kdd2){
  allParameters <- names(parameters)
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(parameters[allParameters[which(zetaXPos == 1)]])
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(parameters[allParameters[which(zetaYPos == 1)]])
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  muX <- parameters$muX
  muY <- parameters$muY
  alpha <- parameters$alpha
  rho1 <- parameters$rhoNor
  rho2 <- parameters$rhoAbn
  
  scndOrderDrvtv <- array(0, dim = c(length(parameters), length(parameters)))
  for (s in 1:nBins){
    for (r in 1:nBins){
      den <- F2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
      - F2(zetaX[r], zetaY[s], muX, muY, alpha, rho2) + F2(zetaX[r], zetaY[s], muX, muY, alpha, rho2)
      for (par in 1:length(parameters)){
        if (par %in% c(1:5)){
          scndOrderDrvtv[5, par] <- scndOrderDrvtv[5, par] +
            switch (par,
                    "1" = 0,
                    "2" = 0,
                    "3" = 0,
                    "4" = 0,
                    "5" = Kdd2 * ((F2drho2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2drho2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                   - F2drho2(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2drho2(zetaX[r], zetaY[s], muX, muY, alpha, rho2))^2
                                  / den)
            )
        }
        if (par == 5 + r){
          scndOrderDrvtv[5, par] <- scndOrderDrvtv[5, par] + Kdd2 * ((F2drho2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2drho2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                                                      - F2drho2(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2drho2(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                       (F2dx(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dx(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2))
                                                                     / den)
        }
        if (par == 5 + length(zetaX) - 2 + s){
          scndOrderDrvtv[5, par] <- scndOrderDrvtv[5, par] + Kdd2 * ((F2drho2(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2drho2(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                                                      - F2drho2(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2drho2(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                       (F2dy(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dy(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2))
                                                                     / den)
        }
        if (any(is.infinite(scndOrderDrvtv))){
          dummyStop <- 1
        }
      }
    }
  }
  scndOrderDrvtv[, 5] <- scndOrderDrvtv[5, ]
  return(-scndOrderDrvtv)
}


ScndZetaX <- function(parameters, nBins, Kdd1, Kdd2){
  allParameters <- names(parameters)
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(parameters[allParameters[which(zetaXPos == 1)]])
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(parameters[allParameters[which(zetaYPos == 1)]])
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  muX <- parameters$muX
  muY <- parameters$muY
  alpha <- parameters$alpha
  rho1 <- parameters$rhoNor
  rho2 <- parameters$rhoAbn
  
  scndOrderDrvtv <- array(0, dim = c(length(parameters), length(parameters)))
  for (r in 1:(nBins - 1)){
    for (s in 1:nBins){
      den11 <- F1(zetaX[r + 1], zetaY[s + 1], rho1) - F1(zetaX[r], zetaY[s + 1], rho1) 
      - F1(zetaX[r + 1], zetaY[s], rho1) + F1(zetaX[r], zetaY[s], rho1)
      den12 <- F1(zetaX[r + 2], zetaY[s + 1], rho1) - F1(zetaX[r + 1], zetaY[s + 1], rho1) 
      - F1(zetaX[r + 2], zetaY[s], rho1) + F1(zetaX[r], zetaY[s], rho1)
      den2 <- F2(zetaX[rp + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2(zetaX[rp], zetaY[s + 1], muX, muY, alpha, rho2) 
      - F2(zetaX[rp], zetaY[s], muX, muY, alpha, rho2) + F2(zetaX[rp], zetaY[s], muX, muY, alpha, rho2)
      for (par in 1:length(parameters)){
        if (par %in% c(1:5)){
          scndOrderDrvtv[5, par] <- scndOrderDrvtv[5, par] +
            switch (par,
                    "1" = 0,
                    "2" = 0,
                    "3" = 0,
                    "4" = 0,
                    "5" = 0
            )
        }
        if (par == 5 + r - 1){
          if (par > 5){
            scndOrderDrvtv[5 + r, par] <- scndOrderDrvtv[5 + r, par] + Kdd2 * ((F2dx(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dx(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                                 (-F2dx(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) + F2dx(zetaX[r], zetaY[s], muX, muY, alpha, rho2))
                                                                               / den2)
            +  Kdd1 * ((F1dx(zetaX[r + 1], zetaY[s + 1], rho1) - F1dx(zetaX[r + 1], zetaY[s], rho1)) * 
                         (-F1dx(zetaX[r], zetaY[s + 1], rho1) + F1dx(zetaX[r], zetaY[s], rho1))
                       / den1)
          }
        }
        if (par == 5 + r){
          scndOrderDrvtv[5 + r, par] <- scndOrderDrvtv[5 + r, par] + Kdd2 * ((F2dx(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dx(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2))^2 / den2)
          + Kdd1 * ((F1dx(zetaX[r + 1], zetaY[s + 1], rho1) - F1dx(zetaX[r + 1], zetaY[s], rho1))^2 / den1)
        }
        
        if (par == 5 + r + 1){
          if (par < 5 + nBins){
            scndOrderDrvtv[5 + r, par] <- scndOrderDrvtv[5 + r, par] + Kdd1 * ((F1dx(zetaX[r + 1], zetaY[s + 1], rho1) - F1dx(zetaX[r + 1], zetaY[s], rho1)) * 
                        (-F1dx(zetaX[r + 2], zetaY[s + 1], rho1) + F1dx(zetaX[r + 2], zetaY[s], rho1))
                      / den1)
            + Kdd2 * ((F2dx(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dx(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2)) * 
                        (-F2dx(zetaX[r + 2], zetaY[s + 1], muX, muY, alpha, rho2) + F2dx(zetaX[r + 2], zetaY[s], muX, muY, alpha, rho2))
                      / den2)
          }
        }
        
        if (par == 5 + nBins - 1 + s - 1){
          if (par > 5 + nBins - 1){
            scndOrderDrvtv[5 + r, par] <- scndOrderDrvtv[5 + r, par] + Kdd2 * ((F2dx(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dx(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                                 (-F2dy(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) + F2dy(zetaX[r], zetaY[s], muX, muY, alpha, rho2))
                                                                               / den2)
            + Kdd1 * ((F1dx(zetaX[r + 1], zetaY[s + 1], rho1) - F1dx(zetaX[r + 1], zetaY[s], rho1)) * 
                        (-F1dy(zetaX[r + 1], zetaY[s + 1], rho1) + F1dy(zetaX[r], zetaY[s], rho1))
                      / den1)
          }
          
        }
        if (par == 5 + nBins - 1 + s){
          scndOrderDrvtv[5 + r, par] <- scndOrderDrvtv[5 + r, par] + Kdd2 * ((F2dy(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dy(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2) 
                                                                               - F2dy(zetaX[r + 1], zetaY[s], muX, muY, alpha, rho2) + F2dy(zetaX[r], zetaY[s], muX, muY, alpha, rho2)) * 
                                                                                (F2dy(zetaX[r + 1], zetaY[s + 1], muX, muY, alpha, rho2) - F2dy(zetaX[r], zetaY[s + 1], muX, muY, alpha, rho2))
                                                                              / den)
        }
        if (any(is.infinite(scndOrderDrvtv))){
          dummyStop <- 1
        }
      }
    }
  }
  
  scndOrderDrvtv[, 6:(5 + nBins - 1)] <- scndOrderDrvtv[6:(5 + nBins - 1), ]
  return(-scndOrderDrvtv)
}