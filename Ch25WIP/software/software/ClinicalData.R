
#########################################################################
## need to run proproc and manually update these values for selected data
# dataset  <- "JOHN"
# dataset  <- "FEDERICA_ALL"
# dataset  <- "TONY"
# for (dataset in c("TONY", "FEDERICA_ALL", "JOHN")){
#   seed <- 123
#   if (dataset == "TONY") {
#     #   AUCProp <- array(dim = c(2, 5))
#     #   AUCProp[1, ] <- c(0.801, 0.895, 0.853, 0.858, 0.891)
#     #   AUCProp[2, ] <- c(0.672, 0.754, 0.793, 0.874, 0.736)
#     fitDataA <- ReadDataFile("./Datasets/ALMLCWBZSZ_Finale_20100402.xlsx")
#     AUCCbmA <- CBMAuc(fitDataA, c(1,2), c(1:5), seed)
#   } else if (dataset == "FEDERICA_ALL") {
#     #   AUCProp <- array(dim = c(5, 4))
#     #   AUCProp[1, ] <- c(0.90771537, 0.83820585, 0.82293475, 0.91688438)
#     #   AUCProp[2, ] <- c(0.88291590, 0.86748922, 0.84645401, 0.90737956)
#     #   AUCProp[3, ] <- c(0.85344896, 0.85076174, 0.76871996, 0.88532448)
#     #   AUCProp[4, ] <- c(0.91055126, 0.86005371, 0.81093809, 0.91300701)
#     #   AUCProp[5, ] <- c(0.85926055, 0.80124530, 0.78650757, 0.88331324)
#     fitDataB <- ReadDataFile("./Datasets/FedericaAll.xlsx")  
#     AUCCbmB <- CBMAuc(fitDataB, c(1:5), c(1:4), seed)
#   } else if (dataset == "JOHN") {
#     #   AUCProp[1, ] <- c(0.86191891, 0.94589972, 0.89672137, 0.76956721, 0.72402757, 0.94051911, 0.93007420, 0.88968301, 0.85230899)
#     #   AUCProp[2, ] <- c(0.89193451, 0.98189765, 0.91123782, 0.86893216, 0.86005077, 0.93143807, 0.97049487, 0.94301712, 0.86941160)
#     fitDataC <- ReadDataFile("./Datasets/PET-CT (10AR).xlsx")   
#     AUCCbmC <- CBMAuc(fitDataC, c(1,2), c(1:9), seed)
#   } else stop("undefined dataset")
#   
#   if (dataset == "TONY"){
#     TONY <- FitRsmRocCurve(fitDataA, 1:2, 1:5, AUCProp = AUCCbmA)  
#   } else if (dataset == "FEDERICA_ALL") {
#     FEDE <- FitRsmRocCurve(fitDataB, 1:5, 1:4, AUCProp = AUCCbmB)  
#   } else if (dataset == "JOHN") {
#     JOHN <- FitRsmRocCurve(fitDataC, 1:2, 1:9, AUCProp = AUCCbmC)  
#   }
# }
# 
# mu <- mean(c(as.numeric(TONY$mu$mu), as.numeric(FEDE$mu$mu), as.numeric(JOHN$mu$mu)))
# lambda <- mean(c(as.numeric(TONY$lambda$lambda), as.numeric(FEDE$lambda$lambda), as.numeric(JOHN$lambda$lambda)))
# nu <- mean(c(as.numeric(TONY$nu$nu), as.numeric(FEDE$nu$nu), as.numeric(JOHN$nu$nu)))

# mu <- 3.090848
