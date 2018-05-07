# fp <- c(rep(5, 99), rep(1, 1))
# tp <- rep(1, 100)
# degData <- ToRJafrocDataset(fp, tp, "ROC")
# EmpiricalOpCharac(degData, 1, 1)
# SaveDataFile(degData, "0.99, 0.lrc", "MRMC")

degData <- ReadDataFile("0.01, 0.lrc", "MRMC")
FitCBM(degData, 1, 1)
