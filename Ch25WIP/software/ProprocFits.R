ProprocFits <- function(filename, pathName2) {
  mrmcFile <- paste0(pathName2,"MRMCRuns/", fileName, "_MRMC proproc area pooled.csv")
  if (!file.exists(mrmcFile)) stop("need to run proproc for this dataset")
  proprocRet <- read.csv(mrmcFile)
  c1 <- matrix(data = proprocRet$c, nrow = length(unique(proprocRet$T)), ncol = length(unique(proprocRet$R)), byrow = TRUE)
  da <- matrix(data = proprocRet$d_a, nrow = length(unique(proprocRet$T)), ncol = length(unique(proprocRet$R)), byrow = TRUE)
  return (list(c1 = c1, da = da))
}
