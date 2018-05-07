# MySyncFile.R
DataFile <- "JS" # "FZ_FROC", "FZ_REAL", "JT", "FR", "VD", "TS", "EK", "JS" 
simSelect <- "CZ_MRMC" # use "CZ_MRMC" or "RM" to select the simulator
if (simSelect == "RM") RMVarStr <- "HH" ############ "HL", "LL", "HH", "LH" ############# and "EST"
Validation <- TRUE # FALSE for Calibration mode, TRUE for Validation mode
overWrite <- FALSE #set to TRUE to overwrite existing files, FALSE otherwise
aucEst <- "RocFit" # DO NOT CHANGE THIS DPC!! 

if (Validation == TRUE) { 
  nhCondition <- FALSE # first run with FALSE flag (to get validation results), then TRUE flag (to get NH rejection rate)
  # set to FALSE for debuggingg; otherwise this MUST BE TRUE
  resampleChi <- TRUE # set to FALSE for debugging/demonstration; otherwise this MUST BE TRUE
  # set to FALSE for debugging/demonstration; otherwise this MUST BE TRUE
  if (nhCondition == TRUE) {
    S <- 2000  
    aucEst <- "Trap"
  } else S <- 200 
  nBoots <- 200
}

selI <- NULL # set NULL to select all modalities
selJ <- NULL
if (DataFile == "JS") selI <- c(3,4) # to include only conventinal and dual energy modalities from JS data
##if (DataFile == "JS") selI <- c(1,3) # to exclude dual energy modalities from JS data
if (DataFile == "JT") selJ <- c(2,3,6,7,8) # to exclude reader 4 from analysis of JT data
if (DataFile == "VD") selJ <- c(1,2,3,5) # to exclude reader 4 from analysis of VD data
simuI <- NULL # NULL for same number as in orginal dataset
simuJ <- NULL # NULL for same number as in orginal dataset
if (DataFile == "JS") simuJ <- 7 # to see if inc. J improves results; was 5 before this
simuK1 <- 50  # NULL for same number as in orginal dataset
simuK2 <- 50  # NULL for same number as in orginal dataset
DesiredNumBins <- 4 # continuous or > than 4 bin integer ratings will be re-binned into 4 ROC bins

if (!is.null((selI)) || !is.null(selJ)){
  if (is.null((selI))){
    Istring <- "ALL"
    Jstring <- paste(selJ, collapse = "_")
  }else if (is.null(selJ)){
    Istring <- paste(selI, collapse = "_")
    Jstring <- "ALL"
  }else{
    Istring <- paste(selI, collapse = "_")
    Jstring <- paste(selJ, collapse = "_")
  }
}else{
  Istring <- "ALL"
  Jstring <- "ALL"
}

UNINITIALIZED <<- -2000 # global
MAX_MARKS <<- 10 # global
a_min <<- 0.001 # clamps on range of allowed values
a_max <<- 6
b_min <<- 0.001 
b_max <<- 6

seed <- 1;set.seed(seed)
#seed <- NULL;set.seed(seed)
alpha <- 0.05

