# ReadTheData.R
ReadTheData <- function(DataFile,selI,selJ)
{
  if (DataFile == "JT") {
    # JT FROC data; convert to ROC
    orgData <- ReadDataFile(fileName = "./datasets/JT.xlsx")
    orgData <- FROC2HrROC(orgData)
  } else if (DataFile == "FZ_FROC") {
    # FZ FROC data; convert to ROC
    orgData <- ReadDataFile(fileName = "./datasets/FZ_FROC.xlsx")
    orgData <- FROC2HrROC(orgData)
  } else if (DataFile == "FZ_REAL") {
    # FZ ROC data! for 1st modality in FROC file
    orgData <- ReadDataFile(fileName = "./datasets/FZ_REAL.xlsx")
    orgData <- FROC2HrROC(orgData)
  } else if (DataFile == "TS") {
    # Tony Svahn Data
    orgData <- ReadDataFile(fileName = "./datasets/TS.xlsx")
    orgData <- FROC2HrROC(orgData)
  } else if (DataFile == "FR") {
    orgData <- ReadDataFile(fileName = "./datasets/FR.lrc", format = "MRMC")
  } else if (DataFile == "JS") {
    orgData <- ReadDataFile(fileName = "./datasets/JS.xlsx")
    orgData <- FROC2HrROC(orgData)
  } else if (DataFile == "VD") {
    orgData <- ReadDataFile(fileName = "./datasets/VD.lrc", format = "MRMC")
  } else if (DataFile == "EK") {
    # krupinski bit depth data
    orgData <- ReadDataFile(fileName = "./datasets/EK.xlsx", format = "MRMC")
  } else {     # this is where you can add additional datasets
               # } else if (DataFile == "EK") {
               # orgData <- ReadDataFile(fileName = "./datasets/EK.xlsx", format = "MRMC")
               # orgData <- FROC2HrROC(orgData) # convert to ROC if necessary
    stop("Incorrect dataset specified") 
  }
  
  # extract the readers and modalities for analysis
  if (!is.null((selI)) || !is.null(selJ)){
    if (is.null((selI))){
      orgData <- ExtractDataset(orgData, rdrs = selJ)
    }else if (is.null(selJ)){
      orgData <- ExtractDataset(orgData, trts = selI)
    }else{
      orgData <- ExtractDataset(orgData, trts = selI, rdrs = selJ)
    }
  }
  
  return(orgData)
}