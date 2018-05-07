rm(list = ls()) # mainDataSummary.R # freeze lines
library(RJafroc)

# included datasets
fileNames <-  c("TONY", "VD", "FR", 
                "FED", "JT", "MAG", 
                "OPT", "PEN", "NICO",
                "RUS", "DOB1", "DOB2", 
                "DOB3", "FZR")

descriptions <- c(
  "Digital breast tomosynthesis vs. mammography",
  "Cine vs. SE MRI for aortic dissection",
  "Digital vs. analog pediatric chest",
  "Image processing in mammography: FROC",
  "Nodule detection in an thorax CT phantom",
  "Tomosynthesis Vs. Radiography Pulmonary Nodules",
  "Calcification detection in digital mammography",
  "Image compression in mammography",
  "Standalone CAD vs. radiologists mammography",
  "Lesion detection in digital mammography",
  "Tomosynthesis, Dual-Energy & Conventional Chest", 
  "Tomosynthesis, Dual-Energy & Conventional Chest",
  "Tomosynthesis, Dual-Energy & Conventional Chest",
  "Image processing in mammography: ROC"
)
cat("Dataset Name, # modalities, # readers, # non-diseased cases, # diseased cases, general description of data\n")
for (f in 1:length(fileNames)){
  fileName <- fileNames[f]
  retFileName <- paste0("allResults", fileName) 
  sysAnalFileName <- system.file("ANALYZED/RSM6", retFileName, package = "RJafroc")
  if (file.exists(sysAnalFileName)){
    load(sysAnalFileName)
    I <- allResults[[1]]$I
    J <- allResults[[1]]$J
    K1 <- allResults[[1]]$K1
    K2 <- allResults[[1]]$K2
    cat(f, fileNames[f], I, J, K1, K2, descriptions[f], "\n")
  }
}
