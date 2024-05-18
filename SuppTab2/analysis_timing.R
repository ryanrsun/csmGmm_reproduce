# do pleiotropy and replication studies
library(data.table)
library(dplyr)
devtools::install.packages("ryanrsun/csmGmm")
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

# which study to do
args <- commandArgs(trailingOnly=TRUE)
num <- as.numeric(args[1])
aID <- as.numeric(args[2])

# open data
setwd("/rsrch3/home/biostatistics/rsun3/summaryStats")
cleanUKB <- fread("bmi_with_overall.txt") 

# output 
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppTab2/output"

# convergence of EM
oldEps <- 0.01
newEps <- 10^(-5)

if (aID == 1) {
  # LC, CAD 
  testDat <- cleanUKB %>% select(Zlc, Zcad, p_LC, p_CAD)
} else if (aID == 2) {
  # LC, BMI
  testDat <- cleanUKB %>% select(Zlc, Zbmi, p_LC, pBMI)
} else if (aID == 3) {
  # CAD, BMI
  testDat <- cleanUKB %>% select(Zcad, Zbmi, p_CAD, pBMI)
} else if (aID == 4) {
  # Overall LC, CAD
  testDat <- cleanUKB %>% select(Zoverall, Zcad, pOverall, p_CAD)
} else if (aID == 5) {
  # Overall LC, BMI
  testDat <- cleanUKB %>% select(Zoverall, Zbmi, pOverall, pBMI)
} else if (aID == 6) {
  # Replication - sqc LC, UKB LC
  setwd("/rsrch3/home/biostatistics/rsun3/summaryStats")
  cleanUKB <- fread("replication_with_lcoverall.txt")
  testDat <- cleanUKB %>% select(Zlc_ilcco, Zlcukb, p_LC_ilcco, pLCukb)
} else if (aID == 7) {
  # Replication - overall LC, UKB LC
  setwd("/rsrch3/home/biostatistics/rsun3/summaryStats")
  cleanUKB <- fread("replication_with_lcoverall.txt")
  testDat <- cleanUKB %>% select(Zoverall, Zlcukb, pOverall, pLCukb)
} else if (aID == 8) {
  # Replication - CAD, UKB CAD
  setwd("/rsrch3/home/biostatistics/rsun3/summaryStats")
  cleanUKB <- fread("cad_for_replication.txt")
  testDat <- cleanUKB %>% select(Zcad_cardio, Zcadukb, p_CAD_cardio, pCADukb)
} else if (aID == 9) {
  # Replication - BMI, UKB BMI
  setwd("/rsrch3/home/biostatistics/rsun3/summaryStats")
  cleanUKB <- fread("bmi_for_replication.txt")
  testDat <- cleanUKB %>% select(Zgiant, ZBMIukb, pGIANT, pBMIukb)
} else if (aID == 10) {
  # Correlation - LC UKB, CAD UKB
  testDat <- cleanUKB %>% select(Zlcukb, Zcadukb, pLC, pCAD)
} else if (aID == 11) {
  # Correlation - LC UKB, bmi
  testDat <- cleanUKB %>% select(Zlcukb, Zbmi, pLC, pBMI)
} else if (aID == 12) {
  # Correlation - CAD UKB, bmi
  testDat <- cleanUKB %>% select(Zcadukb, Zbmi, pCAD, pBMI)
} else if (aID == 13) {
  # three way - sqc, cad, bmi
  testDat <- cleanUKB %>% select(Zlc, Zcad, Zbmi, p_LC, p_CAD, pBMI)
} else if (aID == 14) {
  # three way - overall, cad, bmi
  testDat <- cleanUKB %>% select(Zoverall, Zcad, Zbmi, pOverall, p_CAD, pBMI)
}
testDat <- testDat %>% as.matrix(.)

# adjustment to not get p-values of 0 needed for DACT and HDMT and our own methods
for (col_it in 1:(ncol(testDat)/2)) {
  tooBig <- which(testDat[, col_it] > 8.1)
  tooSmall <- which(testDat[, col_it] < -8.1)
  if (length(tooBig) > 0) {
    testDat[tooBig, col_it] <- 8.1
  }
  if (length(tooSmall) > 0) {
    testDat[tooSmall, col_it] <- -8.1
  }
}
for (col_it in ((ncol(testDat)/2)+1):ncol(testDat)) {
  tooBig <- which(testDat[, col_it] == 1)
  tooSmall <- which(testDat[, col_it] == 0)
  minVal <- min(testDat[which(testDat[, col_it] > 0), col_it]) 
  if (length(tooBig) > 0) {                
    testDat[tooBig, col_it] <- 0.999        
  }        
  if (length(tooSmall) > 0) {                
    testDat[tooSmall, col_it] <- minVal          
  }
}


# start timing
startTime <- proc.time()[[3]]

# 2D cases pleiotropy
if (aID <= 5) {
  initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
  initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                     matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
  newRes <- symm_fit_ind_EM(testStats = testDat[, 1:2], sameDirAlt = FALSE, initMuList = initMuList, initPiList = initPiList, eps=newEps)

} else if (aID > 5 & aID < 10) {
  # replication 
  initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
  initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                     matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
  newRes <- symm_fit_ind_EM(testStats = testDat[, 1:2], sameDirAlt = TRUE, initMuList = initMuList, initPiList = initPiList, eps=newEps)
} else if (aID >= 10 & aID <= 12) {

  # correlated 
  tempCorVal <- cor(testDat[, 1], testDat[, 2])
  tempCorMat <- matrix(data=c(1, tempCorVal, tempCorVal, 1), nrow=2)
  initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
  initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                     matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
  newRes <- symm_fit_cor_EM(testStats = testDat[, 1:2], corMat = tempCorMat, initMuList = initMuList, initPiList = initPiList, eps=newEps)
} else if (aID > 12) {
  # 3D independent 
  initPiList <- list(c(0.82))
  for (i in 2:7) {initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)}
  initPiList[[8]] <- c(0.1)
  tempH <- expand.grid(c(0, 1), c(0, 1), c(0, 1)) %>%
    mutate(s = Var1 + Var2 + Var3) %>%
    arrange(s) %>%
    select(-s) %>%
    as.matrix(.)
  # the symm_fit_ind.R code will add the appropriate 0s to initMuList
  initMuList <- list(matrix(data=0, nrow=3, ncol=1))
  for (i in 2:7) {
    initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
    #initMuList[[i]] <- cbind(tempH[i, ] * 2, tempH[i, ] * 5)
  }
  initMuList[[8]] <- matrix(data=c(8, 8, 8), nrow=3)

  newRes <- symm_fit_ind_EM(testStats = testDat[, 1:3], initMuList = initMuList, initPiList = initPiList, eps = newEps)
  
}

# end timing
endTime <- proc.time()[[3]]
diffTime <- endTime - startTime

# save
setwd(outputDir)
write.table(diffTime, paste0("analysis_timing_num_", num, "_aID", aID, ".txt"), append=F, quote=F, row.names=F, col.names=F) 








