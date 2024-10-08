# summarize timing experiments

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppTab2/summarize_timing.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppTab2/summarize_timing.R")

library(dplyr)
library(data.table)
library(ks)
library(csmGmm)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("SuppTab2", "output/")

# pleiotropy 2d analysis
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "analysis_timing_num_", 1:100, "_aID1.txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
  tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
  fullRes <- rbindlist(tempList)
}
dim(fullRes)
pleio2dAnal <- mean(unlist(fullRes))

# pleiotropy 3d analysis
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "analysis_timing_num_", 1:100, "_aID14.txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
  tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
  fullRes <- rbindlist(tempList)
}
dim(fullRes)
pleio3dAnal <- mean(unlist(fullRes))

# replication 2d analysis
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "analysis_timing_num_", 1:100, "_aID8.txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
  tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
  fullRes <- rbindlist(tempList)
}
dim(fullRes)
rep2dAnal <- mean(unlist(fullRes))

# correlated 2d analysis
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "Fig3C_timing_aID", 1:100, ".txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
  tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
  fullRes <- rbindlist(tempList)
}
dim(fullRes)
cor2dAnal <- mean(unlist(fullRes))

# med 2d simulation fit 1
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "SFig26A_aID", 1:100, "_fit1_timing.txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
  tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
  fullRes <- rbindlist(tempList)
} 
dim(fullRes)
med2dSim_fit1 <- mean(unlist(fullRes$diffTime))

# med 2d simulation fit 3
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "SFig26A_aID", 1:100, "_fit3_timing.txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
    tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
    fullRes <- rbindlist(tempList)
}
dim(fullRes)
med2dSim_fit3 <- mean(unlist(fullRes$diffTime))


# med 2d simulation fit 5
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "SFig26A_aID", 1:100, "_fit5_timing.txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
    tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
    fullRes <- rbindlist(tempList)
}
dim(fullRes)
med2dSim_fit5 <- mean(unlist(fullRes$diffTime))

# med 2d simulation fit 7
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "SFig26A_aID", 1:100, "_fit7_timing.txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
    tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
    fullRes <- rbindlist(tempList)
}
dim(fullRes)
med2dSim_fit7 <- mean(unlist(fullRes$diffTime))

# rep 2d simulation 
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "Fig3D_timing_aID", 1:100, ".txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
  tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
  fullRes <- rbindlist(tempList)
}
dim(fullRes)
rep2dSim <- mean(unlist(fullRes$diffTime))


# cor2d simulation 
allFiles <- list.files(outputDir)
expFiles <- paste0(outputDir, "Fig3C_timing_aID", 1:100, ".txt")
foundFiles <- allFiles[which(allFiles %in% expFiles)]
fullRes <- fread(foundFiles[1])
for (file_it in 2:length(foundFiles)) {
    tempRes <- fread(foundFiles[file_it])
  tempList <- list(fullRes, tempRes)
    fullRes <- rbindlist(tempList)
}
dim(fullRes)
cor2dSim <- mean(unlist(fullRes$diffTime))

# put it in a table
timeDF <- data.frame(Method = c("Mediation", "Pleiotropy 3D", "Replication", "Pleiotropy (Correlated)", 
                                "Mediation", "Mediation", "Mediation", "Mediation",
                                "Replication", "Pleiotropy (Correlated)"),
                     nSNPs = c(rep(6753538, 4), rep(100000, 6)),
                     Model = c("csmGmm", "csmGmm", "r-csmGmm", "c-csmGmm", "csmGmm-1", 
                               "csmGmm-3", "csmGmm-5", "csmGmm-7", "r-cspGmm", "c-csmGmm"),
                     Time = c(pleio2dAnal, pleio3dAnal, rep2dAnal, cor2dAnal, med2dSim_fit1,
                              med2dSim_fit3, med2dSim_fit5, med2dSim_fit7, rep2dSim, cor2dSim))

library(xtable)
print(xtable(timeDF, digits=0), include.rownames=F)












