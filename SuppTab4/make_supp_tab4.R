# Supplementary Table 4

# load libraries
library(data.table)
library(dplyr)
library(devtools)
library(ks)
devtools::install.packages("ryanrsun/csmGmm")
library(csmGmm)
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)
#------------------------------------------------------------------#
# parameters to be changed
# set output directory 
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/Fig4/output"
fnameOut <- "processed_ukb_data"
rejectOutRoot <- "reject_bmi_with_overall_neg5_reject_aID"
#------------------------------------------------------------------#


# how the raw output files are named
fnameRoot <- paste0("Fig4_data_aID", 1:9, ".txt")
fnameDACT <- paste0(fnameRoot, "_DACTp.txt")
fnameHDMT <- paste0(fnameRoot, "_hdmt.txt")
fnameKernel <- paste0(fnameRoot, "_kernel.txt")
fname7 <- paste0(fnameRoot, "_df7.txt")
fname50 <- paste0(fnameRoot, "df50.txt")
fnameNew <- c(paste0(fnameRoot, "_newlfdr.txt"))

selections <- list()
selections[[1]] <- c("Zcad", "Zbmi")
selections[[2]] <- c("Zoverall", "Zcad")
selections[[3]] <- c("Zoverall", "Zbmi")
selections[[4]] <- c("Zoverall", "Zlcukb")
selections[[5]] <- c("Zcad_cardio", "Zcadukb")
selections[[6]] <- c("Zoverall", "Zcad", "Zbmi")
selections[[7]] <- c("Zlcukb", "Zcadukb")
selections[[8]] <- c("Zlcukb", "Zbmi")
selections[[9]] <- c("Zcadukb", "Zbmi")


# results
allResults <- c()
for (file_it in 6:6) {

  # open data
  cleanZ <- fread("bmi_with_overall.txt")
  fdrLimitHDMT <- fdrLimitHDMTi
  fdrLimitDACT <- fdrLimitDACTi
  fdrLimitKernel <- fdrLimitKerneli
  fdrLimit50 <- fdrLimit50i
  fdrLimit7 <- fdrLimit7i  
 
  # hold temporary results
  tempRes <- data.frame(Method=c("New", "Kernel", "df50", "df7", "HDMT", "DACT"), numReject=NA)

  # new 
  tempNew <- fread(fnameNew[file_it], header=T, data.table=F)
  tempDat <- cleanZ %>% select(all_of(selections[[file_it]]), chrpos) %>%
    mutate(origIdx = 1:nrow(.)) %>%   
    mutate(newLfdr = tempNew$x) %>%
    arrange(newLfdr) %>%
    mutate(cumNew = cummean(newLfdr)) %>%
    mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[1] <- sum(tempDat$rejNew)
  
  # kernel
  tempKernel <- fread(fnameKernel[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(kernelLfdr = tempKernel$x) %>%
    arrange(kernelLfdr) %>%
    mutate(cumKernel = cummean(kernelLfdr)) %>%
    mutate(rejKernel = ifelse(cumKernel < fdrLimitKernel, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[2] <- sum(tempDat$rejKernel)

  # df50
  tempdf50 <- fread(fname50[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(df50Lfdr = tempdf50$x) %>%
    arrange(df50Lfdr) %>%
    mutate(cumdf50 = cummean(df50Lfdr)) %>%
    mutate(rejdf50 = ifelse(cumdf50 < fdrLimit50, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[3] <- sum(tempDat$rejdf50)

  # df7
  tempdf7 <- fread(fname7[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(df7Lfdr = tempdf7$x) %>%
    arrange(df7Lfdr) %>%
    mutate(cumdf7 = cummean(df7Lfdr)) %>%
    mutate(rejdf7 = ifelse(cumdf7 < fdrLimit7, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[4] <- sum(tempDat$rejdf7)
  
  # any rejection
  tempDat <- tempDat %>% mutate(rejAny = ifelse(rejdf7 == 1 |  rejdf50 == 1 | rejKernel == 1 | rejNew == 1, 1, 0))
  tempDat <- tempDat %>% mutate(minlfdr = pmin(kernelLfdr, df50Lfdr, df7Lfdr))
  
}

# display
tempDat <- tempDat %>% arrange(minlfdr) %>%
  select(chrpos, Zoverall, Zcad, Zbmi, kernelLfdr, df50Lfdr, df7Lfdr)
head(tempDat, 10)

