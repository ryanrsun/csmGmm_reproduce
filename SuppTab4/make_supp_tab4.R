# Supplementary Table 4

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppTab4/make_supp_tab4.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppTab4/make_supp_tab4.R")

# load libraries
library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)
library(here)

outputDir <- here::here("Fig4", "output")
dataDir <- here::here("Data")

fnameRoot <- paste0(outputDir, "/Fig4_data_aID", 1:9)
rejectOutRoot <- "reject_bmi_with_overall_neg5_reject_S_1_aID"

# how the raw output files are named
fnameDACT <- paste0(fnameRoot, "_DACTp.txt")
fnameHDMT <- paste0(fnameRoot, "_hdmt.txt")
fnameKernel <- paste0(fnameRoot, "_kernel.txt")
fname7 <- paste0(fnameRoot, "_df7.txt")
fname50 <- paste0(fnameRoot, "_df50.txt")
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
  cleanZ <- fread(paste0(dataDir, "/bmi_with_overall.txt"))
  fdrLimitHDMT <- 0.01
  fdrLimitDACT <- 0.01
  fdrLimitKernel <- 0.01
  fdrLimit50 <- 0.01
  fdrLimit7 <- 0.01  
  fdrLimitNew <- 0.1
 
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

