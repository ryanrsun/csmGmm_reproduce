# Process raw analysis of UKB data

# load libraries
library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

#------------------------------------------------------------------#
# parameters to be changed

# source the .R scripts from the supportingCode/ folder in the csmGmm_reproduce repository
setwd('/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SupportingCode/')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

# set output directory 
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/test/output"
fnameOut <- "processed_ukb_data"
rejectOutRoot <- "reject_bmi_with_overall_neg5_reject_aID"
#------------------------------------------------------------------#

# nominal fdr
if (Snum == 1) {
  fdrLimitHDMTi <- 0.1
  fdrLimitDACTi <- 0.01
  fdrLimitKerneli <- 0.01
  fdrLimit7i <- 0.01
  fdrLimit50i <- 0.01
  fdrLimitNew <- 0.1
  fdrLimitHDMTr <- 0.01
  fdrLimitDACTr <- 0.01
  fdrLimitKernelr <- 0.01
  fdrLimit7r <- 0.01
  fdrLimit50r <- 0.01
} else {
  fdrLimitHDMTi <- 0.1
  fdrLimitDACTi <- 0.1
  fdrLimitKerneli <- 0.1
  fdrLimit7i <- 0.1
  fdrLimit50i <- 0.1
  fdrLimitNew <- 0.1
  fdrLimitHDMTr <- 0.1
  fdrLimitDACTr <- 0.1
  fdrLimitKernelr <- 0.1
  fdrLimit7r <- 0.1
  fdrLimit50r <- 0.1
}

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
for (file_it in 1:9) {

  setwd(outputDir)
  if (file_it == 4) {
    cleanZ <- fread("replication_with_lcoverall.txt")
    fdrLimitHDMT <- fdrLimitHDMTr
    fdrLimitDACT <- fdrLimitDACTr
    fdrLimitKernel <- fdrLimitKernelr
    fdrLimit50 <- fdrLimit50r
    fdrLimit7 <- fdrLimit7r
  } else if (file_it == 5) {
    cleanZ <- fread("cad_for_replication.txt")
    fdrLimitHDMT <- fdrLimitHDMTr
    fdrLimitDACT <- fdrLimitDACTr
    fdrLimitKernel <- fdrLimitKernelr
    fdrLimit50 <- fdrLimit50r
    fdrLimit7 <- fdrLimit7r
  } else {
    cleanZ <- fread("bmi_with_overall.txt")
    fdrLimitHDMT <- fdrLimitHDMTi
    fdrLimitDACT <- fdrLimitDACTi
    fdrLimitKernel <- fdrLimitKerneli
    fdrLimit50 <- fdrLimit50i
    fdrLimit7 <- fdrLimit7i  
  }
 
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
  
  if (file_it >= 7) {
    # any rejection
    rejectDat <- tempDat %>% filter(rejNew == 1)

    # allResults
    allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))

    # save
    setwd(outputDir)
    write.table(rejectDat, paste0(rejectOutRoot, file_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    next
  }

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

  if (file_it >= 4 & file_it <= 6) {
    # any rejection
    tempDat <- tempDat %>% mutate(rejAny = ifelse(rejdf7 == 1 |  rejdf50 == 1 | rejKernel == 1 | rejNew == 1, 1, 0))
    rejectDat <- tempDat %>% filter(rejAny == 1)

    # allResults
    allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))
    
    # save
    setwd(outputDir)
    write.table(rejectDat, paste0(rejectOutRoot, file_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    next
  }

  # HDMT
  tempHDMT <- fread(fnameHDMT[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(hdmtFdr = tempHDMT$fixedFdr) %>%
    mutate(rejHDMT = ifelse(hdmtFdr < fdrLimitHDMT, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[5] <- sum(tempDat$rejHDMT)

  # DACT
  tempDACT <- fread(fnameDACT[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(DACTp = tempDACT$x) %>%
    arrange(DACTp) %>%
    mutate(rankedIdxP = 1:nrow(.)) %>%
    mutate(km = 1:nrow(.) / nrow(.)) %>%
    mutate(RHS = km * fdrLimitDACT)
  rejected <- which(tempDat$DACTp <= tempDat$RHS)
  if (length(rejected) == 0) {
    maxIdx <- 0
  } else {maxIdx <- max(rejected)}
  tempDat <- tempDat %>% mutate(rejDACT = ifelse(rankedIdxP <= maxIdx, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[6] <- sum(tempDat$rejDACT)

  # any rejection
  tempDat <- tempDat %>% mutate(rejAny = ifelse(rejDACT == 1 | rejHDMT == 1 | rejdf7 == 1 | 
                                                rejdf50 == 1 | rejKernel == 1 | rejNew == 1, 1, 0))

  rejectDat <- tempDat %>% filter(rejAny == 1) 

  # allResults
  allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))

  # save rejections
  setwd(outputDir)
  write.table(rejectDat, paste0(rejectOutRoot, file_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')

  cat(file_it)
}

# save
setwd(outputDir)
write.table(allResults, paste0(fnameOut, "_S", Snum, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')






