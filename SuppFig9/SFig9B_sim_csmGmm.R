# For Supp Fig 9B

# load libraries
library(mvtnorm)
library(data.table)
library(bindata)
library(dplyr)
library(magrittr)
library(devtools)
devtools::install.packages("ryanrsun/csmGmm")
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

#------------------------------------------------------------------#
# parameters to be changed
# set output directory 
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig9/output"
outName <- paste0("sim_n1k_j100k_ind5d_changepi0_new_aID", aID, ".txt")

# option to save or load intermediate data to save time,
# set as FALSE for first run and then TRUE thereafter
loadData <- TRUE
saveData <- FALSE
# the name will be [testStatsName]_[betaStart]_S[Snum]_aID[aID].txt
testStatsName <- "allZ_changepi0"
betaName <- "allBeta_changepi0"
#-------------------------------------------------------------------#

# parameters
doHDMT <- FALSE
doDACT <- FALSE
doKernel <- FALSE
do50df <- FALSE
do7df <- FALSE
doNew <- TRUE
qvalue <- 0.1
nSNPs <- 10^5 
setSize <- 1000
n <- 1000
nDims <- 5
nSets <- nSNPs / setSize
nSims <- 1
margprob <- rep(0.3, setSize)
simsPerEffSize <- 100
effSizeMult <- ceiling(aID / simsPerEffSize)
betaMin <- rep(0.38, nDims)
betaMax <- rep(0.38, nDims)
beta0 <- -1

# determines how many signals there are
pi1 <- 0.01 * effSizeMult
pi2 <- 0.001 * effSizeMult
pi3 <- pi2 / 3
pi0 <- 1 - pi1 - pi2 - pi3 - pi3 - pi3
sProp <- c(pi0, pi1, pi2, pi3, pi3, pi3)
hMat <- expand.grid(rep(list(c(-1, 0, 1)), nDims)) %>%
  as.matrix(.) %>%
  cbind(., rowSums(abs(.))) %>%
  as.data.frame(.) %>%
  set_colnames(c(paste0("Var", 1:(ncol(.)-1)), "s")) %>%
  arrange(s)
number <- c()
for (s_it in 0:max(hMat$s)) {
  numRows <- length(which(hMat$s == s_it))
  number <- c(number, rep(sProp[s_it + 1] * nSNPs / numRows, numRows))
}
hMat <- hMat %>% mutate(number = number)

# record results here
powerRes <- data.frame(nCausal=rep(NA, nSims),  minEff1=betaMin[1],
                       seed=NA, pi0aTrue=NA, pi0bTrue=NA,
                       powerDACT=NA, powerHDMT=NA, powerKernel=NA, power7df=NA,
                       power50df=NA, powerNew=NA, fdpDACT=NA, fdpHDMT=NA, fdpKernel=NA, fdp7df=NA, fdp50df=NA, fdpNew=NA,
                       inconKernel=NA, incon7df=NA, incon50df=NA, inconNew=NA)
# each loop is one simulation iteration
for (sim_it in 1:nSims) {

  # set the seed 
  set.seed(aID * 10^5 + sim_it)
  powerRes$seed[sim_it] <- aID * 10^5 + sim_it

  # load or save data
  if (loadData) {
    setwd(outputDir)
    allZ <- fread(paste0(testStatsName, "_S", Snum, "_aID", aID,  ".txt"), data.table=F)
    allBeta <- fread(paste0(betaName, "_S", Snum, "_aID", aID, ".txt"), data.table=F)
  } else {
    # hold test statistics and signals
    allZ <- matrix(data=NA, nrow=nSNPs, ncol=nDims)
    allBeta <- matrix(data=NA, nrow=nSNPs, ncol=nDims)

    # select signal locations
    sigLocsMat <- allocate_sigs(hMat = hMat, nSNPs = nSNPs, nDims=nDims)
  
    # generate data in multiple sets - faster than all at once
    for (set_it in 1:nSets)  {
      # save test statistics      
      statsMat <- matrix(data=NA, nrow=setSize, ncol=nDims)

      # generate coefficient matrix
      coefMat <- set_beta(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize, 
                       betaMin=betaMin, betaMax=betaMax)

      # loop through each dimension
      for (dimension_it in 1:nDims) {
        tempG <- sapply(X=margprob, FUN=rbinom, n=n, size=2)
        tempCoef <- coefMat[, dimension_it]

        # outcome
        tempEta <- sweep(tempG, MARGIN=2, STATS=tempCoef, FUN="*") + matrix(data=beta0, nrow=nrow(tempG), ncol=ncol(tempG))
        tempMu <- rje::expit(as.numeric(tempEta))
        tempY <- rbinom(n=length(tempMu), size=1, prob=tempMu)
        # put back in matrix form
        yMat <- matrix(data=tempY, nrow=nrow(tempEta), ncol=ncol(tempEta), byrow=FALSE)

        # loop through each SNP
        for (test_it in 1:ncol(yMat)) {
          tempMod <- glm(yMat[, test_it] ~ tempG[, test_it], family=binomial)
          statsMat[test_it, dimension_it] <- summary(tempMod)$coefficients[2, 3]
        }
        cat("simulated dimension ", dimension_it, "\n")
      }

      # record data
      startIdx <- (set_it - 1) * setSize + 1
      endIdx <- set_it * setSize
      allZ[startIdx:endIdx, ] <- statsMat
      allBeta[startIdx:endIdx, ] <- coefMat
   
      # checkpoint 
      if(set_it%%1 == 0) {cat(set_it)}
    } # done generating data

    # save it
    if (saveData) { 
      setwd(outputDir)
      write.table(allZ, paste0(testStatsName, "_S", Snum, "_aID", aID, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
      write.table(allBeta, paste0(betaName, "_S", Snum, "_aID", aID, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    } 
  }

  # number of signals and causal SNPs
  causalVec <- as.numeric(allBeta[, 1] != 0 & allBeta[, 2] != 0)
  powerRes$nCausal[sim_it] <- sum(causalVec)
  powerRes$pi0aTrue[sim_it] <- length(which(allBeta[, 1] != 0))
  powerRes$pi0bTrue[sim_it] <- length(which(allBeta[, 2] != 0))

  # adjustment to not get p-values of 0 needed for DACT and HDMT 
  for (col_it in 1:ncol(allZ)) {
    tooBig <- which(allZ[, col_it] > 8.1)
    tooSmall <- which(allZ[, col_it] < -8.1)
    if (length(tooBig) > 0) {
      allZ[tooBig, col_it] <- 8.1
    }
    if (length(tooSmall) > 0) {
      allZ[tooSmall, col_it] <- -8.1
    }
  }
  # p-value matrix
  allP <- 1- pchisq(as.matrix(allZ)^2, df=1)
 
  # hold the results
  totOut <- data.frame(X1 = allP[, 1], X2 = allP[, 2], origIdx = 1:nrow(allP), causal=causalVec) %>%
          mutate(pmax = pmax(X1, X2))

  # analyze it 
  # start with HDMT
  nullprop <- tryCatch(null_estimation(allP), error=function(e) e, warning=function(w) w)
  if (class(nullprop)[1] == "list" & doHDMT) {
    hdmtRes <- tryCatch(fdr_est_orig(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,allP),
                     error = function(e) e, warning = function(w) w)
    if (class(hdmtRes)[1] == "numeric") {

      totOut <- totOut %>% mutate(hdmtRes = hdmtRes)

      # calculate power and fdp for HDMT
      powerRes$fdpHDMT[sim_it] <- length(which(totOut$hdmtRes < qvalue & totOut$causal == 0)) /  length(which(totOut$hdmtRes < qvalue))
      powerRes$powerHDMT[sim_it] <- length(which(totOut$hdmtRes < qvalue & totOut$causal == 1)) / sum(causalVec)
      powerRes$nRejHDMT[sim_it] <- length(which(totOut$hdmtRes < qvalue)) 
    } else {totOut <- totOut %>% mutate(hdmtRes = NA)}
  } else {totOut <- totOut %>% mutate(hdmtRes = NA)}
  
  # DACT
  if (class(nullprop)[1] != "list") {
    nullprop <- NULL
  }
  if (doDACT) {
    DACTout <- tryCatch(DACT_noEst(p_a = allP[, 1], p_b = allP[, 2], nullEst=nullprop, correction="JC"),
                        error = function(e) e, warning=function(w) w)
    if (class(DACTout)[1] %in% c("simpleError", "simpleWarning")) {
      totOut <-  totOut %>% mutate(DACTp = NA)
    } else {
      DACTdf <- data.frame(DACTp = DACTout$pval, origIdx = 1:length(DACTout$pval)) %>%
        arrange(DACTp) %>%
        mutate(rankedIdxP = 1:nrow(.)) %>%
        mutate(km = 1:nrow(.)/ nrow(.)) %>%
        mutate(RHS = km * qvalue) 
      rejected <- which(DACTdf$DACTp <= DACTdf$RHS)
      if (length(rejected) == 0) {
        maxIdx <- 0
      } else {maxIdx <- max(rejected)}
      DACTdf <- DACTdf %>% mutate(reject = ifelse(rankedIdxP <= maxIdx, 1, 0)) %>%
        arrange(origIdx)
      # append results
      totOut <- totOut %>% mutate(DACTp = DACTdf$DACTp, DACTrej = DACTdf$reject)

      # power and FDP for DACT
      powerRes$fdpDACT[sim_it] <- length(which(totOut$DACTrej == 1 & totOut$causal == 0)) / length(which(totOut$DACTrej == 1))
      powerRes$powerDACT[sim_it] <- length(which(totOut$DACTrej == 1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$nRejDACT[sim_it] <- length(which(totOut$DACTrej == 1)) 
    }
  } else {totOut <-  totOut %>% mutate(DACTp = NA, DACTrej = NA)}

  # old method - kernel
  if (doKernel) { 
    oldResKernel <- emp_bayes_framework(summary_tab = allZ, kernel = TRUE, joint=FALSE, ind = TRUE, 
                                        dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    if (class(oldResKernel)[1] == "list") {
      totOut <- totOut %>% mutate(kernelLfdr = oldResKernel$lfdrVec) %>%
        arrange(kernelLfdr) %>%
        mutate(kernelAvg = cummean(kernelLfdr)) %>% # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
      powerRes$fdpKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue & totOut$causal == 0)) / length(which(totOut$kernelAvg < qvalue))
      powerRes$powerKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue & totOut$causal == 1)) / sum(causalVec)
      powerRes$nRejKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue)) 
    } else {
    totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)
    }
  } else {totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)}
  
  # old method - 7 df
  if (do7df) { 
    oldRes7df <- emp_bayes_framework(summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE, 
                                      dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    if (class(oldRes7df)[1] == "list") {
      totOut <- totOut %>% mutate(df7Lfdr = oldRes7df$lfdrVec) %>%
        arrange(df7Lfdr) %>%
        mutate(df7Avg = cummean(df7Lfdr)) %>% 
        # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
      powerRes$fdp7df[sim_it] <- length(which(totOut$df7Avg < qvalue & totOut$causal == 0)) /  length(which(totOut$df7Avg < qvalue))
      powerRes$power7df[sim_it] <- length(which(totOut$df7Avg < qvalue & totOut$causal == 1)) / sum(causalVec)
      powerRes$nRej7df[sim_it] <- length(which(totOut$df7Avg < qvalue))
    } else {
      totOut <- totOut %>% mutate(df7Lfdr = NA, df7Avg=NA)
    } 
  } else {totOut <- totOut %>% mutate(df7Lfdr = NA, df7Avg=NA)}
  
  # old method - 50 df
  if (do50df) {
    oldRes50df <- emp_bayes_framework(summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE, 
                                      dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    if (class(oldRes50df)[1] == "list") {
      totOut <- totOut %>% mutate(df50Lfdr = oldRes50df$lfdrVec) %>%
        arrange(df50Lfdr) %>%
        mutate(df50Avg = cummean(df50Lfdr)) %>% 
        # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
      powerRes$fdp50df[sim_it] <- length(which(totOut$df50Avg < qvalue & totOut$causal == 0)) /  length(which(totOut$df50Avg < qvalue))
      powerRes$power50df[sim_it] <- length(which(totOut$df50Avg < qvalue & totOut$causal == 1)) / sum(causalVec)
      powerRes$nRej50df[sim_it] <- length(which(totOut$df50Avg < qvalue)) 
    } else {
      totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)
    }
  } else {totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)}
  
  # new method
  if (doNew) {
    # timing
    startTimeNew <- Sys.time()

    initPiList <- list(c(0.82))
    for (i in 2:(2^nDims)) {initPiList[[i]] <- 0.18 / (2^nDims - 1)}
    # the symm_fit_ind.R code will add the appropriate 0s to initMuList
    initMuList <- list(matrix(data=rep(0, nDims), nrow=nDims, ncol=1))
    for (i in 2:(2^nDims)) {
      initMuList[[i]] <- matrix(data=rep(3, nDims), nrow=nDims, ncol=1)
    }
    newRes <- symm_fit_ind_EM(testStats = allZ, initMuList = initMuList, initPiList = initPiList, eps=10^(-1))
  
    # record
    endTimeNew <- Sys.time() 
    totOut <- totOut %>% mutate(newLfdr = newRes$lfdrResults) %>%
      arrange(newLfdr) %>%
      mutate(newAvg = cummean(newLfdr)) %>% 
      # don't forget to set it back to original index for further additions!!
      arrange(origIdx)
    powerRes$fdpNew[sim_it] <- length(which(totOut$newAvg < qvalue & totOut$causal == 0)) /  length(which(totOut$newAvg < qvalue))
    powerRes$powerNew[sim_it] <- length(which(totOut$newAvg < qvalue & totOut$causal == 1)) / sum(causalVec)
    powerRes$nRejNew[sim_it] <- length(which(totOut$newAvg < qvalue)) 
    powerRes$timeNew[sim_it] <- endTimeNew - startTimeNew 
  }
  cat('\n Done with ', sim_it, '\n')
}

setwd(outputDir)
write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')


