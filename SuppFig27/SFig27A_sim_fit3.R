# For Supp Figure 27A

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig27/SFig27A_sim_fit3.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig27/SFig27A_sim_fit3.R")

# load libraries
library(mvtnorm)
library(data.table)
library(bindata)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("SuppFig27", "output")
outName <- paste0(outputDir, "/SFig27A_aID", aID, "_fit3.txt")

# option to save or load intermediate data to save time
loadData <- FALSE
saveData <- TRUE
testStatsName <- here::here(outputDir, "SFig27Af3_allZ")
betaName <- here::here(outputDir, "SFig27Af3_allBeta")

# parameters
doHDMT <- FALSE
doDACT <- FALSE
doKernel <- FALSE
do50df <- FALSE
do7df <- FALSE
doNew <- TRUE
qvalue <- 0.1
nSNPs <- 3 * 10^5 
setSize <- 1000
n <- 1000
nDims <- 2
nSets <- nSNPs / setSize
nSims <- 2
margprob <- rep(0.3, setSize)
simsPerEffSize <- 50
effSizeMult <- ceiling(aID / simsPerEffSize)
betaOptsNon <- cbind(c(0 + effSizeMult * 0.01, 0.02 + effSizeMult * 0.01, 0.04 + effSizeMult * 0.01, 0.06 + effSizeMult * 0.01, 0.08 + effSizeMult * 0.01),
                     c(0 + effSizeMult * 0.01, 0.02 + effSizeMult * 0.01, 0.04 + effSizeMult * 0.01, 0.06 + effSizeMult * 0.01, 0.08 + effSizeMult * 0.01))
betaOptsCausal <- cbind(c(0.08 + effSizeMult * 0.01, 0.1 + effSizeMult * 0.01, 0.12 + effSizeMult * 0.01, 0.14 + effSizeMult * 0.01, 0.16 + effSizeMult * 0.01),
                        c(0.08 + effSizeMult * 0.01, 0.1 + effSizeMult * 0.01, 0.12 + effSizeMult * 0.01, 0.14 + effSizeMult * 0.01, 0.16 + effSizeMult * 0.01))
beta0 <- -1

# allocate snps new function for set mbl
allocate_sigs_mbl <- function(nSNPs) {

  # we're going to do 200k SNPs
  # 20% null, the rest non-null
  totSigs <- nSNPs * 0.8
  # sample the locations
  tempPos <- sample(x=1:nSNPs, size=totSigs, replace=F)

  # will put the signal locations here
  sigLocsMat <- matrix(data=0, nrow=length(tempPos), ncol=nDims * 2)
  # 20% empty, 10% for each of other configurations, divide by 3 for three true means
  twoPct <- nSNPs / 50
  # first the (1,0) and (-1, 0)
  sigLocsMat[1:(10 * twoPct), 1] <- tempPos[1:(10 * twoPct)]
  sigLocsMat[1:twoPct, 3] <- 1
  sigLocsMat[(1 * twoPct + 1):(2 * twoPct), 3] <- 2
  sigLocsMat[(2 * twoPct + 1):(3 * twoPct), 3] <- 3
  sigLocsMat[(3 * twoPct + 1):(4 * twoPct), 3] <- 4
  sigLocsMat[(4 * twoPct + 1):(5 * twoPct), 3] <- 5
  sigLocsMat[(5 * twoPct + 1):(6 * twoPct), 3] <- -1
  sigLocsMat[(6 * twoPct + 1):(7 * twoPct), 3] <- -2
  sigLocsMat[(7 * twoPct + 1):(8 * twoPct), 3] <- -3
  sigLocsMat[(8 * twoPct + 1):(9 * twoPct), 3] <- -4
  sigLocsMat[(9 * twoPct + 1):(10 * twoPct), 3] <- -5
  # now the (0, 1) and (0, -1)
  sigLocsMat[(10 * twoPct + 1):(20 * twoPct), 2] <- tempPos[(10 * twoPct + 1):(20 * twoPct)]
  sigLocsMat[(10 * twoPct + 1):(11 * twoPct), 4] <- 1
  sigLocsMat[(11 * twoPct + 1):(12 * twoPct), 4] <- 2
  sigLocsMat[(12 * twoPct + 1):(13 * twoPct), 4] <- 3
  sigLocsMat[(13 * twoPct + 1):(14 * twoPct), 4] <- 4
  sigLocsMat[(14 * twoPct + 1):(15 * twoPct), 4] <- 5
  sigLocsMat[(15 * twoPct + 1):(16 * twoPct), 4] <- -1
  sigLocsMat[(16 * twoPct + 1):(17 * twoPct), 4] <- -2
  sigLocsMat[(17 * twoPct + 1):(18 * twoPct), 4] <- -3
  sigLocsMat[(18 * twoPct + 1):(19 * twoPct), 4] <- -4
  sigLocsMat[(19 * twoPct + 1):(20 * twoPct), 4] <- -5

  # now the (1, 1) and (1, -1)
  sigLocsMat[(20 * twoPct + 1):(30 * twoPct), 1:2] <- tempPos[(20 * twoPct + 1):(30 * twoPct)]
  sigLocsMat[(20 * twoPct + 1):(21 * twoPct), 3:4] <- 1
  sigLocsMat[(21 * twoPct + 1):(22 * twoPct), 3:4] <- 2
  sigLocsMat[(22 * twoPct + 1):(23 * twoPct), 3:4] <- 3
  sigLocsMat[(23 * twoPct + 1):(24 * twoPct), 3:4] <- 4
  sigLocsMat[(24 * twoPct + 1):(25 * twoPct), 3:4] <- 5
  sigLocsMat[(25 * twoPct + 1):(26 * twoPct), 3] <- 1
  sigLocsMat[(25 * twoPct + 1):(26 * twoPct), 4] <- -1
  sigLocsMat[(26 * twoPct + 1):(27 * twoPct), 3] <- 2
  sigLocsMat[(26 * twoPct + 1):(27 * twoPct), 4] <- -2
  sigLocsMat[(27 * twoPct + 1):(28 * twoPct), 3] <- 3
  sigLocsMat[(27 * twoPct + 1):(28 * twoPct), 4] <- -3
  sigLocsMat[(28 * twoPct + 1):(29 * twoPct), 3] <- 4
  sigLocsMat[(28 * twoPct + 1):(29 * twoPct), 4] <- -4
  sigLocsMat[(29 * twoPct + 1):(30 * twoPct), 3] <- 5
  sigLocsMat[(29 * twoPct + 1):(30 * twoPct), 4] <- -5

  # (now the (-1, 1) and (-1, -1)
  sigLocsMat[(30 * twoPct + 1):(40 * twoPct), 1:2] <- tempPos[(30 * twoPct + 1):(40 * twoPct)]
  sigLocsMat[(30 * twoPct + 1):(31 * twoPct), 3:4] <- -1
  sigLocsMat[(31 * twoPct + 1):(32 * twoPct), 3:4] <- -2
  sigLocsMat[(32 * twoPct + 1):(33 * twoPct), 3:4] <- -3
  sigLocsMat[(33 * twoPct + 1):(34 * twoPct), 3:4] <- -4
  sigLocsMat[(34 * twoPct + 1):(35 * twoPct), 3:4] <- -5
  sigLocsMat[(35 * twoPct + 1):(36 * twoPct), 3] <- -1
  sigLocsMat[(35 * twoPct + 1):(36 * twoPct), 4] <- 1
  sigLocsMat[(36 * twoPct + 1):(37 * twoPct), 3] <- -2
  sigLocsMat[(36 * twoPct + 1):(37 * twoPct), 4] <- 2
  sigLocsMat[(37 * twoPct + 1):(38 * twoPct), 3] <- -3
  sigLocsMat[(37 * twoPct + 1):(38 * twoPct), 4] <- 3
  sigLocsMat[(38 * twoPct + 1):(39 * twoPct), 3] <- -4
  sigLocsMat[(38 * twoPct + 1):(39 * twoPct), 4] <- 4
  sigLocsMat[(39 * twoPct + 1):(40 * twoPct), 3] <- -5
  sigLocsMat[(39 * twoPct + 1):(40 * twoPct), 4] <- 5


  sigLocsMat[(18 * twoPct + 1):(19 * twoPct), 3:4] <- -1
  sigLocsMat[(19 * twoPct + 1):(20 * twoPct), 3:4] <- -2
  sigLocsMat[(20 * twoPct + 1):(21 * twoPct), 3:4] <- -3
  sigLocsMat[(21 * twoPct + 1):(22 * twoPct), 3] <- -1
  sigLocsMat[(21 * twoPct + 1):(22 * twoPct), 4] <- 1
  sigLocsMat[(22 * twoPct + 1):(23 * twoPct), 3] <- -2
  sigLocsMat[(22 * twoPct + 1):(23 * twoPct), 4] <- 2
  sigLocsMat[(23 * twoPct + 1):(24 * twoPct), 3] <- -3
  sigLocsMat[(23 * twoPct + 1):(24 * twoPct), 4] <- 3
  sigLocsMat <- data.frame(sigLocsMat)
  colnames(sigLocsMat) <- c("Pos1", "Pos2", "Eff1", "Eff2")
  return(sigLocsMat)
}

# place signals new function for set mbl
set_beta_mbl <- function(sigLocsMat, set_it, setSize, randomBeta, betaOptsNon=NULL, betaOptsCausal=NULL, betaMin=NULL, betaMax=NULL) {

  # return this 
  betaMat <- matrix(data=0, nrow=setSize, 2)

  # indices of this set 
  tempIdx <- ((set_it - 1) * setSize + 1):(set_it * setSize)

  # causal and noncausal
  oneOnly <- which(tempIdx %in% sigLocsMat[, 1] & !(tempIdx %in% sigLocsMat[, 2]))
  twoOnly <- which(tempIdx %in% sigLocsMat[, 2] & !(tempIdx %in% sigLocsMat[, 1]))
  bothIn <- which(tempIdx %in% sigLocsMat[, 1] & tempIdx %in% sigLocsMat[, 2])

  if (length(oneOnly) > 0) {
    tempSigLoc <- data.frame(Pos1 = tempIdx[oneOnly])
    tempMerge <- merge(tempSigLoc, sigLocsMat, by="Pos1") %>% arrange(Pos1)
    tempDir <- sign(tempMerge$Eff1)
    tempMagIdx <- abs(tempMerge$Eff1)
    tempMag <- betaOptsNon[tempMagIdx, 1]
    betaMat[oneOnly, 1] <- tempMag * tempDir
  }

  if (length(twoOnly) > 0) {
    tempSigLoc <- data.frame(Pos2 = tempIdx[twoOnly])
    tempMerge <- merge(tempSigLoc, sigLocsMat, by="Pos2")
    tempDir <- sign(tempMerge$Eff2)
    tempMagIdx <- abs(tempMerge$Eff2)
    tempMag <- betaOptsNon[tempMagIdx, 2]
    betaMat[twoOnly, 2] <- tempMag * tempDir
  }

  if (length(bothIn) > 0) {
    tempSigLoc <- data.frame(Pos1 = tempIdx[bothIn])
    tempMerge <- merge(tempSigLoc, sigLocsMat, by="Pos1")
    tempDir1 <- sign(tempMerge$Eff1)
    tempDir2 <- sign(tempMerge$Eff2)
    tempMagIdx1 <- abs(tempMerge$Eff1)
    tempMagIdx2 <- abs(tempMerge$Eff2)
    tempMag1 <- betaOptsCausal[tempMagIdx1, 1]
    tempMag2 <- betaOptsCausal[tempMagIdx2, 2]
    betaMat[bothIn, 1] <- tempMag1 * tempDir1
    betaMat[bothIn, 2] <- tempMag2 * tempDir2
  }
  # return
  return(betaMat)
}



# record results here
powerRes <- data.frame(nCausal=rep(NA, nSims),  minEff1=betaOptsNon[1,1],
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
    allZ <- fread(paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
    allBeta <- fread(paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
  } else {
    # hold test statistics and signals
    allZ <- matrix(data=NA, nrow=nSNPs, ncol=2)
    allBeta <- matrix(data=NA, nrow=nSNPs, ncol=2)

    # select signal locations
    sigLocsMat <- allocate_sigs_mbl(nSNPs = nSNPs)
    
    # generate data in multiple sets - faster than all at once
    for (set_it in 1:nSets)  {

      # generate coefficient matrix
      coefMat <- set_beta_mbl(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize, randomBeta = randomBeta, betaOptsNon=betaOptsNon, betaOptsCausal=betaOptsCausal)

      # two dimensional mediation case data generation for k=1 dimension
      tempG <- sapply(X=margprob, FUN=rbinom, n=n, size=2)
      adjG <- sweep(tempG, MARGIN=2, STATS=apply(tempG, 2, mean), FUN="-")
      tempAlpha <- coefMat[, 1] 
      tempM <- sweep(tempG, MARGIN=2, STATS=tempAlpha, FUN="*") + matrix(data = rnorm(n=n*nrow(coefMat)), nrow=n, ncol=nrow(coefMat))
      adjM <- sweep(tempM, MARGIN=2, STATS=apply(tempM, 2, mean), FUN="-")
      sigSqHat <- apply(adjM, 2, myvar_fun)
 
      # calculate test statistics  for k=1
      statsMat <- matrix(data=NA, nrow=setSize, ncol=2) 
      tempNum <- apply(adjG * adjM, 2, sum)
      tempDenom <- sqrt(apply(adjG^2, 2, sum) * sigSqHat)
      statsMat[, 1] <- tempNum / tempDenom

      # mediation data generation for k=2 dimension 
      tempBeta <- coefMat[, 2]
      tempEta <- sweep(tempM, MARGIN=2, STATS=tempBeta, FUN="*") + matrix(data=beta0, nrow=nrow(tempM), ncol=ncol(tempM))
      tempMu <- rje::expit(as.numeric(tempEta))
      tempY <- rbinom(n=length(tempMu), size=1, prob=tempMu)
      yMat <- matrix(data=tempY, nrow=nrow(tempEta), ncol=ncol(tempEta), byrow=FALSE)

      # calculate test statistics for k=2 
      for (test_it in 1:ncol(yMat)) {
        tempMod <- glm(yMat[, test_it] ~ tempG[, test_it] + tempM[, test_it], family=binomial)
        statsMat[test_it, 2] <- summary(tempMod)$coefficients[3, 3]
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
      write.table(allZ, paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
      write.table(allBeta, paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
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
      powerRes$inconKernel[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec))
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
      powerRes$incon7df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes7df$lfdrVec))
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
      powerRes$incon50df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes50df$lfdrVec))
      powerRes$nRej50df[sim_it] <- length(which(totOut$df50Avg < qvalue)) 
    } else {
      totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)
    }
  } else {totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)}
  
  # new method
  if (doNew) {
    initPiList <- list(c(0.2), c(0.2/3, 0.2/3, 0.2/3), c(0.2/3, 0.2/3, 0.2/3), c(0.4/3, 0.4/3, 0.4/3))
    initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 1, 0, 2, 0, 3), nrow=2),
                     matrix(data=c(1, 0, 2, 0, 3, 0), nrow=2), matrix(data=c(3, 3, 4, 4, 5, 5), nrow=2)) 
    newRes <- symm_fit_ind_EM(testStats = allZ, initMuList = initMuList, initPiList = initPiList, eps=10^(-5))
    # record
    totOut <- totOut %>% mutate(newLfdr = newRes$lfdrResults) %>%
      arrange(newLfdr) %>%
      mutate(newAvg = cummean(newLfdr)) %>% 
      # don't forget to set it back to original index for further additions!!
      arrange(origIdx)
    powerRes$fdpNew[sim_it] <- length(which(totOut$newAvg < qvalue & totOut$causal == 0)) /  length(which(totOut$newAvg < qvalue))
    powerRes$powerNew[sim_it] <- length(which(totOut$newAvg < qvalue & totOut$causal == 1)) / sum(causalVec)
    powerRes$inconNew[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = newRes$lfdrResults))
    powerRes$nRejNew[sim_it] <- length(which(totOut$newAvg < qvalue)) 
  }
  cat('\n Done with ', sim_it, '\n')
}

write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')


