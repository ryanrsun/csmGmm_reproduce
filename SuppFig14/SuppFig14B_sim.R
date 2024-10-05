# For Supp Fig 14B

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig14/SuppFig14B_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig14/SuppFig14B_sim.R")

# load libraries
library(mvtnorm)
library(data.table)
library(bindata)
library(dplyr)
library(magrittr)
library(rje)
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
outputDir <- here::here("SuppFig14", "output")
outName <- paste0(outputDir, "/SFig14B_aID", aID, ".txt")

# option to save or load intermediate data to save time,
# set as FALSE for first run and then TRUE thereafter
loadData <- FALSE
saveData <- TRUE
testStatsName <- here::here(outputDir, "SFig14B_allZ")
betaName <- here::here(outputDir, "SFig14B_allBeta")

# parameters
doHDMT <- TRUE
doDACT <- TRUE
doKernel <- TRUE
do50df <- TRUE
do7df <- TRUE
doNew <- TRUE
qvalue <- 0.1
nSNPs <- 10^5
setSize <- 1000
n <- 1000
nDims <- 3
nSets <- nSNPs / setSize
nSims <- 2
margprob <- rep(0.01, setSize)
simsPerEffSize <- 50
effSizeMult <- ceiling(aID / simsPerEffSize)
betaMin <- c(1.7, 1.7, 1.7)
betaMax <- c(1.7, 1.7, 1.7)
beta0 <- -1

# determines how many signals there are
pi1 <- 0.01 * effSizeMult
pi11 <- 3 * pi1^2
pi111 <- pi11 / 2
pi00 <- 1 - 3 * pi1 - 3 * pi11 - pi111
sProp <- c(pi00, 3 * pi1, 3 * pi11, pi111)
hMat <- expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)) %>%
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
hMat <- hMat %>% mutate(number = number) %>%
  mutate(number = round(number))

# record results here
powerRes <- data.frame(nCausal=rep(NA, nSims), minEff1=pi1,
                       seed=NA, pi0aTrue=NA, pi0bTrue=NA,
                       nRejDACT=NA, nRejHDMT=NA, nRejKernel=NA, nRej7df=NA, nRej50df=NA, nRejNew=NA,
                       powerDACT=NA, powerHDMT=NA, powerKernel=NA, power7df=NA,
                       power50df=NA, powerNew=NA, fdpDACT=NA, fdpHDMT=NA, fdpKernel=NA, fdp7df=NA, fdp50df=NA, fdpNew=NA,
                       inconKernel=NA, incon7df=NA, incon50df=NA, inconNew=NA)
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
    allZ <- matrix(data=NA, nrow=nSNPs, ncol=nDims)
    allBeta <- matrix(data=NA, nrow=nSNPs, ncol=nDims)

    # select signal locations
    sigLocsMat <- allocate_sigs(hMat = hMat, nSNPs = nSNPs, nDims=nDims)
  
    # generate data in multiple sets - faster than all at once 
    for (set_it in 1:nSets)  {

      # information we need to save
      statsMat <- matrix(data=NA, nrow=setSize, ncol=nDims)

      # generate coefficient matrix    
      coefMat <- set_beta(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize, 
                       betaMin=betaMin, betaMax=betaMax) 
  
      ########################################### 
      # three dimensional pleiotropy
      tempG1 <-  sapply(X=margprob, FUN=rbinom, n=n, size=2)
      tempG2 <-  sapply(X=margprob, FUN=rbinom, n=n, size=2)
      tempG3 <- sapply(X=margprob, FUN=rbinom, n=n, size=2) 
      tempAlpha <- coefMat[, 1]
      tempBeta <- coefMat[, 2]
      tempGamma <- coefMat[, 3]
    
      # generate first outcome
      tempEta1 <- sweep(tempG1, MARGIN=2, STATS=tempAlpha, FUN="*") + matrix(data=beta0, nrow=nrow(tempG1), ncol=ncol(tempG1))
      tempMu1 <- rje::expit(as.numeric(tempEta1))
      tempY1 <- rbinom(n=length(tempMu1), size=1, prob=tempMu1)
      yMat1 <- matrix(data=tempY1, nrow=nrow(tempEta1), ncol=ncol(tempEta1), byrow=FALSE)
      # generate second outcome 
      tempEta2 <- sweep(tempG2, MARGIN=2, STATS=tempBeta, FUN="*") + matrix(data=beta0, nrow=nrow(tempG2), ncol=ncol(tempG2))
      tempMu2 <- rje::expit(as.numeric(tempEta2))
      tempY2 <- rbinom(n=length(tempMu2), size=1, prob=tempMu2)
      yMat2 <- matrix(data=tempY2, nrow=nrow(tempEta2), ncol=ncol(tempEta2), byrow=FALSE)
      # generate third outcome
      tempEta3 <- sweep(tempG3, MARGIN=2, STATS=tempGamma, FUN="*") + matrix(data=beta0, nrow=nrow(tempG3), ncol=ncol(tempG3))
      tempMu3 <- rje::expit(as.numeric(tempEta3))
      tempY3 <- rbinom(n=length(tempMu3), size=1, prob=tempMu3)
      yMat3 <- matrix(data=tempY3, nrow=nrow(tempEta3), ncol=ncol(tempEta3), byrow=FALSE)

      # calculate all test statistics 
      for (test_it in 1:ncol(yMat1)) {
        tempMod1 <- glm(yMat1[, test_it] ~ tempG1[, test_it], family=binomial)
        tempMod2 <- glm(yMat2[, test_it] ~ tempG2[, test_it], family=binomial)
        tempMod3 <- glm(yMat3[, test_it] ~ tempG3[, test_it], family=binomial)
        statsMat[test_it, 1] <- summary(tempMod1)$coefficients[2, 3]
        statsMat[test_it, 2] <- summary(tempMod2)$coefficients[2, 3]
        statsMat[test_it, 3] <- summary(tempMod3)$coefficients[2, 3]
      } 
   
      # record
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
  causalVec <- as.numeric(allBeta[, 1] != 0 & allBeta[, 2] != 0 & allBeta[, 3] != 0) 
  powerRes$nCausal[sim_it] <- sum(causalVec)
  powerRes$pi0aTrue[sim_it] <- length(which(allBeta[, 1] != 0))
  powerRes$pi0bTrue[sim_it] <- length(which(allBeta[, 2] != 0))
  powerRes$pi0cTrue[sim_it] <- length(which(allBeta[, 3] != 0))

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
  totOut <- data.frame(X1 = allP[, 1], X2 = allP[, 2], X3 = allP[, 3], origIdx = 1:nrow(allP), causal=causalVec)
 
  # old method - kernel
  if (doKernel) {
    oldResKernel <- emp_bayes_framework(summary_tab = allZ, kernel = TRUE, joint=FALSE, ind = TRUE,
                                        dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    if (class(oldResKernel)[1] == "list") {
      totOut <- totOut %>% mutate(kernelLfdr = oldResKernel$lfdrVec) %>%
        arrange(kernelLfdr) %>%
        mutate(kernelAvg = cummean(kernelLfdr)) %>% 
        arrange(origIdx)
      powerRes$fdpKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 0)) / length(which(totOut$kernelAvg < 0.1))
      powerRes$powerKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$inconKernel[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec))
      powerRes$nRejKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1)) 
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
      powerRes$fdp7df[sim_it] <- length(which(totOut$df7Avg < 0.1 & totOut$causal == 0)) /  length(which(totOut$df7Avg < 0.1))
      powerRes$power7df[sim_it] <- length(which(totOut$df7Avg < 0.1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$incon7df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes7df$lfdrVec))
      powerRes$nRej7df[sim_it] <- length(which(totOut$df7Avg < 0.1))
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
    powerRes$fdp50df[sim_it] <- length(which(totOut$df50Avg < 0.1 & totOut$causal == 0)) /  length(which(totOut$df50Avg < 0.1))
    powerRes$power50df[sim_it] <- length(which(totOut$df50Avg < 0.1 & totOut$causal == 1)) / sum(causalVec)
    powerRes$incon50df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes50df$lfdrVec))
    powerRes$nRej50df[sim_it] <- length(which(totOut$df50Avg < 0.1))
    } else {
      totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)
    }
  } else {totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)}
  # new method
  if (doNew) {
    initPiList <- list(c(0.82))
    for (i in 2:7) {initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)}
    initPiList[[8]] <- c(0.1)
    # the csmGmm package will add the appropriate 0s to initMuList
    initMuList <- list(matrix(data=0, nrow=3, ncol=1))
    for (i in 2:7) {
      initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
    }
    initMuList[[8]] <- matrix(data=c(8, 8, 8), nrow=3)
    newRes <- symm_fit_ind_EM(testStats = allZ, initMuList = initMuList, initPiList = initPiList, eps=10^(-5))

    # record
    totOut <- totOut %>% mutate(newLfdr = newRes$lfdrResults) %>%
      arrange(newLfdr) %>%
      mutate(newAvg = cummean(newLfdr)) %>%
      # don't forget to set it back to original index for further additions!!
      arrange(origIdx)
    powerRes$fdpNew[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$newAvg < 0.1))
    powerRes$powerNew[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
    powerRes$inconNew[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = newRes$lfdrResults))
    powerRes$nRejNew[sim_it] <- length(which(totOut$newAvg < 0.1))
  }
  cat('\n Done with ', sim_it, '\n')
  cat('\n \n \n \n \n \n \n \n \n \n \n')
}

write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')


