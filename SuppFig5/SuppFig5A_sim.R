# Supp Fig 5A

# Note that we cannot share the real data, so we have generated fake data to show
# how this script works. If the fake data is replaced with the real data, then 
# our results will be reproduced.

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig5/SuppFig5A_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig5/SuppFig5A_sim.R")

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
outputDir <- here::here("SuppFig5", "output")
outName <- paste0(outputDir, "/SFig5A_aID", aID, ".txt")

# real genotypes location and names
genotypeDir <- here::here("Data") 
genotypeNames <- paste0(genotypeDir, rep("/cleanG_set1_dataset1.txt", 2000))

# option to save or load intermediate data to save time,
# set as FALSE for first run and then TRUE thereafter
loadData <- FALSE
saveData <- TRUE
# these names are for if saveData <- TRUE
testStatsName <- here::here(outputDir, "SFig5A_allZ")
betaName <- here::here(outputDir, "SFig5A_allBeta")

# parameters
doHDMT <- TRUE
doDACT <- TRUE
doKernel <- TRUE
do50df <- TRUE
do7df <- TRUE
doNew <- TRUE
qvalue <- 0.1
nSNPs <- 10^5
SNPstoGen <- nSNPs + 3 * 10^4 # some SNPs are monomorphic in smaller samples
setSize <- 100
setSizeMult <- 4
setsPerFile <- 2000 / (setSizeMult * setSize) # our genotype files have n=100k people and J=2k SNPs
nDims <- 2
n <- 1000
nSets <- SNPstoGen / setSize
mafCut <- 0.005
nSims <- 10
simsPerEffSize <- 20
effSizeMult <- ceiling(aID / simsPerEffSize)
betaMin <- c(0 + 0.01 * effSizeMult, 0.04 + 0.01 * effSizeMult) 
betaMax <- c(0.1 + 0.01 * effSizeMult, 0.14 + 0.01 * effSizeMult)

# determines how many signals there are
sProp <- c(0.9588, 0.04, 0.0012)
hMat <- expand.grid(c(-1, 0, 1), c(-1, 0, 1)) %>%
  as.data.frame(.) %>%
  mutate(s = abs(Var1) + abs(Var2)) %>%
  arrange(s)
number <- c()
for (s_it in 0:max(hMat$s)) {
  numRows <- length(which(hMat$s == s_it))
  number <- c(number, rep(sProp[s_it + 1] * nSNPs / numRows, numRows))
}
hMat <- hMat %>% mutate(number = number)

# record results here
powerRes <- data.frame(nCausal=rep(NA, nSims),  nSig2=NA, minEff1=betaMin[1],
                       seed=NA, pi0aTrue=NA, pi0bTrue=NA,
                       powerDACT=NA, powerHDMT=NA, powerKernel=NA, power7df=NA,
                       power50df=NA, powerNew=NA, fdpDACT=NA, fdpHDMT=NA, fdpKernel=NA, fdp7df=NA, fdp50df=NA, fdpNew=NA,
                       inconKernel=NA, incon7df=NA, incon50df=NA, inconNew=NA,
                       sig2PowDACT=NA, sig2PowDACTb=NA, sig2PowHDMT=NA, sig2PowKernel=NA, sig2Pow7df=NA, sig2Pow50df=NA, sig2PowNew=NA,
                       sig2fdpDACT=NA, sig2fdpDACTb=NA, sig2fdpHDMT=NA, sig2fdpKernel=NA, sig2fdp7df=NA, sig2fdp50df=NA, sig2fdpNew=NA)
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
    allZ <- matrix(data=NA, nrow=SNPstoGen, ncol=2)
    allNCP <- matrix(data=NA, nrow=SNPstoGen, ncol=2)
    allBeta <- matrix(data=NA, nrow=SNPstoGen, ncol=2)

    # select signal locations
    sigLocsMat <- allocate_sigs(hMat = hMat, nSNPs = nSNPs, nDims=nDims)
  
    # generate data in sets
    for (set_it in 1:nSets)  {

      # information we need to ave
      betaMat <- c()
      statsMat <- c()
      ncpMat <- c() 
  
      # open real genotypes file 
      if (set_it%%setsPerFile == 1) {
        if (set_it > 1) {rm(genoFile)}
        genoFileNum <- ceiling(set_it / setsPerFile)
        genoFile <- fread(genotypeNames[genoFileNum], nrows=nDims*n, header=F) %>% as.matrix(.)
      }

      # cut genotype out of the file
      genoStart <- ((set_it%%setsPerFile) - 1) * (setSize * setSizeMult) + 1
      genoEnd <- (set_it%%setsPerFile) * (setSize * setSizeMult)
      tempG1 <- genoFile[1:n, genoStart:genoEnd]
      tempG2 <- genoFile[(n+1):(2*n), genoStart:genoEnd]
      margprob1 <- apply(tempG1, 2, mean) / 2
      margprob2 <- apply(tempG2, 2, mean) / 2
      margprob <- apply(cbind(margprob1, margprob2), 1, mean)

      # can't do the super rare variants
      keepIdx <- which(margprob > mafCut)
      if (length(keepIdx) < setSize) {next}
      tempG1 <- tempG1[, keepIdx[1:setSize]]
      tempG2 <- tempG2[, keepIdx[1:setSize]]
      margprob = margprob[keepIdx[1:setSize]]
   
      # center genotypes
      adjG1 <- sweep(tempG1, MARGIN=2, STATS=apply(tempG1, 2, mean), FUN="-")
      adjG2 <- sweep(tempG2, MARGIN=2, STATS=apply(tempG2, 2, mean), FUN="-")

      # estimate covariance each time 
      tempCov <- cov(adjG1)
    
      # generate beta matrix
      coefMat <- set_beta(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize, 
                       betaMin=betaMin, betaMax=betaMax) 

      # outcome
      tempY1 <- adjG1 %*% coefMat[, 1] + rnorm(n=n)
      adjY1 <- tempY1 - mean(tempY1)
      tempY2 <- adjG2 %*% coefMat[, 2] + rnorm(n=n)
      adjY2 <- tempY2 - mean(tempY2)

      # test statistics
      tempStats1 <- t(adjG1) %*% adjY1 / sqrt(apply(adjG1, 2, myvar_fun) * n)
      tempStats2 <- t(adjG2) %*% adjY2 / sqrt(apply(adjG2, 2, myvar_fun) * n)
      statsMat <- cbind(tempStats1, tempStats2)

      # ncp
      tempNCP1 <- sqrt(n) * as.numeric(tempCov %*% coefMat[, 1]) / sqrt(apply(adjG1, 2, myvar_fun))
      tempNCP2 <- sqrt(n) * as.numeric(tempCov %*% coefMat[, 2]) / sqrt(apply(adjG2, 2, myvar_fun))
      ncpMat <- cbind(tempNCP1, tempNCP2)

      # record
      startIdx <- (set_it - 1) * setSize + 1
      endIdx <- set_it * setSize
      allZ[startIdx:endIdx, ] <- statsMat
      allNCP[startIdx:endIdx, ] <- ncpMat
      allBeta[startIdx:endIdx, ] <- coefMat
   
      # checkpoint 
      if(set_it%%100 == 0) {cat(set_it)}
    } # done generating data

    # save it 
    if (saveData) {
      write.table(allZ, paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
      write.table(allBeta, paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    }
  }
  
  # check number of non NAs
  goodIdx <- which(!is.na(allZ[, 1]) & !is.na(allZ[, 2]))
  if (length(goodIdx) < nSNPs) {
    next
  }  else {
    allZ <- allZ[goodIdx[1:nSNPs], ]
    allNCP <- allNCP[goodIdx[1:nSNPs], ]
    allBeta <- allBeta[goodIdx[1:nSNPs], ]
  }

  # number of signals and causal SNPs
  causalVec <- as.numeric(allBeta[, 1] != 0 & allBeta[, 2] != 0)
  sig2Vec <- as.numeric(abs(allNCP[, 1]) > 2 & abs(allNCP[, 2]) > 2) 

  powerRes$nCausal[sim_it] <- sum(causalVec)
  powerRes$nSig2[sim_it] <- sum(sig2Vec)
  powerRes$pi0aTrue[sim_it] <- length(which(allNCP[, 1] != 0))
  powerRes$pi0bTrue[sim_it] <- length(which(allNCP[, 2] != 0))

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
      powerRes$sig2PowHDMT[sim_it] <- length(which(totOut$hdmtRes < qvalue & totOut$sig2 == 1)) / sum(sig2Vec)
      powerRes$sig2fdpHDMT[sim_it] <- length(which(totOut$hdmtRes < qvalue & totOut$sig2 == 0)) / length(which(totOut$hdmtRes < qvalue))  
    } else {totOut <- totOut %>% mutate(hdmtRes = NA)}
  } else {totOut <- totOut %>% mutate(hdmtRes = NA)}

  # DACT
  if (class(nullprop)[1] != "list") {
    nullprop <- NULL
  }
  if (doDACT) {
    DACTout <- tryCatch(DACT_noEst(p_a = allP[, 1], p_b = allP[, 2], nullEst = nullprop, correction="JC"),
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
      powerRes$sig2PowDACT[sim_it] <- length(which(totOut$DACTrej < qvalue & totOut$sig2 == 1)) / sum(sig2Vec)
      powerRes$sig2fdpDACT[sim_it] <- length(which(totOut$DACTrej < qvalue & totOut$sig2 == 0)) / length(which(totOut$DACTrej < qvalue)) 
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
      powerRes$sig2PowKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue & totOut$sig2 == 1)) / sum(sig2Vec)
      powerRes$sig2fdpKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue & totOut$sig2 == 0)) / length(which(totOut$kernelAvg < qvalue)) 
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
      powerRes$sig2Pow7df[sim_it] <- length(which(totOut$df7Avg < qvalue & totOut$sig2 == 1)) / sum(sig2Vec)
      powerRes$sig2fdp7df[sim_it] <- length(which(totOut$df7Avg < qvalue & totOut$sig2 == 0)) / length(which(totOut$df7Avg < qvalue)) 
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
      powerRes$sig2Pow50df[sim_it] <- length(which(totOut$df50Avg < qvalue & totOut$sig2 == 1)) / sum(sig2Vec)
      powerRes$sig2fdp50df[sim_it] <- length(which(totOut$df50Avg < qvalue & totOut$sig2 == 0)) / length(which(totOut$df50Avg < qvalue))
    } else {
      totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)
    }
  } else {totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)}

  # new method
  if (doNew) {
    initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
    initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                       matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
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
    powerRes$sig2PowNew[sim_it] <- length(which(totOut$newAvg < qvalue & totOut$sig2 == 1)) / sum(sig2Vec)
    powerRes$sig2fdpNew[sim_it] <- length(which(totOut$newAvg < qvalue & totOut$sig2 == 0)) / length(which(totOut$newAvg < qvalue))
  }

  cat('\n Done with ', sim_it, '\n')
  cat('\n \n \n \n \n \n \n \n \n \n')
}

write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')


