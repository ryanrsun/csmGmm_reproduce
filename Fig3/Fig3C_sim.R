# Figure 3C

# load libraries
library(mvtnorm)
library(data.table)
library(bindata)
library(dplyr)
library(rje)
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/Fig3/output"
outName <- paste0("Fig3C_aID", aID, ".txt")

# option to save or load intermediate data to save time,
# set as FALSE for first run and then TRUE thereafter
loadData <- FALSE
saveData <- TRUE
testStatsName <- "Fig3C_allZ"
betaName <- "Fig3C_allBeta"
#-------------------------------------------------------------------#

# parameters
outcomeCor <- 0.1
doHDMT <- TRUE
doDACT <- TRUE
doKernel <- TRUE
do50df <- TRUE
do7df <- TRUE
doNew <- TRUE
qvalue <- 0.1
nSNPs <- 10^5
setSize <- 100
n <- 1000
nDims <- 2
nSets <- nSNPs / setSize
nSims <- 5
margprob <- rep(0.3, setSize)
simsPerEffSize <- 40
effSizeMult <- ceiling(aID / simsPerEffSize)
betaOpts <- cbind(c(0.08 + 0.01 * effSizeMult, 0.18 + 0.01 * effSizeMult),
                  c(0.08 + 0.01 * effSizeMult, 0.18 + 0.01 * effSizeMult))
betaMin <- betaOpts[1, ]
betaMax <- betaOpts[2, ]
beta0 <- -1

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

# tells rmvbin how to generate correlated binary outcomes
sigmaValsMat <- matrix(data=0, nrow=99, ncol=99)
for (row_it in 1:99) {
  for (col_it in row_it:99) {
    mp1 <- row_it / 100
    mp2 <- col_it / 100
    tempCor <- min(c(mp1, mp2, outcomeCor))
    tempBC <- matrix(data=c(1, tempCor, tempCor, 1), nrow=2)
    tempCP <- bincorr2commonprob(margprob=c(mp1, mp2), bincorr = tempBC)
    tempSig <- tryCatch(commonprob2sigma(commonprob = tempCP), error=function(e) e)
    if (class(tempSig)[1] == "simpleError") {next}
    sigmaValsMat[row_it, col_it] <- tempSig[1, 2]
    sigmaValsMat[col_it, row_it] <- tempSig[1, 2]
  }
}

# generate two correlated outcomes for one subject
gen_one_subject <- function(x, sigmaValsMat) {
  # x should be rounded
  tempVal <- sigmaValsMat[100 * x[1], 100 * x[2]]
  tempSig <- matrix(data=c(1, tempVal, tempVal, 1), nrow=2)
  return(rmvbin(n=1, margprob=x, sigma=tempSig))
}

# record results here
powerRes <- data.frame(nCausal=rep(NA, nSims), minEff1=betaMin[1],
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
    setwd(outputDir)
    allZ <- fread(paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
    allBeta <- fread(paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
  } else { 
    # hold test statistics and signals
    allZ <- matrix(data=NA, nrow=nSNPs, ncol=2)
    allBeta <- matrix(data=NA, nrow=nSNPs, ncol=2)

    # select signal locations
    sigLocsMat <- allocate_sigs(hMat = hMat, nSNPs = nSNPs, nDims=nDims)

    # generate data in multiple sets - faster than all at once
    for (set_it in 1:nSets)  {

      # generate coefficient matrix
      statsMat <- matrix(data=NA, nrow=setSize, ncol=nDims) 
      coefMat <- set_beta(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize,
                       betaMin=betaMin, betaMax=betaMax)
      
      # generate genotypes
      tempG <- sapply(X=margprob, FUN=rbinom, n=n, size=2)
      etaMat <- tempG %*% coefMat + matrix(data=beta0, nrow=nrow(tempG), ncol=ncol(coefMat))
      muMat <- rje::expit(etaMat)
      # make the outcome
      roundedMu <- round(muMat, digits=2)
      for (dim_it in 1:nDims) {
        tooSmall <- which(roundedMu[, dim_it] == 0)
        if (length(tooSmall) > 0) {
          roundedMu[tooSmall, dim_it] <- 0.01
        }
        tooBig <- which(roundedMu[, dim_it] == 1)
        if (length(tooBig) > 0) {
          roundedMu[tooBig, dim_it] <- 0.99
       }
      }
      yMat <- t(apply(roundedMu, 1, gen_one_subject, sigmaValsMat = sigmaValsMat))

      # calculate test statistics
      for (dim_it in 1:nDims) {
        for (snp_it in 1:setSize) {
          tempMod <- glm(yMat[, dim_it] ~ tempG[, snp_it], family=binomial)
          statsMat[snp_it, dim_it] <- summary(tempMod)$coefficients[2, 3]
        }
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
      powerRes$fdpHDMT[sim_it] <- length(which(totOut$hdmtRes < 0.1 & totOut$causal == 0)) /  length(which(totOut$hdmtRes < 0.1))
      powerRes$powerHDMT[sim_it] <- length(which(totOut$hdmtRes < 0.1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$nRejHDMT[sim_it] <- length(which(totOut$hdmtRes < 0.1))

    } else {totOut <- totOut %>% mutate(hdmtRes = NA)}
  } else {totOut <- totOut %>% mutate(hdmtRes = NA)}

  # DACT
  if (class(nullprop)[1] != "list") {
    nullprop <- NULL
  }
  if (doDACT) {
    DACTout <- tryCatch(DACT_noEst(p_a = allP[, 1], p_b = allP[, 2], nullEst = nullprop, correction="JC"),
                        error = function(e) e, warning=function(w) w)
    if (class(DACTout)[1] %in% c("simpleWarning", "simpleError")) {
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
      powerRes$fdpKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 0)) / length(which(totOut$kernelAvg < 0.1))
      powerRes$powerKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$inconKernel[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec))
      powerRes$nRejKernel[sim_it] <-  length(which(totOut$kernelAvg < 0.1)) 
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
  
  # new method accounting for correlation
  estCorAll <- cor(allZ)
  initPiListCor <- list(c(0.82), c(0.04), c(0.04), c(0.1))
  initMuListCor <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3), nrow=2),
                     matrix(data=c(3, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
  corRes <- symm_fit_cor_EM(testStats = as.matrix(allZ), corMat = estCorAll, initMuList = initMuListCor, initPiList = initPiListCor, eps=10^(-5))
  powerRes$estCorAll[sim_it] <- estCorAll[1, 2]
  # record
  totOut <- totOut %>% mutate(corLfdr = corRes$lfdrResults) %>%
    arrange(corLfdr) %>%
    mutate(corAvg = cummean(corLfdr)) %>%
    # don't forget to set it back to original index for further additions!!
    arrange(origIdx)
  powerRes$nRejNew[sim_it] <- length(which(totOut$corAvg < 0.1))
  powerRes$fdpNew[sim_it] <- length(which(totOut$corAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$corAvg < 0.1))
  powerRes$powerNew[sim_it] <- length(which(totOut$corAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)

  cat('\n Done with ', sim_it, '\n')
}

setwd(outputDir)
write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')


