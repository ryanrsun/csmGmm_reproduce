# mediation style, M=\betaG + eps and Y = \alphaM + eps
# n=1k subjects, J=100k SNPs
# just independent SNPs and mediators, random effect sizes.

# only difference is same effect both dimensions, instead of 0.18, 0.28

# load libraries
library(mvtnorm)
library(data.table)
library(bindata)
library(dplyr)
setwd('/rsrch3/home/biostatistics/rsun3/github/ancillaryFunctions')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)
setwd('/rsrch3/home/biostatistics/rsun3/github/altSymmMix/R')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

# input
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# output
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/Fig3/output"
outName <- paste0("sim_n1k_j100k_bincor_changeeff_raiseAlt_one_aID", aID, ".txt")

# parameters
outcomeCor <- 0.1
nDims <- 2
nSNPs <- 10^5 
setSize <- 100
n <- 1000
nSets <- nSNPs / setSize
nSims <- 5
margprob <- rep(0.3, setSize)
beta0 <- -1

# random beta
randomBeta <- TRUE
simsPerEffSize <- 40
effSizeMult <- ceiling(aID / simsPerEffSize)
betaOpts <- cbind(c(0.08 + 0.01 * effSizeMult, 0.18 + 0.01 * effSizeMult),
                  c(0.08 + 0.01 * effSizeMult, 0.18 + 0.01 * effSizeMult))
betaMin <- betaOpts[1, ]
betaMax <- betaOpts[2, ]

# make signals matrix
# change pi00
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

# get rid of the n-1 in base R
myvar_fun <- function(x) {
  sum(x^2) / length(x)
}

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
  #cat(row_it)
}


# apply over n*2 matrix of marginal probabilities of 1
gen_one_subject <- function(x, sigmaValsMat) {
  # x should be rounded
  tempVal <- sigmaValsMat[100 * x[1], 100 * x[2]]
  tempSig <- matrix(data=c(1, tempVal, tempVal, 1), nrow=2)
  return(rmvbin(n=1, margprob=x, sigma=tempSig))
}


# function to allocate the signals.
allocate_sigs <- function(hMat, nSNPs) {

  # total number of signal locations
  hSigsDF <- hMat %>% filter(s != 0)
  totSigs <- sum(hSigsDF$number)

  # sample the locations
  tempPos <- sample(x=1:nSNPs, size=totSigs, replace=F)
 
  # will put the signal locations here
  sigLocsMat <- matrix(data=NA, nrow=length(tempPos), ncol=nDims)

  # loop through all types of signals
  counter <- 1
  for (row_it in 1:nrow(hSigsDF)) {
    tempNum <- hSigsDF$number[row_it]
    if (tempNum == 0) {next} 
    for (col_it in 1:nDims) {
      sigLocsMat[counter:(counter + tempNum - 1), col_it] <- hSigsDF[row_it, col_it] * tempPos[counter:(counter + tempNum - 1)]
    }
    # increment counter
    counter <- counter + tempNum 
  }

  return(sigLocsMat)
}

# place signals
set_beta <- function(sigLocsMat, set_it, setSize, randomBeta, betaOpts=NULL, betaMin=NULL, betaMax=NULL) {
 
  # return this 
  betaMat <- matrix(data=0, nrow=setSize, ncol=ncol(sigLocsMat))  
 
  # indices of this set 
  tempIdx <- ((set_it - 1) * setSize + 1):(set_it * setSize)
 
  # loop through dimensions to add signals 
  for (col_it in 1:ncol(sigLocsMat)) {
    whichPos <- which(tempIdx %in% sigLocsMat[, col_it])
    whichNeg <- which((-1 * tempIdx) %in% sigLocsMat[, col_it])
    if (length(whichPos) > 0) {
      betaMat[whichPos, col_it] <- rep(1, length(whichPos))
    }
    if (length(whichNeg) > 0) {
      betaMat[whichNeg, col_it] <- rep(-1, length(whichNeg))
    }
  
    # random beta
    if (randomBeta) {
      betaMat[, col_it] <- betaMat[, col_it] * runif(n = setSize, min=betaMin[col_it], max=betaMax[col_it])
    } else {
      betaMat[, col_it] <- betaMat[, col_it] * betaOpts[sample(1:nrow(betaOpts), size=setSize, replace=TRUE), col_it]
    }

  }

  # return
  return(betaMat) 
}

# loop through simulations - takes 3-5 min to generate a set of 100k snps with blocks of 500 snps
powerRes <- data.frame(nCausal=rep(NA, nSims), nSig1=NA, nSig15=NA, nSig2=NA, minEff1=betaMin[1], minEff2=betaMin[1],
                       seed=NA, pi0aTrue=NA, pi0bTrue=NA, nCase1=NA, nCase2=NA, pi0aHDMT=NA, pi0bHDMT=NA,
                       pi0aKernel=NA, pi0bKernel=NA, pi0aNew=NA, pi0bNew=NA, pi0aDf7=NA, pi0bDf7=NA, pi0aDf50=NA, pi0bDf50=NA,
                       rej1DACT=NA, rej1DACTb=NA, rej1HDMT=NA, rej1Kernel=NA, rej1df7=NA, rej1df50=NA, rej1New=NA,
                       rej2DACT=NA, rej2DACTb=NA, rej2HDMT=NA, rej2Kernel=NA, rej2df7=NA, rej2df50=NA, rej2New=NA,
                       powerDACT=NA, powerDACTb=NA, powerHDMT=NA, powerKernel=NA, power7df=NA,
                       power50df=NA, powerNew=NA, fdpDACT=NA, fdpDACTb=NA, fdpHDMT=NA, fdpKernel=NA, fdp7df=NA, fdp50df=NA, fdpNew=NA,
                       nRejDACT=NA, nRejDACTb=NA, nRejHDMT=NA, nRejKernel=NA, nRej7df=NA, nRej50df=NA, nRejNew=NA,
                       nRejHDMTorig=NA, powerHDMTorig=NA, fdpHDMTorig=NA, fwerHDMT=NA, nRejDACTorig=NA, powerDACTorig=NA, fdpDACTorig=NA, fwerDACT=NA,
                       inconKernel=NA, incon7df=NA, incon50df=NA, inconNew=NA, inconSigKernel=NA, inconSig7df=NA, inconSig50df=NA,
                       estCor=NA, estCorAll=NA, pi0aCor=NA, pi0bCor=NA, rej1Cor=NA, rej2Cor=NA, powerCor=NA, nRejCor=NA, incorCor=NA)
for (sim_it in 1:nSims) {
 
  set.seed(aID * 10^5 + sim_it)
  powerRes$seed[sim_it] <- aID * 10^5 + sim_it

  # hold test statistics and signals
  allZ <- matrix(data=NA, nrow=nSNPs, ncol=2)
  allBeta <- matrix(data=NA, nrow=nSNPs, ncol=2)

  # select signal locations
  sigLocsMat <- allocate_sigs(hMat = hMat, nSNPs = nSNPs)
  
  # generate data in "sets" - easier for both fake and real
  for (set_it in 1:nSets)  {

    # information we need to save
    statsMat <- matrix(data=NA, nrow=setSize, ncol=nDims)
   
    # generate beta matrix
    coefMat <- set_beta(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize, randomBeta = randomBeta, 
                       betaMin=betaMin, betaMax=betaMax) 
    # generate genotypes
    tempG <- sapply(X=margprob, FUN=rbinom, n=n, size=2)
    # it's not a new outcome for each G, it's a new outcome for each block of 100 G
    # for logistic regression hopefully it's not a big deal? 
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

    # loop through dimensions
    for (dim_it in 1:nDims) { 
      for (snp_it in 1:setSize) {
        tempMod <- glm(yMat[, dim_it] ~ tempG[, snp_it], family=binomial)
        statsMat[snp_it, dim_it] <- summary(tempMod)$coefficients[2, 3]
      }
    }
   
    # record
    startIdx <- (set_it - 1) * setSize + 1
    endIdx <- set_it * setSize
    
    allZ[startIdx:endIdx, ] <- statsMat
    allBeta[startIdx:endIdx, ] <- coefMat
   
    # checkpoint 
    if(set_it%%1 == 0) {cat(set_it)}
  } # done generating data

  setwd(outputDir)
  write.table(allZ, paste0("old_allZ_aID", aID, "_sim", sim_it, ".txt"))
  write.table(allBeta, paste0("old_allBeta_aID", aID, "_sim", sim_it, ".txt"))

  # number of signals and causal SNPs
  causalVec <- as.numeric(allBeta[, 1] != 0 & allBeta[, 2] != 0)
  case1Vec <- as.numeric(abs(allBeta[, 1]) > 0 & allBeta[, 2] == 0)
  case2Vec <- as.numeric(allBeta[, 1] == 0 & abs(allBeta[, 2]) > 0)
  case12Vec <- case1Vec + case2Vec

  powerRes$nCausal[sim_it] <- sum(causalVec)
  powerRes$nCase1[sim_it] <- sum(case1Vec)
  powerRes$nCase2[sim_it] <- sum(case2Vec)
  powerRes$pi0aTrue[sim_it] <- length(which(allBeta[, 1] != 0))
  powerRes$pi0bTrue[sim_it] <- length(which(allBeta[, 2] != 0))

  # adjustment to not get p-values of 0 needed for DACT and HDMT and our own methods
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
  allP <- 1- pchisq(allZ^2, df=1)
  
  # analyze it 
  # start with HDMT
  nullprop <- tryCatch(null_estimation(allP), error=function(e) e, warning=function(w) w)
  if (class(nullprop)[1] == "list") {
    hdmtOut <- tryCatch(fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,allP),
                     error = function(e) e, warning = function(w) w)
    if (class(hdmtOut)[1] == "data.frame") {
      # estimated proportions of null
      powerRes$pi0aHDMT[sim_it] <- nullprop$alpha1
      powerRes$pi0bHDMT[sim_it] <- nullprop$alpha2

      # append causal indicator vector!
      hdmtOut <- hdmtOut %>% mutate(causal = causalVec, sig1=NA, case1 = case1Vec, case2=case2Vec, sig15 = NA, sig2 = NA)
      # calculate power and fdp for HDMT
      powerRes$rej1HDMT[sim_it] <- length(which(hdmtOut$fixedFdr < 0.1 & hdmtOut$case1 == 1))
      powerRes$rej2HDMT[sim_it] <- length(which(hdmtOut$fixedFdr < 0.1 & hdmtOut$case2 == 1))
      powerRes$nRejHDMT[sim_it] <- length(which(hdmtOut$fixedFdr < 0.1))
      powerRes$fdpHDMT[sim_it] <- length(which(hdmtOut$fixedFdr < 0.1 & hdmtOut$causal == 0)) /  length(which(hdmtOut$fixedFdr < 0.1))
      powerRes$powerHDMT[sim_it] <- length(which(hdmtOut$fixedFdr < 0.1 & hdmtOut$causal == 1)) / sum(causalVec)
    } else {hdmtOut <- data.frame(X1 = allP[, 1], X2 = allP[, 2], pmax=NA, origIdx = 1:nrow(allP), avgRej=NA, fdr=NA, causal=causalVec, sig1=NA, sig15 =NA, sig2=NA)}
  } else {
    hdmtOut <- data.frame(X1 = allP[, 1], X2 = allP[, 2], pmax=NA, origIdx = 1:nrow(allP), avgRej=NA, fdr=NA, causal=causalVec, sig1=NA, case12 = case12Vec, sig15=NA, sig2=NA)
  }
  
  # DACT
  if (class(nullprop)[1] != "list") {
    nullprop <- NULL
  }
  DACTout <- tryCatch(DACT(p_a = allP[, 1], p_b = allP[, 2], nullEst=nullprop, correction="JC", fdr=TRUE),
                      error = function(e) e, warning=function(w) w)
  if (class(DACTout)[1] != "list") {
    hdmtDACT <-  hdmtOut %>% mutate(DACTp = NA, DACTlfdr=NA, DACTrej = NA)
  } else {
    DACTdf <- data.frame(DACTp = DACTout$p_dact, DACTlfdr = DACTout$lfdr, origIdx = 1:length(DACTout$p_dact)) %>%
      arrange(DACTp) %>%
      mutate(rankedIdxP = 1:nrow(.)) %>%
      mutate(km = 1:nrow(.)/ nrow(.)) %>%
      mutate(RHS = km * 0.1) 
    rejected <- which(DACTdf$DACTp <= DACTdf$RHS)
    if (length(rejected) == 0) {
      maxIdx <- 0
    } else {maxIdx <- max(rejected)}
    DACTdf <- DACTdf %>% mutate(reject = ifelse(rankedIdxP <= maxIdx, 1, 0)) %>%
      arrange(origIdx)
    # append results
    hdmtDACT <- hdmtOut %>% mutate(DACTp = DACTdf$DACTp, DACTrej = DACTdf$reject)
    # power and FDP for DACT
    powerRes$rej1DACT[sim_it] <- length(which(hdmtDACT$DACTrej == 1 & hdmtDACT$case1 == 1))
    powerRes$rej2DACT[sim_it] <- length(which(hdmtDACT$DACTrej == 1 & hdmtDACT$case2 == 1)) 
    powerRes$nRejDACT[sim_it] <- length(which(hdmtDACT$DACTrej == 1))
    powerRes$fdpDACT[sim_it] <- length(which(hdmtDACT$DACTrej == 1 & hdmtDACT$causal == 0)) / length(which(hdmtDACT$DACTrej == 1))
  
    # for DACT FDR
    DACTdf <- DACTdf %>% arrange(DACTlfdr) %>%
      mutate(cumlfdrDACT = cummean(DACTlfdr)) %>%
      mutate(DACTrejB = ifelse(cumlfdrDACT < 0.1, 1, 0)) %>%
      arrange(origIdx)
    hdmtDACT <- hdmtDACT %>% mutate(DACTlfdr = DACTdf$DACTlfdr, cumlfdrDACT = DACTdf$cumlfdrDACT, DACTrejB = DACTdf$DACTrejB)
    powerRes$rej1DACTb[sim_it] <- length(which(hdmtDACT$DACTrejB == 1 & hdmtDACT$case1 == 1))
    powerRes$rej1DACTb[sim_it] <- length(which(hdmtDACT$DACTrejB == 1 & hdmtDACT$case2 == 1)) 
    powerRes$nRejDACTb[sim_it] <- length(which(hdmtDACT$DACTrejB == 1))
    powerRes$fdpDACTb[sim_it] <- length(which(hdmtDACT$DACTrejB == 1 & hdmtDACT$causal == 0)) / length(which(hdmtDACT$DACTrejB == 1))
    powerRes$powerDACTb[sim_it] <- length(which(hdmtDACT$DACTrejB == 1 & hdmtDACT$causal == 1)) / sum(causalVec)
  }

  # old method - kernel
  oldResKernel <- emp_bayes_framework(summary_tab = allZ, kernel = TRUE, joint=FALSE, ind = TRUE, 
                                dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldResKernel)[1] == "list") {
    totOut <- hdmtDACT %>% mutate(kernelLfdr = oldResKernel$lfdrVec) %>%
      arrange(kernelLfdr) %>%
      mutate(kernelAvg = cummean(kernelLfdr)) %>% 
      # don't forget to set it back to original index for further additions!!
      arrange(origIdx)
    checkIncon <- totOut %>%
      mutate(Z1 = allZ[,1], Z2 = allZ[,2]) %>%
      filter(kernelAvg < 0.1) %>%
      select(kernelLfdr, Z1, Z2) %>% 
      as.matrix(.)
    if (nrow(checkIncon) == 0) {
      inconKernelSig <- c()
    } else {
      inconKernelSig <- check_incongruous(zMatrix = checkIncon[, 2:3], lfdr = checkIncon[,1])
    }
    powerRes$pi0aKernel[sim_it] <- oldResKernel$pi0_estim[1]
    powerRes$pi0bKernel[sim_it] <- oldResKernel$pi0_estim[2] 
    powerRes$rej1Kernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$case1 == 1))
    powerRes$rej2Kernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$case2 == 1))    
    powerRes$nRejKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1))
    powerRes$fdpKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 0)) / length(which(totOut$kernelAvg < 0.1))
    powerRes$powerKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
    powerRes$inconKernel[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec))
    powerRes$inconSigKernel[sim_it] <- length(inconKernelSig)
  } else {
    totOut <- hdmtDACT %>% mutate(kernelLfdr = NA, kernelAvg=NA)
  }
  
  # old method - 7 df
  oldRes7df <- emp_bayes_framework(summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE, 
                                      dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldRes7df)[1] == "list") {
    totOut <- totOut %>% mutate(df7Lfdr = oldRes7df$lfdrVec) %>%
      arrange(df7Lfdr) %>%
      mutate(df7Avg = cummean(df7Lfdr)) %>% 
      # don't forget to set it back to original index for further additions!!
      arrange(origIdx)
    checkIncon <- totOut %>%  
      mutate(Z1 = allZ[,1], Z2 = allZ[,2]) %>%
      filter(df7Avg < 0.1) %>%
      select(df7Lfdr, Z1, Z2) %>%
      as.matrix(.)
    if (nrow(checkIncon) == 0) {
      incondf7Sig <- c()
    } else {
      incondf7Sig <- check_incongruous(zMatrix = checkIncon[, 2:3], lfdr = checkIncon[,1])
    }
    powerRes$pi0aDf7[sim_it] <- oldRes7df$pi0_estim[1]
    powerRes$pi0bDf7[sim_it] <- oldRes7df$pi0_estim[2]
    powerRes$rej1df7[sim_it] <- length(which(totOut$df7Avg < 0.1 & totOut$case1 == 1))
    powerRes$rej2df7[sim_it] <- length(which(totOut$df7Avg < 0.1 & totOut$case2 == 1)) 
    powerRes$nRej7df[sim_it] <- length(which(totOut$df7Avg < 0.1))
    powerRes$fdp7df[sim_it] <- length(which(totOut$df7Avg < 0.1 & totOut$causal == 0)) /  length(which(totOut$df7Avg < 0.1))
    powerRes$power7df[sim_it] <- length(which(totOut$df7Avg < 0.1 & totOut$causal == 1)) / sum(causalVec)
    powerRes$incon7df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes7df$lfdrVec))
    powerRes$inconSig7df[sim_it] <- length(incondf7Sig) 
  } else {
    totOut <- totOut %>% mutate(df7Lfdr = NA, df7Avg=NA)
  }
  
  # old method - 50 df
  oldRes50df <- emp_bayes_framework(summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE, 
                                   dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldRes50df)[1] == "list") {
    totOut <- totOut %>% mutate(df50Lfdr = oldRes50df$lfdrVec) %>%
      arrange(df50Lfdr) %>%
      mutate(df50Avg = cummean(df50Lfdr)) %>% 
      # don't forget to set it back to original index for further additions!!
      arrange(origIdx)
    checkIncon <- totOut %>% 
      mutate(Z1 = allZ[,1], Z2 = allZ[,2]) %>%
      filter(df50Avg < 0.1) %>%
      select(df50Lfdr, Z1, Z2) %>% 
      as.matrix(.)
    if (nrow(checkIncon) == 0) {
      incondf50Sig <- c()
    } else {
      incondf50Sig <- check_incongruous(zMatrix = checkIncon[, 2:3], lfdr = checkIncon[,1])
    }
    powerRes$pi0aDf50[sim_it] <- oldRes50df$pi0_estim[1]
    powerRes$pi0bDf50[sim_it] <- oldRes50df$pi0_estim[2] 
    powerRes$rej1df50[sim_it] <- length(which(totOut$df50Avg < 0.1 & totOut$case1 == 1))
    powerRes$rej2df50[sim_it] <- length(which(totOut$df50Avg < 0.1 & totOut$case2 == 1))
    powerRes$nRej50df[sim_it] <- length(which(totOut$df50Avg < 0.1))
    powerRes$fdp50df[sim_it] <- length(which(totOut$df50Avg < 0.1 & totOut$causal == 0)) /  length(which(totOut$df50Avg < 0.1))
    powerRes$power50df[sim_it] <- length(which(totOut$df50Avg < 0.1 & totOut$causal == 1)) / sum(causalVec)
    powerRes$incon50df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes50df$lfdrVec))
    powerRes$inconSig50df[sim_it] <- length(incondf50Sig) 
  } else {
    totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)
  }
  
  # new method
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
  checkIncon <- totOut %>%
    mutate(Z1 = allZ[,1], Z2 = allZ[,2]) %>%
    filter(newAvg < 0.1) %>%
    select(newLfdr, Z1, Z2) %>%
    as.matrix(.)
  powerRes$pi0aNew[sim_it] <- sum(newRes$piInfo[[2]]) + newRes$piInfo[[1]]
  powerRes$pi0bNew[sim_it] <- sum(newRes$piInfo[[3]]) + newRes$piInfo[[1]]  
  powerRes$rej1New[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$case1 == 1)) 
  powerRes$rej2New[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$case2 == 1))
  powerRes$nRejNew[sim_it] <- length(which(totOut$newAvg < 0.1))
  powerRes$fdpNew[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$newAvg < 0.1))
  powerRes$powerNew[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
  powerRes$inconNew[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = newRes$lfdrResults))

  # new method accounting for correlation
  estCorAll <- cor(allZ)
  tempZ <- allZ[which(abs(allZ[, 1]) < 2 & abs(allZ[, 2]) < 2), ]
  estCor <- cor(tempZ)
  initPiListCor <- list(c(0.82), c(0.04), c(0.04), c(0.1))
  initMuListCor <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3), nrow=2),
                     matrix(data=c(3, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
  corRes <- symm_fit_cor_EM(testStats = allZ, corMat = estCorAll, initMuList = initMuListCor, initPiList = initPiListCor, eps=10^(-5))
  powerRes$estCor[sim_it] <- estCor[1, 2]
  powerRes$estCorAll[sim_it] <- estCorAll[1, 2] 
  # record
  totOut <- totOut %>% mutate(corLfdr = corRes$lfdrResults) %>%
    arrange(corLfdr) %>%
    mutate(corAvg = cummean(corLfdr)) %>%
    # don't forget to set it back to original index for further additions!!
    arrange(origIdx)
  powerRes$pi0aCor[sim_it] <- sum(corRes$piInfo[[2]]) + corRes$piInfo[[1]]
  powerRes$pi0bCor[sim_it] <- sum(corRes$piInfo[[3]]) + corRes$piInfo[[1]]
  powerRes$rej1Cor[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$case1 == 1))
  powerRes$rej2Cor[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$case2 == 1)) 
  powerRes$nRejCor[sim_it] <- length(which(totOut$corAvg < 0.1))
  powerRes$fdpCor[sim_it] <- length(which(totOut$corAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$corAvg < 0.1))
  powerRes$powerCor[sim_it] <- length(which(totOut$corAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
  powerRes$inconCor[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = corRes$lfdrResults))

  cat('\n Done with ', sim_it, '\n')
  cat('\n \n \n \n \n \n \n \n \n \n')
}

setwd(outputDir)
write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')


