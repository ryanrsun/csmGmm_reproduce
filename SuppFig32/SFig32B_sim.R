# For Supp Figure 32B

# load libraries
library(mvtnorm)
library(magrittr)
library(data.table)
library(bindata)
library(dplyr)
library(ks)
devtools::install.packages("ryanrsun/csmGmm")
library(csmGmm)
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig32/output"
outName <- paste0("SFig32B_aID", aID, ".txt")

# option to save or load intermediate data to save time,
# set as FALSE for first run and then TRUE thereafter
loadData <- FALSE
saveData <- FALSE
testStatsName <- "SFig32B_allZ"
betaName <- "SFig32B_allBeta"
#-------------------------------------------------------------------#

# parameters
doHDMT <- TRUE
doDACT <- TRUE
doKernel <- TRUE
do50df <- TRUE
do7df <- TRUE
doNew <- TRUE
qvalue <- 0.1
nSNPs <- 10^5 
rho <- 0.4
corMat <- matrix(data=c(1, rho, rho, 1), nrow=2)
nSims <- 5
simsPerEffSize <- 20
effSizeMult <- ceiling(aID / simsPerEffSize)
mu1 <- 3.5
mu2 <- 3.5
mu3 <- 3.5

# make signals matrix
pi1 <- 0.01 * effSizeMult
pi11 <- 3 * pi1^2
pi00 <- 1 - 2*pi1 - pi11
sProp <- c(pi00, 2 * pi1, pi11)
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

# function to generate the test statistics
genData <- function(hMat, nSNPs, corMat, mu1, mu2, mu3) {

  n1 <- sum(hMat$number[which(abs(hMat$Var1) == 0 & abs(hMat$Var2) == 1)]) / 2
  n2 <- sum(hMat$number[which(abs(hMat$Var1) == 1 & abs(hMat$Var2) == 0)]) / 2
  n3 <- round(sum(hMat$number[which(abs(hMat$Var1) == 1 & abs(hMat$Var2) == 1)]) / 4) 
  numSigs <- 2*n1 + 2*n2 + 4*n3
  sigLocs <- sample(x=1:nSNPs, size=numSigs, replace=FALSE)
  sigMat <- matrix(data=0, nrow=nSNPs, ncol=2)
  sigMat[sigLocs[1:n1], 2] <- runif(n=n1, min=mu1 - 0.2, max=mu1 + 0.2)
  sigMat[sigLocs[(n1+1):(2*n1)], 2] <- -runif(n=n1, min=mu1 - 0.2, max=mu1 + 0.2)
  sigMat[sigLocs[(2*n1+1):(2*n1+n2)], 1] <- runif(n=n2, min=mu2 - 0.2, max=mu2 + 0.2)
  sigMat[sigLocs[(2*n1+n2+1):(2*n1+2*n2)], 1] <- -runif(n=n2, min=mu2 - 0.2, max=mu2 + 0.2)
  sigMat[sigLocs[(2*n1+2*n2+1):(2*n1+2*n1+n3)], 1] <- runif(n=n3, min=mu3 - 0.2, max=mu3 + 0.2)
  sigMat[sigLocs[(2*n1+2*n2+1):(2*n1+2*n1+n3)], 2] <- runif(n=n3, min=mu3 - 0.2, max=mu3 + 0.2)
  sigMat[sigLocs[(2*n1+2*n1+n3+1):(2*n1+2*n1+2*n3)], 1] <- runif(n=n3, min=mu3 - 0.2, max=mu3 + 0.2)
  sigMat[sigLocs[(2*n1+2*n1+n3+1):(2*n1+2*n1+2*n3)], 2] <- -runif(n=n3, min=mu3 - 0.2, max=mu3 + 0.2)
  sigMat[sigLocs[(2*n1+2*n1+2*n3+1):(2*n1+2*n1+3*n3)], 1] <- -runif(n=n3, min=mu3 - 0.2, max=mu3 + 0.2)
  sigMat[sigLocs[(2*n1+2*n1+2*n3+1):(2*n1+2*n1+3*n3)], 2] <- runif(n=n3, min=mu3 - 0.2, max=mu3 + 0.2)
  sigMat[sigLocs[(2*n1+2*n1+3*n3+1):(2*n1+2*n1+4*n3)], 1] <- -runif(n=n3, min=mu3 - 0.2, max=mu3 + 0.2)
  sigMat[sigLocs[(2*n1+2*n1+3*n3+1):(2*n1+2*n1+4*n3)], 2] <- -runif(n=n3, min=mu3 - 0.2, max=mu3 + 0.2)

  # first set of signals
  sigK1 <- matrix(data=sigMat[, 1], ncol=2, byrow=T)
  sigK2 <- matrix(data=sigMat[, 2], ncol=2, byrow=T)

  # function to generate Z
  genFullLikZ <- function(x, corMat) {
      rmvnorm(n=1, mean=x, sigma=corMat)
  }

  # generate Z
  Zcol1 <- t(apply(sigK1, 1, genFullLikZ, corMat=corMat))
  Zcol2 <- t(apply(sigK2, 1, genFullLikZ, corMat=corMat))

  # put back in original form
  testStats <- cbind(as.numeric(t(Zcol1)), as.numeric(t(Zcol2)))

  return(list(testStats = testStats, sigMat = sigMat))
}

# record results here
powerRes <- data.frame(nCausal=rep(NA, nSims),  minEff1=pi1,
                       seed=NA, pi0aTrue=NA, pi0bTrue=NA,
                       powerDACT=NA, powerHDMT=NA, powerKernel=NA, power7df=NA,
                       power50df=NA, powerNew=NA, fdpDACT=NA, fdpHDMT=NA, fdpKernel=NA, fdp7df=NA, fdp50df=NA, fdpNew=NA,
                       inconKernel=NA, incon7df=NA, incon50df=NA, inconNew=NA)
# each loop is one simulation iteration
for (sim_it in 1:nSims) {

  # set the seed 
  set.seed(aID * 10^5 + sim_it)
  powerRes$seed[sim_it] <- aID * 10^5 + sim_it

	# hold test statistics and signals
  dat <- genData(hMat = hMat, nSNPs = nSNPs, corMat = corMat, mu1 = mu1, mu2 = mu2, mu3 = mu3)
  allZ <- dat$testStats
  allBeta <- dat$sigMat

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
  }
  
  # full likelihood
  tempDat <- rbind(matrix(data=allZ[, 1], ncol=2, byrow=TRUE),
                  matrix(data=allZ[, 2], ncol=2, byrow=TRUE))
  estCor <- cor(tempDat)[1, 2]
  estCorMat <- matrix(data=c(1, estCor, estCor, 1), nrow=2)
  initPiListFull <- list(c(0.82), c(0.04), c(0.04), c(0.1))
  initMuListFull <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3), nrow=2),
                     matrix(data=c(3, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
  fullRes <- symm_fit_cor_EM_fulllik(testStats = allZ, corMat = estCorMat, initMuList = initMuListFull, initPiList = initPiListFull, eps=10^(-5))
  # record
  totOut <- totOut %>% mutate(fullLfdr = fullRes$lfdrResults) %>%
    arrange(fullLfdr) %>%
    mutate(fullAvg = cummean(fullLfdr)) %>%
    # don't forget to set it back to original index for further additions!!
    arrange(origIdx)
  powerRes$nRejDACT[sim_it] <- length(which(totOut$fullAvg < 0.1))
  powerRes$fdpDACT[sim_it] <- length(which(totOut$fullAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$fullAvg < 0.1))
  powerRes$powerDACT[sim_it] <- length(which(totOut$fullAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
  powerRes$inconDACT[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = fullRes$lfdrResults)) 
  
  cat('\n Done with ', sim_it, '\n')
}

setwd(outputDir)
write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')


