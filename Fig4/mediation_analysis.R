# Perform mediation analysis to get raw data for Table 1

# load libraries
library(mvtnorm)
library(data.table)
library(dplyr)
library(magrittr)
library(devtools)
library(ks)
devtools::install_packages("ryanrsun/csmGmm")
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/Fig4/output"
# please place the files starting with "med_analysis_aID" (see DATA folder) in the
# output folder as well
outRoot <- paste0("med_analysis_aID", aID)
# please place this file (see DATA folder) in output folder as well
twasFname <- "scc_lung_addchr1.csv"
#-------------------------------------------------------------------#

# convergence of EM
oldEps <- 0.01
newEps <- 10^(-5)

# list of SNPs
tab2DF <- data.frame(RS = c("rs71658797", "rs6920364", "rs11780471", 
                            "rs55781567", "rs56113850", "rs13080835",
                            "rs7705526", "rs4236709", "rs885518", "rs11591710",
                            "rs1056562", "rs77468143", "rs41309931", "rs116822326",
                            "rs7953330"),
                     Gene = c("FUBP1", "RNASET2", "CHRNA2", 
                              "CHRNA5", "CYP2A6", "TP63", "TERT", "NRG1", "CDKN2A",
                              "OBFC1", "AMICA1", "SECISBP2L", "RTEL1", "HCP5", "RAD52"),
                     BP = c(77967507, 167376466, 27344719, 78857986,
                            41353107, 189357199, 1285974, 32410110, 21830157, 105687632,
                            118125625, 49376624, 62326579, 31434111, 998819),
                     Chr = c(1, 6, 8, 15, 19, 3, 5, 8, 9, 10, 11, 15, 20, 6, 12)) %>%
  slice(aID:aID)

# read the lung twas results
setwd(outputDir)
  twasRes <- read.csv(twasFname) %>%
  select(gene_name, zscore, pvalue) %>%
  set_colnames(c("Gene", "Z_twas", "p_twas")) %>%
  filter(!is.na(Z_twas))

# loop through 15 SNPs
for (snp_it in 1:nrow(tab2DF)) {
  setwd(outputDir) 
  tempFname <- paste0(tab2DF$RS[snp_it], "_", tab2DF$Gene[snp_it], "_", tab2DF$Chr[snp_it], "_", tab2DF$BP[snp_it], ".txt")
  tempEqtl <- fread(tempFname) %>% 
    select(Gene, testStat, pval) %>%
    set_colnames(c("Gene", "Z_eqtl", "p_eqtl"))

  # merge
  allDat <- tempEqtl %>% merge(twasRes, by="Gene") 
  testDat <- allDat %>%
    select(Z_twas, Z_eqtl, p_twas, p_eqtl) %>%
    as.matrix(.)

  # adjustment to not get p-values of 0 needed for DACT and HDMT and our own methods
  for (col_it in 1:(ncol(testDat)/2)) {
    tooBig <- which(testDat[, col_it] > 8.1)
    tooSmall <- which(testDat[, col_it] < -8.1)
    if (length(tooBig) > 0) {
      testDat[tooBig, col_it] <- 8.1
    }
    if (length(tooSmall) > 0) {
      testDat[tooSmall, col_it] <- -8.1
    }
  }
  for (col_it in ((ncol(testDat)/2)+1):ncol(testDat)) {
    tooBig <- which(testDat[, col_it] == 1)
    tooSmall <- which(testDat[, col_it] == 0)
    minVal <- min(testDat[which(testDat[, col_it] > 0), col_it])
    if (length(tooBig) > 0) {
      testDat[tooBig, col_it] <- 0.999
    }
    if (length(tooSmall) > 0) {
      testDat[tooSmall, col_it] <- minVal
    }
  }

  # save the test statistics
  setwd(outputDir)
  write.table(allDat, paste0(outRoot, "_dat.txt"), append=F, quote=F, row.names=F, col.names=T)

  # run HDMT 
  nullprop <- tryCatch(null_estimation(testDat[, 3:4]), error=function(e) e, warning=function(w) w)
  if (class(nullprop)[1] == "list") {
    hdmtOut <- tryCatch(fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,testDat[,3:4]),
                        error = function(e) e, warning = function(w) w)
  } else {hdmtOut <- c()}
  if (class(nullprop)[1] != "list" | class(hdmtOut)[1] != "data.frame") {
    hdmtOut <- rep(NA, nrow(testDat))
  }
  # save
  setwd(outputDir)
  write.table(hdmtOut, paste0(outRoot, "_hdmt.txt"), append=F, quote=F, row.names=F, col.names=T)

  # run the DACT
  if (class(nullprop)[1] != "list") {
    nullprop <- NULL
  }
  DACTout <- tryCatch(DACT_noEst(p_a = testDat[, 3], p_b = testDat[, 4], nullEst=nullprop, correction="JC"),
                      error = function(e) e, warning=function(w) w)
  if (class(DACTout)[1] != "list") {
    DACTfreqp <- rep(NA, nrow(testDat))
  } else {
    DACTfreqp <- DACTout$pval
  }
  # save
  setwd(outputDir)
  write.table(DACTfreqp, paste0(outRoot, "_DACTp.txt"), append=F, quote=F, row.names=F, col.names=T)

  # kernel
  oldResKernel <- emp_bayes_framework(summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = FALSE, kernel = TRUE, joint=FALSE, ind = TRUE,
                                      dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE) 
  if (class(oldResKernel)[1] != "list") {
      kernelLfdr <- rep(NA, nrow(testDat))
  } else {
      kernelLfdr <- oldResKernel$lfdrVec
  }
  # save
  setwd(outputDir)
  write.table(kernelLfdr, paste0(outRoot, "_kernel.txt"), append=F, quote=F, row.names=F, col.names=T)

  # 7 df
  oldRes7df <- emp_bayes_framework(summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = FALSE, kernel = FALSE, joint=FALSE, ind = TRUE,
                                   dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldRes7df)[1] != "list") {
        df7Lfdr <- rep(NA, nrow(testDat))
  } else {
        df7Lfdr <- oldRes7df$lfdrVec
  }
  # save
  setwd(outputDir)
  write.table(df7Lfdr, paste0(outRoot, "_df7.txt"), append=F, quote=F, row.names=F, col.names=T)

  # 50 df
  oldRes50df <- emp_bayes_framework(summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = FALSE, kernel = FALSE, joint=FALSE, ind = TRUE,
                                    dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)

  if (class(oldRes50df)[1] != "list") {
    df50Lfdr <- rep(NA, nrow(testDat))
  } else {
    df50Lfdr <- oldRes50df$lfdrVec
  }
  # save
  setwd(outputDir)
  write.table(df50Lfdr, paste0(outRoot, "df50.txt"), append=F, quote=F, row.names=F, col.names=T)

  # new
  initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
  initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                     matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
  newRes <- symm_fit_ind_EM(testStats = testDat[, 1:2], sameDirAlt = FALSE, initMuList = initMuList, initPiList = initPiList, eps=newEps)

  # save
  setwd(outputDir)
  newOut <- allDat %>% mutate(newLfdr = newRes$lfdrResults)
  write.table(newOut, paste0(outRoot, "_newlfdr.txt"), append=F, quote=F, row.names=F, col.names=T)
  write.table(do.call(cbind, newRes$muInfo), paste0(outRoot, "_muInfo.txt"), append=F, quote=F, row.names=F, col.names=T)
  write.table(do.call(cbind, newRes$piInfo), paste0(outRoot, "_piInfo.txt"), append=F, quote=F, row.names=F, col.names=T)
}


