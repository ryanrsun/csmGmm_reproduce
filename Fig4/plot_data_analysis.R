# Make Figure 4 and Tables 1 and 2
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(data.table)
library(xtable)
library(devtools)
devtools::install.packages("ryanrsun/csmGmm")
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

# 1 is CAD and BMI
# 2 is ILCCO overall and Cardiogram CAD
# 3 is ILCCO overall and UKB BMI
# 4 is replication ILCCO overall and UKB lc
# 5 is replication CAD
# 6 is three way ILCCO overall, Cardiogram CAD, UKB BMI

#------------------------------------------------------------------#
# parameters to be changed
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/Fig4/output"
origOutputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/Fig4/origOutput"
#------------------------------------------------------------------#


# for colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# plot manhattan function
plotManhattan <- function(plotRes, chrCounts, colValues, shapeValues, ylimits, legName) {
  # arrange data by chromosome
  plotRes <- plotRes %>% arrange(Chr)
  uniqueChrs <- sort(unique(plotRes$Chr))

  # add true positions
  truePos <- rep(NA, nrow(plotRes))
  counter <- 1
  for (chr_it in 1:length(uniqueChrs)) {
    tempChr <- uniqueChrs[chr_it]
    tempDat <- plotRes %>% filter(Chr == tempChr)
    truePos[counter:(counter + nrow(tempDat) - 1)] <- rep(sum(chrCounts[1:tempChr]), nrow(tempDat)) + tempDat$BP
    counter <- counter + nrow(tempDat)
  }

  # plot
  xBreaks <- cumsum(chrCounts[-1])
  xBreaksLabs <- 1:22
  xBreaksLabs[c(9, 11, 13, 15, 16, 17, 19, 20, 21)] <- ""

  plotDat <- plotRes %>% mutate(truePos = truePos)

  returnPlot <- ggplot(plotDat, aes(x=truePos, y=-log10(newLfdr), color=as.factor(cat), shape=as.factor(cat))) +
    geom_point() +
    xlab("Chromosome") + ylab("-log10(lfdr)") +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(name="Chr", breaks=xBreaks, labels=xBreaksLabs) +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=18),
          legend.title=element_text(size=18), legend.text=element_text(size=18)) +
    guides(colour = guide_legend(override.aes = list(size=4)))


  return(returnPlot)
}

# add position information to data
setwd(outputDir)
s1 <- fread("reject_bmi_with_overall_neg5_reject_aID1.txt")
s2 <- fread("reject_bmi_with_overall_neg5_reject_aID2.txt")
s3 <- fread("reject_bmi_with_overall_neg5_reject_aID3.txt")
s4 <- fread("reject_bmi_with_overall_neg5_reject_aID4.txt")
s5 <- fread("reject_bmi_with_overall_neg5_reject_aID5.txt")
s1new <- s1 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))
s2new <- s2 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))
s3new <- s3 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))
s4new <- s4 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))
s5new <- s5 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))

# for plotting axes
setwd(outputDir)
allZ <- fread("bmi_with_overall.txt") %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP))
chrCounts <- rep(0, 23)
for (chr_it in 1:22) {
  tempDat <- allZ %>% filter(Chr == chr_it)
  maxPos <- max(tempDat$BP)
  chrCounts[chr_it + 1] <- maxPos
}

# data for manhattan plot - 2 way pleiotropy
manDataTwo <- rbind(s1new %>% select(chrpos, Chr, BP, newLfdr) %>% mutate(cat = "CAD,BMI"),
                 s2new %>% select(chrpos, Chr, BP, newLfdr) %>% mutate(cat = "LC,CAD"),
                 s3new %>% select(chrpos, Chr, BP, newLfdr) %>% mutate(cat = "LC,BMI")) %>%
  as.data.frame(.) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP)) %>%
  arrange(newLfdr) %>%
  # distinct() keeps the first one, so we arrange first
  distinct(., chrpos, .keep_all = TRUE)

manPlotTwo <- plotManhattan(plotRes = manDataTwo, chrCounts,
                            colValues = gg_color_hue(3), shapeValues=c(8,16,17), ylimits=c(0, 12.5), legName="Pleiotropy")
manPlotTwo


# data for manhattan plot - replication
manDataRep <- rbind(s5new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(pheno = "CAD"),
                    s4new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(pheno = "LC"),
                    s2new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(pheno = "Both")) %>%
  as.data.frame(.) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP)) %>%
  # overlap with LC, CAD pleiotropy
  mutate(Pleio = ifelse(chrpos %in% s2new$chrpos, 1, 0)) %>%
  mutate(RepCAD = ifelse(chrpos %in% s5new$chrpos, 1, 0)) %>%
  mutate(RepLC = ifelse(chrpos %in% s4new$chrpos, 1, 0)) %>%
  mutate(cat = ifelse(Pleio == 1 & RepCAD == 1 & RepLC == 1, "Pleio + Rep", "Rep Only")) %>%
  mutate(cat = ifelse(Pleio == 1 & (RepCAD == 0 | RepLC == 0), "Pleio Only", cat)) %>%  
  arrange(pheno) %>%
  # distinct() keeps the first one, so we arrange first - want to show Pleio lfdr if in pleio
  distinct(., chrpos, .keep_all = TRUE)

manPlotRep <- plotManhattan(plotRes = manDataRep, chrCounts,
                            colValues=c(gg_color_hue(3)[3], "black", "darkorange"), shapeValues=c(17, 18, 15),
                            ylimits=c(0, 12.5), legName="Lung Cancer")
manPlotRep

#----------------------------------------------------------------------------#
# Table 2

# read summary of pleiotropy analysis
setwd(outputDir)
adjAnal <- fread("processed_ukb_data_S1.txt")
origAnal <- fread("processed_ukb_data_S2.txt")

tab2 <- origAnal %>% filter(aID == 6) %>% select(Method, numReject) %>%
  mutate(a4 = origAnal %>% filter(aID == 2) %>% select(numReject) %>% unlist(.)) %>%
  mutate(a4new = adjAnal %>% filter(aID == 2) %>% select(numReject) %>% unlist(.)) %>%
  mutate(a5 = origAnal %>% filter(aID == 3) %>% select(numReject) %>% unlist(.)) %>%
  mutate(a5new = adjAnal %>% filter(aID == 3) %>% select(numReject) %>% unlist(.)) %>%
  mutate(a3 = origAnal %>% filter(aID == 1) %>% select(numReject) %>% unlist(.)) %>%
  mutate(a3new = adjAnal %>% filter(aID == 1) %>% select(numReject) %>% unlist(.))

tab2final <- tab2 %>%
  set_colnames(c("Method", "All3", "Orig", "Adj", "Orig ", "Adj ", "Orig  ", "Adj  ")) %>%
  mutate(Method = c("csmGmm", "Kernel", "locfdr50", "locfdr7", "HDMT", "DACT"))
tab2final

#----------------------------------------------------------------------------------------#


# Table 1
setwd(outputDir)
qval <- 0.1

initDat <- fread("med_analysis_aID1_newlfdr.txt") %>%
  select(Gene, Z_eqtl, p_eqtl, Z_twas, p_twas)

# hold results
allResList <- list()
rejTab <- c()
for(snp_it in 1:15) {
  
  tempRoot <- paste0("med_analysis_aID", snp_it)
  tempStatsDat <- fread(paste0(tempRoot, "_dat.txt")) %>%
    mutate(SNP = snp_it)

  temphdmtDat <- fread(paste0(tempRoot, "_hdmt.txt")) %>% select(fdr, fixedFdr) %>% set_colnames(c("hdmtLfdr", "hdmtFixed")) %>%
      mutate(rejHDMT = ifelse(hdmtFixed < qval, 1, 0)) 

  if (file.exists(paste0(tempRoot, "_newlfdr.txt"))) { 
    tempNewDat <- fread(paste0(tempRoot, "_newlfdr.txt")) %>%
      mutate(idx = 1:nrow(.))  %>%
      arrange(newLfdr) %>%
      mutate(avgNew = cummean(newLfdr)) %>%
      mutate(rejNew = ifelse(avgNew < qval, 1, 0)) %>%
      arrange(idx) %>%
      select(-idx, -Gene, -Z_eqtl, -p_eqtl, -Z_twas, -p_twas)
  } else {tempNewDat <- data.frame(newLfdr=NA, rejNew=rep(0, nrow(tempStatsDat)), avgNew=NA)}
  
  if (file.exists(paste0(tempRoot, "_kernel.txt"))) {
    tempKernelDat <- fread(paste0(tempRoot, "_kernel.txt")) %>% set_colnames(c("kernelLfdr")) %>%
      mutate(idx = 1:nrow(.)) %>% 
      arrange(kernelLfdr) %>%
      mutate(avgKernel = cummean(kernelLfdr)) %>%
      mutate(rejKernel = ifelse(avgKernel < qval, 1, 0)) %>%
      arrange(idx) %>%
      select(-idx) 
  } else {tempKernelDat <- data.frame(kernelLfdr=NA, rejKernel=rep(0, nrow(tempStatsDat)), avgKernel=NA)}

  if (file.exists(paste0(tempRoot, "df50.txt"))) { 
    tempDf50Dat <- fread(paste0(tempRoot, "df50.txt")) %>% set_colnames(c("Df50Lfdr")) %>%
      mutate(idx = 1:nrow(.)) %>%      
      arrange(Df50Lfdr) %>%
      mutate(avgDf50 = cummean(Df50Lfdr)) %>%
      mutate(rejDf50 = ifelse(avgDf50 < qval, 1, 0)) %>%
      arrange(idx) %>%
      select(-idx) 
  } else {tempDf50Dat <- data.frame(Df50Lfdr=NA, rejDf50=rep(0, nrow(tempStatsDat)), avgDf50=NA)}

  if (file.exists(paste0(tempRoot, "_df7.txt"))) {
    tempDf7Dat <- fread(paste0(tempRoot, "_df7.txt")) %>% set_colnames(c("Df7Lfdr")) %>%
      mutate(idx = 1:nrow(.)) %>%      
      arrange(Df7Lfdr) %>%
      mutate(avgDf7 = cummean(Df7Lfdr)) %>%
      mutate(rejDf7 = ifelse(avgDf7 < qval, 1, 0)) %>%
      arrange(idx) %>%
      select(-idx)   
  } else {tempDf7Dat <- data.frame(Df7Lfdr=NA, rejDf7=rep(0, nrow(tempStatsDat)), avgDf7=NA)} 
  
  tempDACTDat <- fread(paste0(tempRoot, "_DACTp.txt")) %>% set_colnames("DACTp") %>%
    mutate(idx = 1:nrow(.)) %>%  
    arrange(DACTp) %>%
    mutate(rankedIdxP = 1:nrow(.)) %>%
    mutate(km = 1:nrow(.) / nrow(.)) %>%
    mutate(RHS = km * qval)

  # for DACT do BH correction
  rejDACT <- which(tempDACTDat$DACTp <= tempDACTDat$RHS)
  if (length(rejDACT) == 0) {
    maxIdx <- 0
  } else {maxIdx <- max(rejDACT)}
  tempDACTDat <- tempDACTDat %>% mutate(rejDACT = ifelse(rankedIdxP <= maxIdx, 1, 0)) %>%
    arrange(idx)

  # put in list
  fullDat <- cbind(tempStatsDat, temphdmtDat, tempNewDat, tempKernelDat, tempDf50Dat, tempDf7Dat, tempDACTDat) 
  allResList[[snp_it]] <- fullDat

  # rejections
  tempRej <- fullDat %>% filter(rejDACT == 1 | rejHDMT == 1 | rejDf50 == 1 | rejDf7 == 1 | rejNew == 1 | rejKernel == 1)
  rejTab <- rbind(rejTab, tempRej)

  cat(snp_it)
}

# look at rejections
allRej <- rejTab %>% mutate(numRej = rejDACT+ rejHDMT +rejDf50 + rejDf7+ rejNew+rejKernel, na.rm=TRUE) %>%
  arrange(desc(numRej)) %>%
  slice(1:5) %>%
  select(Gene, Z_eqtl, Z_twas, numRej, SNP)

# merge with information
tab2DF <- data.frame(RS = c("rs71658797", "rs6920364", "rs11780471",
                            "rs55781567", "rs56113850", "rs13080835",
                            "rs7705526", "rs4236709", "rs885518", "rs11591710",
                            "rs1056562", "rs77468143", "rs41309931", "rs116822326",
                            "rs7953330"),
                     Locus = c("FUBP1", "RNASET2", "CHRNA2",
                              "CHRNA5", "CYP2A6", "TP63", "TERT", "NRG1", "CDKN2A",
                              "OBFC1", "AMICA1", "SECISBP2L", "RTEL1", "HCP5", "RAD52"),
                     BP = c(77967507, 167376466, 27344719, 78857986,
                            41353107, 189357199, 1285974, 32410110, 21830157, 105687632,
                            118125625, 49376624, 62326579, 31434111, 998819),
                     Chr = c(1, 6, 8, 15, 19, 3, 5, 8, 9, 10, 11, 15, 20, 6, 12), SNP=1:15)

mergedRej <- merge(allRej, tab2DF, by="SNP") %>%
  select(RS, Chr, BP, Locus, Gene, Z_eqtl, Z_twas, numRej)



