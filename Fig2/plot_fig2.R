# Collect results and plot Figure 2

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/Fig1/plot_fig2.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Fig2/plot_fig2.R")

library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(here)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)


# get output file names
outputDir <- here::here("Fig2", "output")
names2aq1 <- here::here(outputDir, paste0("Power_correction2A_S1_aID", 1:400, ".txt"))
names2aq2 <- here::here(outputDir, paste0("Power_correction2A_S2_aID", 1:400, ".txt"))
names2aq3 <- here::here(outputDir, paste0("Power_correction2A_S3_aID", 1:400, ".txt"))
names2aq4 <- here::here(outputDir, paste0("Power_correction2A_S4_aID", 1:400, ".txt"))
names2aq5 <- here::here(outputDir, paste0("Power_correction2A_S5_aID", 1:400, ".txt"))
names2aq6 <- here::here(outputDir, paste0("Power_correction2A_S6_aID", 1:400, ".txt"))
names2aq7 <- here::here(outputDir, paste0("Power_correction2A_S7_aID", 1:400, ".txt"))
names2aq8 <- here::here(outputDir, paste0("Power_correction2A_S8_aID", 1:400, ".txt"))
names2aq9 <- here::here(outputDir, paste0("Power_correction2A_S9_aID", 1:400, ".txt"))
names2aList <- list(names2aq1, names2aq2, names2aq3, names2aq4, names2aq5, names2aq6,
                    names2aq7, names2aq8, names2aq9)
names2bq1 <- here::here(outputDir, paste0("Power_correction2B_S1_aID", 1:400, ".txt"))
names2bq2 <- here::here(outputDir, paste0("Power_correction2B_S2_aID", 1:400, ".txt"))
names2bq3 <- here::here(outputDir, paste0("Power_correction2B_S3_aID", 1:400, ".txt"))
names2bq4 <- here::here(outputDir, paste0("Power_correction2B_S4_aID", 1:400, ".txt"))
names2bq5 <- here::here(outputDir, paste0("Power_correction2B_S5_aID", 1:400, ".txt"))
names2bq6 <- here::here(outputDir, paste0("Power_correction2B_S6_aID", 1:400, ".txt"))
names2bq7 <- here::here(outputDir, paste0("Power_correction2B_S7_aID", 1:400, ".txt"))
names2bq8 <- here::here(outputDir, paste0("Power_correction2B_S8_aID", 1:400, ".txt"))
names2bq9 <- here::here(outputDir, paste0("Power_correction2B_S9_aID", 1:400, ".txt"))
names2bList <- list(names2bq1, names2bq2, names2bq3, names2bq4, names2bq5, names2bq6,
                    names2bq7, names2bq8, names2bq9)

# read raw output files
res2a <- c()
for (file_it in 1:length(names2aq1)) {
  for (q_it in 1:9) {
    tempRes <- tryCatch(fread(names2aList[[q_it]][file_it]), error=function(e) e)
    if (class(tempRes)[1] == "simpleError") {next}
    tempList <- list(res2a, tempRes %>% mutate(q = q_it))
    res2a <- rbindlist(tempList)
  }
}

res2b <- c()
for (file_it in 1:length(names2bq1)) {
  for (q_it in 1:9) {
    tempRes <- tryCatch(fread(names2bList[[q_it]][file_it]), error=function(e) e)
    if (class(tempRes)[1] == "simpleError") {next}
    tempList <- list(res2b, tempRes %>% mutate(q = q_it))
    res2b <- rbindlist(tempList)
  }
}

# summarize
setwd(outputDir)
for (q_it in 1:9) {
  tempResa <- res2a %>% filter(q == q_it)
  tempSummarya <- summarize_raw(tempResa)
  tempSummaryaName <- here::here(outputDir, paste0("Fig2a_power", q_it, "_summary.txt"))
  write.table(tempSummarya, append=F, quote=F, row.names=F, col.names=T, sep='\t')

  tempResb <- res2b %>% filter(q == q_it)
  tempSummaryb <- summarize_raw(tempResb)
  tempSummarybName <- here::here(outputDir, paste0("Fig2b_power", q_it, "_summary.txt"))
  write.table(tempSummaryb, tempSummarybName, append=F, quote=F, row.names=F, col.names=T, sep='\t')
}

#------------------------------------------------#
# plotting starts

# define colors
gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"


# load all the files needed for corrected power
tempSummaryaName <- here::here(outputDir, paste0("Fig2a_power1_summary.txt"))
Fig2a_fixq <- fread(tempSummaryaName, data.table=F)  %>% mutate(q=0.01)
for (q_it in 2:9) {
  tempSummaryaName <- here::here(outputDir, paste0("Fig2a_power", q_it, "_summary.txt"))
  tempFile <- fread(tempSummaryaName, data.table=F) %>% mutate(q = q_it / 100)
  Fig2a_fixq <- rbind(Fig2a_fixq, tempFile)
}
Fig2a_fixq <- Fig2a_fixq %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  rbind(., fread("Fig1c_summary.txt") %>% mutate(q = 0.1) %>% 
        mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
        mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
        mutate(Method = ifelse(Method == "df50", "locfdr50df", Method))) %>%
  filter(minEff1 < 0.2)

# correct the power
Fig2a_correctedPow <- Fig2a_fixq %>% filter(q == 0.1) %>%
  select(minEff1, Method, Power) %>%
  filter(minEff1 < 0.2) %>%
  mutate(qval = NA, actualPower=NA, actualfdp=NA)
allEffs <- unique(Fig2a_fixq$minEff1)
allMets <- unique(Fig2a_fixq$Method)
for (eff_it in 1:length(allEffs)) {
  tempEff <- allEffs[eff_it]
  for (met_it in 1:length(allMets)) {
    tempMet <- allMets[met_it]
    if (tempMet %in% c("HDMT", "csmGmm")) {next}
    tempDat <- Fig2a_fixq %>% filter(minEff1 == tempEff & Method == tempMet) %>%
      arrange(q)
    minRow <- max(which(tempDat$FDP < 0.1))
    fillRow <- which(Fig2a_correctedPow$minEff1 == tempEff & Fig2a_correctedPow$Method == tempMet)
    Fig2a_correctedPow$qval[fillRow] <- tempDat$q[minRow]
    Fig2a_correctedPow$actualPower[fillRow] <- tempDat$Power[minRow]
    Fig2a_correctedPow$actualfdp[fillRow] <- tempDat$FDP[minRow]
  }
}
Fig2a_correctedPow <- Fig2a_correctedPow %>%
  filter(minEff1 <= 0.19) %>%
  mutate(actualPower = ifelse(Method %in% c("HDMT", "csmGmm"), Power, actualPower)) 

# plot
Fig2a_plot <- ggplot(data=Fig2a_correctedPow,
             aes(x=minEff1, y=actualPower, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") + ylab("Power (2D Mediation)") +
  ylim(c(0, 0.8)) + xlim(c(0, 0.2)) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# Fig 2b
tempSummarybName <- here::here(outputDir, paste0("Fig2b_power1_summary.txt"))
Fig2b_fixq <- fread(tempSummarybName, data.table=F)  %>% mutate(q = 0.01)
for (q_it in 2:9) {
  tempSummarybName <- here::here(outputDir, paste0("Fig2b_power", q_it, "_summary.txt"))
  tempFile <- fread(tempSummarybName, data.table=F) %>% mutate(q = q_it / 100)
  Fig2b_fixq <- rbind(Fig2b_fixq, tempFile)
}
Fig2b_fixq <- Fig2b_fixq %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  rbind(., fread("Fig1d_summary.txt") %>% mutate(q = 0.1) %>%
        mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
        mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
        mutate(Method = ifelse(Method == "df50", "locfdr50df", Method))) %>%
  filter(minEff1 < 0.08)
 
# correct the power
Fig2b_correctedPow <- Fig2b_fixq %>% filter(q == 0.1) %>%
  select(minEff1, Method, Power) %>%
  mutate(qval = NA, actualPower=NA, actualfdp=NA)
allEffs <- unique(Fig2b_fixq$minEff1)
allMets <- unique(Fig2b_fixq$Method)
for (eff_it in 1:length(allEffs)) {
  tempEff <- allEffs[eff_it]
  for (met_it in 1:length(allMets)) {
    tempMet <- allMets[met_it]
    if (tempMet %in% c("HDMT", "csmGmm")) {next}
    tempDat <- Fig2b_fixq %>% filter(minEff1 == tempEff & Method == tempMet) %>%
      arrange(q)
    minRow <- which.min(abs(tempDat$FDP - 0.1))
    fillRow <- which(Fig2b_correctedPow$minEff1 == tempEff & Fig2b_correctedPow$Method == tempMet)
    Fig2b_correctedPow$qval[fillRow] <- tempDat$q[minRow]
    Fig2b_correctedPow$actualPower[fillRow] <- tempDat$Power[minRow]
    Fig2b_correctedPow$actualfdp[fillRow] <- tempDat$FDP[minRow]
  }
}
# restrict
Fig2b_correctedPow <- Fig2b_correctedPow %>%
  mutate(actualPower = ifelse(Method %in% c("HDMT", "csmGmm"), Power, actualPower))

# plot
Fig2b_plot <- ggplot(data=Fig2b_correctedPow,
             aes(x=minEff1, y=actualPower, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("Power (2D Mediation)") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

# Fig 2C
summary1c <- fread(here::here("Fig1", "output", "Fig1c_summary.txt"))
Fig2c_plot <- ggplot(data=summary1c %>% 
        mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
        mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
        mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
        mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                                  "locfdr50df", "HDMT"))), 
        aes(x=minEff1, y=Incongruous, group=Method)) + 
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") +
  ylab("Num Incongruous") +
  xlim(c(0, 0.20)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

# Fig 2D
summary1d <- fread(here::here("Fig1", "output", "Fig1d_summary.txt"))
Fig2d_plot <- ggplot(data=summary1d  %>% 
        mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
        mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
        mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
        mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                                  "locfdr50df", "HDMT"))),
  aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("Num Incongruous") +
  xlim(c(0, 0.08)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))





