# Plot Supp Fig 2
library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(devtools)
devtools::install.packages("ryanrsun/csmGmm")
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#-----------------------------------------#
# change to where the output files are stored
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/Fig2/origOutput"
#-----------------------------------------#

# define colors
gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"

# Supp Fig 2a
setwd(outputDir)
SFig2a_dat <- fread("Fig2a_power5_summary.txt", data.table=F)  %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot Supp Fig 2a
SFig2a_plot <- ggplot(data=SFig2a_dat %>% filter(minEff1 >= 0.07 & minEff1 <= 0.26),  aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylab("FDP (Mediation, With Alt, q=0.05)") +
  xlab("Min Effect Size") + ylim(c(0, 0.2)) +
  geom_hline(yintercept = 0.05, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))


# Supp Fig 2b
setwd(outputDir)
SFig2b_dat <- fread("Fig2b_power5_summary.txt", data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot supp fig 2b
SFig2b_plot <- ggplot(data=SFig2b_dat, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylab("FDP (Mediation, With Alt, q=0.05)") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  geom_hline(yintercept = 0.05, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))


# Supp Fig 2c
setwd(outputDir)
SFig2c_fixq <- fread("Fig2a_power1_summary.txt", data.table=F) %>%
  mutate(qval=0.01)
for (q_it in 2:5) {
  tempFile <- fread(paste0("Fig2a_power", q_it, "_summary.txt"), data.table=F) %>% mutate(qval = q_it * 0.01)
  SFig2c_fixq <- rbind(SFig2c_fixq, tempFile)
}
SFig2c_fixq <- SFig2c_fixq %>%  
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# correct the power
SFig2c_correctedPos <- SFig2c_fixq %>% filter(qval == 0.05) %>%
  select(minEff1, Method, Power) %>%
  mutate(qval = NA, actualPower=NA, actualfdp=NA)
allEffs <- unique(SFig2c_fixq$minEff1)
allMets <- unique(SFig2c_fixq$Method)
for (eff_it in 1:length(allEffs)) {
  tempEff <- allEffs[eff_it]
  for (met_it in 1:length(allMets)) {
    tempMet <- allMets[met_it]
    if (tempMet %in% c("HDMT", "csmGmm")) {next}
    tempDat <- SFig2c_fixq %>% filter(minEff1 == tempEff & Method == tempMet) %>%
      arrange(qval)
    minRow <- max(which(tempDat$FDP < 0.05))
    fillRow <- which(SFig2c_correctedPos$minEff1 == tempEff & SFig2c_correctedPos$Method == tempMet)
    SFig2c_correctedPos$qval[fillRow] <- tempDat$qval[minRow]
    SFig2c_correctedPos$actualPower[fillRow] <- tempDat$Power[minRow]
    SFig2c_correctedPos$actualfdp[fillRow] <- tempDat$FDP[minRow]
  }
}
# restrict
SFig2c_correctedPos <- SFig2c_correctedPos %>%
  mutate(actualPower = ifelse(Method %in% c("HDMT", "csmGmm"), Power, actualPower))

# plot Supp fig 2c
SFig2c_plot <- ggplot(data=SFig2c_correctedPos,
                                        aes(x=minEff1, y=actualPower, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  xlab("Min Effect Size") + ylab("Power (at q=0.05)") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# Supp Fig 2d
setwd(outputDir)
Sfig2d_fixq <- fread("Fig2b_power1_summary.txt", data.table=F) %>%
  mutate(qval=0.01)
for (q_it in 2:5) {
  tempFile <- fread(paste0("Fig2b_power", q_it, "_summary.txt"), data.table=F) %>% mutate(qval = q_it * 0.01)
  Sfig2d_fixq <- rbind(Sfig2d_fixq, tempFile)
}
Sfig2d_fixq <- Sfig2d_fixq %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# correct the power
Sfig2d_fixq_correctedPow <- Sfig2d_fixq %>% filter(qval == 0.05) %>%
  select(minEff1, Method, Power) %>%
  mutate(qval = NA, actualPower=NA, actualfdp=NA)
allEffs <- unique(Sfig2d_fixq$minEff1)
allMets <- unique(Sfig2d_fixq$Method)
for (eff_it in 1:length(allEffs)) {
  tempEff <- allEffs[eff_it]
  for (met_it in 1:length(allMets)) {
    tempMet <- allMets[met_it]
    if (tempMet %in% c("HDMT", "csmGmm")) {next}
    tempDat <- Sfig2d_fixq %>% filter(minEff1 == tempEff & Method == tempMet) %>%
      arrange(qval)
    minRow <- which.min(abs(tempDat$FDP - 0.05))
    fillRow <- which(Sfig2d_fixq_correctedPow$minEff1 == tempEff & Sfig2d_fixq_correctedPow$Method == tempMet)
    Sfig2d_fixq_correctedPow$qval[fillRow] <- tempDat$qval[minRow]
    Sfig2d_fixq_correctedPow$actualPower[fillRow] <- tempDat$Power[minRow]
    Sfig2d_fixq_correctedPow$actualfdp[fillRow] <- tempDat$FDP[minRow]
  }
}
# restrict
Sfig2d_fixq_correctedPow <- Sfig2d_fixq_correctedPow %>%
  mutate(actualPower = ifelse(Method %in% c("HDMT", "csmGmm"), Power, actualPower)) %>%
  filter(minEff1 <= 0.08)

# plot ind2d changeeff corrected power
Sfig2d_plot <- ggplot(data=Sfig2d_fixq_correctedPow,
                                        aes(x=minEff1, y=actualPower, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylab("Power (at q=0.05)") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))








