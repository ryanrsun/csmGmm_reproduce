# Collect results and plot Supp Figure 10

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig10/plot_sfig10.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig10/plot_sfig10.R")

library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(devtools)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("SuppFig10", "output/")
names10aq1 <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff_powerS1_aID", 1:2000, ".txt")
names10aq2 <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff_powerS2_aID", 1:2000, ".txt")
names10aq3 <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff_powerS3_aID", 1:2000, ".txt")
names10aq4 <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff_powerS4_aID", 1:2000, ".txt")
names10aq5 <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff_powerS5_aID", 1:2000, ".txt")
names10aq6 <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff_powerS6_aID", 1:2000, ".txt")
names10aq7 <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff_powerS7_aID", 1:2000, ".txt")
names10aq8 <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff_powerS8_aID", 1:2000, ".txt")
names10aq9 <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff_powerS9_aID", 1:2000, ".txt")
names10aList <- list(names10aq1, names10aq2, names10aq3, names10aq4, names10aq5, names10aq6,
                    names10aq7, names10aq8, names10aq9)
names10bq1 <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff_powerS1_aID", 1:2000, ".txt")
names10bq2 <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff_powerS2_aID", 1:2000, ".txt")
names10bq3 <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff_powerS3_aID", 1:2000, ".txt")
names10bq4 <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff_powerS4_aID", 1:2000, ".txt")
names10bq5 <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff_powerS5_aID", 1:2000, ".txt")
names10bq6 <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff_powerS6_aID", 1:2000, ".txt")
names10bq7 <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff_powerS7_aID", 1:2000, ".txt")
names10bq8 <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff_powerS8_aID", 1:2000, ".txt")
names10bq9 <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff_powerS9_aID", 1:2000, ".txt")
names10bList <- list(names10bq1, names10bq2, names10bq3, names10bq4, names10bq5, names10bq6,
                    names10bq7, names10bq8, names10bq9)

#-----------------------------------------#

# read raw output files
res10a <- c()
for (file_it in 1:length(names10aq1)) {
  for (q_it in 1:9) {
    tempRes <- tryCatch(fread(names10aList[[q_it]][file_it]), error=function(e) e)
    if (class(tempRes)[1] == "simpleError") {next}
    tempList <- list(res10a, tempRes %>% mutate(q = q_it))
    res10a <- rbindlist(tempList)
  }
}

res10b <- c()
for (file_it in 1:length(names10bq1)) {
  for (q_it in 1:9) {
    tempRes <- tryCatch(fread(names10bList[[q_it]][file_it]), error=function(e) e)
    if (class(tempRes)[1] == "simpleError") {next}
    tempList <- list(res10b, tempRes %>% mutate(q = q_it))
    res10b <- rbindlist(tempList)
  }
}

# summarize
for (q_it in 1:9) {
  tempResa <- res10a %>% filter(q == q_it)
  tempSummarya <- summarize_raw(tempResa)
  write.table(tempSummarya, paste0(outputDir, "Fig10a_power", q_it, "_summary.txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')

  tempResb <- res10b %>% filter(q == q_it)
  tempSummaryb <- summarize_raw(tempResb)
  write.table(tempSummaryb, paste0(outputDir, "Fig10b_power", q_it, "_summary.txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
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

# 5d power
# load all the files needed for corrected power
ind5d_changeeff <- fread(paste0(outputDir, "Fig10a_power1_summary.txt")) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method))
ind5d_changeeff_fixq <- fread(paste0(outputDir, "Fig10a_power1_summary.txt"), data.table=F) %>%
  mutate(qval=0.01)
for (q_it in 2:9) {
  tempFile <- fread(paste0(outputDir, "Fig10a_power", q_it, "_summary.txt"), data.table=F) %>% mutate(qval = q_it * 0.01)
  ind5d_changeeff_fixq <- rbind(ind5d_changeeff_fixq, tempFile)
}
ind5d_changeeff_fixq <- ind5d_changeeff_fixq %>% select(-Incongruous) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  select(minEff1, Method, nRej, Power, FDP, numNA, qval) %>%
  rbind(., ind5d_changeeff %>% select(minEff1, Method, nRej, Power, FDP, numNA) %>% mutate(qval = 0.1))

# correct the power
ind5d_changeeff_correctedPow <- ind5d_changeeff_fixq %>% filter(qval == 0.1) %>%
  select(minEff1, Method, Power) %>%
  mutate(qval = NA, actualPower=NA, actualfdp=NA)
allEffs <- unique(ind5d_changeeff_fixq$minEff1)
allMets <- unique(ind5d_changeeff_fixq$Method)
for (eff_it in 1:length(allEffs)) {
  tempEff <- allEffs[eff_it]
  for (met_it in 1:length(allMets)) {
    tempMet <- allMets[met_it]
    if (tempMet %in% c("HDMT")) {next}
    tempDat <- ind5d_changeeff_fixq %>% filter(minEff1 == tempEff & Method == tempMet) %>%
      arrange(qval)
    minRow <- max(which(tempDat$FDP < 0.1))
    #minRow <- which.min(abs(tempDat$FDP - 0.1))
    fillRow <- which(ind5d_changeeff_correctedPow$minEff1 == tempEff & ind5d_changeeff_correctedPow$Method == tempMet)
    ind5d_changeeff_correctedPow$qval[fillRow] <- tempDat$qval[minRow]
    ind5d_changeeff_correctedPow$actualPower[fillRow] <- tempDat$Power[minRow]
    ind5d_changeeff_correctedPow$actualfdp[fillRow] <- tempDat$FDP[minRow]
  }
}
# restrict
ind5d_changeeff_correctedPow <- ind5d_changeeff_correctedPow %>%
  filter(minEff1 <= 1) %>%
  mutate(actualPower = ifelse(Method %in% c("HDMT", "csmGmm"), Power, actualPower))

# plot ind6d changeeff corrected power
ind5d_changeeff_correctedPowPlot <- ggplot(data=ind5d_changeeff_correctedPow %>% filter(!is.na(actualPower)),
                                           aes(x=minEff1, y=actualPower, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") + ylab("Power (5D Pleiotropy)") +
  # ylim(c(0, 0.8)) + xlim(c(0, 0.2)) +
  scale_color_manual(values=mycols[c(1, 3:6)]) +
  scale_linetype_manual(values=c(1, 3:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind5d_changeeff_correctedPowPlot


# 6d power
# load all the files needed for corrected power
ind6d_changeeff <- fread(paste0(outputDir, "med2d_changeeff_ind6d_v2.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method))
ind6d_changeeff_fixq <- fread(paste0(outputDir, "Fig10b_power1_summary.txt"), data.table=F) %>%
  mutate(qval=0.01)
for (q_it in 2:9) {
  tempFile <- fread(paste0(outputDir, "Fig10b_power", q_it, "_summary.txt"), data.table=F) %>% mutate(qval = q_it * 0.01)
  ind6d_changeeff_fixq <- rbind(ind6d_changeeff_fixq, tempFile)
}
ind6d_changeeff_fixq <- ind6d_changeeff_fixq %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  select(minEff1, Method, nRej, Power, FDP, numNA, qval) %>%
  rbind(., ind6d_changeeff %>% select(minEff1, Method, nRej, Power, FDP, numNA) %>% mutate(qval=0.1))  

# correct the power
ind6d_changeeff_correctedPow <- ind6d_changeeff_fixq %>% filter(qval == 0.1) %>%
  select(minEff1, Method, Power) %>%
  mutate(qval = NA, actualPower=NA, actualfdp=NA)
allEffs <- unique(ind6d_changeeff_fixq$minEff1)
allMets <- unique(ind6d_changeeff_fixq$Method)
for (eff_it in 1:length(allEffs)) {
  tempEff <- allEffs[eff_it]
  for (met_it in 1:length(allMets)) {
    tempMet <- allMets[met_it]
    if (tempMet %in% c("HDMT")) {next}
    tempDat <- ind6d_changeeff_fixq %>% filter(minEff1 == tempEff & Method == tempMet) %>%
      arrange(qval)
    minRow <- max(which(tempDat$FDP < 0.1))
    #minRow <- which.min(abs(tempDat$FDP - 0.1))
    fillRow <- which(ind6d_changeeff_correctedPow$minEff1 == tempEff & ind6d_changeeff_correctedPow$Method == tempMet)
    ind6d_changeeff_correctedPow$qval[fillRow] <- tempDat$qval[minRow]
    ind6d_changeeff_correctedPow$actualPower[fillRow] <- tempDat$Power[minRow]
    ind6d_changeeff_correctedPow$actualfdp[fillRow] <- tempDat$FDP[minRow]
  }
}
# restrict
ind6d_changeeff_correctedPow <- ind6d_changeeff_correctedPow %>%
  filter(minEff1 <= 1) %>%
  mutate(actualPower = ifelse(Method %in% c("HDMT", "csmGmm"), Power, actualPower)) %>%
  filter(!(Method %in% c("DACT", "HDMT")))

# plot ind6d changeeff corrected power
ind6d_changeeff_correctedPowPlot <- ggplot(data=ind6d_changeeff_correctedPow,
                                           aes(x=minEff1, y=actualPower, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") + ylab("Power (6D Pleiotropy)") +
 # ylim(c(0, 0.8)) + xlim(c(0, 0.2)) +
  scale_color_manual(values=mycols[c(1, 3:6)]) +
  scale_linetype_manual(values=c(1, 3:6)) +
  #scale_x_continuous(breaks=seq(from=0.25, to=0.35, by=0.05)) +
  #xlim(c(0.25, 0.35)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind6d_changeeff_correctedPowPlot


# save power
ind56d_power_plot <- plot_grid(ind5d_changeeff_correctedPowPlot + theme(legend.position = "none"),
                               ind6d_changeeff_correctedPowPlot + theme(legend.position = "none"),
                             labels=c("A", "B", "C", "D"), nrow=1, label_size=22)
ind56d_power_legend <- get_legend(ind6d_changeeff_correctedPowPlot +  theme(legend.direction="horizontal",
                                                                  legend.justification="center",
                                                                  legend.box.just="bottom"))
plot_grid(ind56d_power_plot, ind56d_power_legend, ncol=1, rel_heights=c(1, 0.15))
#ggsave('ind56d_power.pdf', width=20, height=12)



