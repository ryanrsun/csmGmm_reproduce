# Collect results and plot Figure S4

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig4/plot_figS4.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig4/plot_figS4.R")

library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(devtools)
library(csmGmm)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("SuppFig4", "output/")

names4aq1 <- paste0(outputDir, "sim_n1k_j100k_power3d_changeeff_power1_raiseAlt_aID", 1:600, ".txt")
names4aq2 <- paste0(outputDir, "sim_n1k_j100k_power3d_changeeff_power2_raiseAlt_aID", 1:600, ".txt")
names4aq3 <- paste0(outputDir, "sim_n1k_j100k_power3d_changeeff_power3_raiseAlt_aID", 1:600, ".txt")
names4aq4 <- paste0(outputDir, "sim_n1k_j100k_power3d_changeeff_power4_raiseAlt_aID", 1:600, ".txt")
names4aq5 <- paste0(outputDir, "sim_n1k_j100k_power3d_changeeff_power5_raiseAlt_aID", 1:600, ".txt")
names4aq6 <- paste0(outputDir, "sim_n1k_j100k_power3d_changeeff_power6_raiseAlt_aID", 1:600, ".txt")
names4aq7 <- paste0(outputDir, "sim_n1k_j100k_power3d_changeeff_power7_raiseAlt_aID", 1:600, ".txt")
names4aq8 <- paste0(outputDir, "sim_n1k_j100k_power3d_changeeff_power8_raiseAlt_aID", 1:600, ".txt")
names4aq9 <- paste0(outputDir, "sim_n1k_j100k_power3d_changeeff_power9_raiseAlt_aID", 1:600, ".txt")
names4aList <- list(names4aq1, names4aq2, names4aq3, names4aq4, names4aq5, names4aq6,
                    names4aq7, names4aq8, names4aq9)
names4bq1 <- paste0(outputDir, "Power_correctionS4B_S1_aID", 1:400, ".txt")
names4bq2 <- paste0(outputDir, "Power_correctionS4B_S2_aID", 1:400, ".txt")
names4bq3 <- paste0(outputDir, "Power_correctionS4B_S3_aID", 1:400, ".txt")
names4bq4 <- paste0(outputDir, "Power_correctionS4B_S4_aID", 1:400, ".txt")
names4bq5 <- paste0(outputDir, "Power_correctionS4B_S5_aID", 1:400, ".txt")
names4bq6 <- paste0(outputDir, "Power_correctionS4B_S6_aID", 1:400, ".txt")
names4bq7 <- paste0(outputDir, "Power_correctionS4B_S7_aID", 1:400, ".txt")
names4bq8 <- paste0(outputDir, "Power_correctionS4B_S8_aID", 1:400, ".txt")
names4bq9 <- paste0(outputDir, "Power_correctionS4B_S8_aID", 1:400, ".txt")
names4bList <- list(names4bq1, names4bq2, names4bq3, names4bq4, names4bq5, names4bq6,
                    names4bq7, names4bq8, names4bq9)

#-----------------------------------------#

# read raw output files
res4a <- c()
for (file_it in 1:length(names4aq1)) {
  for (q_it in 1:9) {
    tempRes <- tryCatch(fread(names4aList[[q_it]][file_it]), error=function(e) e)
    if (class(tempRes)[1] == "simpleError") {next}
    tempList <- list(res4a, tempRes %>% mutate(q = q_it))
    res4a <- rbindlist(tempList)
  }
}

res4b <- c()
for (file_it in 1:length(names4bq1)) {
  for (q_it in 1:9) {
    tempRes <- tryCatch(fread(names4bList[[q_it]][file_it]), error=function(e) e)
    if (class(tempRes)[1] == "simpleError") {next}
    tempList <- list(res4b, tempRes %>% mutate(q = q_it))
    res4b <- rbindlist(tempList)
  }
}

# summarize
for (q_it in 1:9) {
  tempResa <- res4a %>% filter(q == q_it)
  tempSummarya <- summarize_raw(tempResa)
  write.table(tempSummarya, paste0(outputDir, "FigS4a_power", q_it, "_summary.txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')

  tempResb <- res4b %>% filter(q == q_it)
  tempSummaryb <- summarize_raw(tempResb)
  write.table(tempSummaryb, paste0(outputDir, "FigS4b_power", q_it, "_summary.txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
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
FigS4a_fixq <- fread(paste0(outputDir, "FigS4a_power1_summary.txt"), data.table=F)  %>% mutate(q=0.01)
for (q_it in 2:9) {
  tempFile <- fread(paste0(outputDir, "FigS4a_power", q_it, "_summary.txt"), data.table=F) %>% mutate(q = q_it / 100)
  FigS4a_fixq <- rbind(FigS4a_fixq, tempFile)
}
FigS4a_fixq <- FigS4a_fixq %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  rbind(., fread(here::here("Fig3/output/Fig3a_summary.txt")) %>% mutate(q = 0.1) %>% 
        mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
        mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
        mutate(Method = ifelse(Method == "df50", "locfdr50df", Method))) 

# correct the power
FigS4a_correctedPow <- FigS4a_fixq %>% filter(q == 0.1) %>%
  select(minEff1, Method, Power) %>%
  mutate(qval = NA, actualPower=NA, actualfdp=NA) %>%
  filter(minEff1 >= 0.25 & minEff1 < 0.45)
allEffs <- unique(FigS4a_fixq$minEff1)
allMets <- unique(FigS4a_fixq$Method)
for (eff_it in 1:length(allEffs)) {
  tempEff <- allEffs[eff_it]
  for (met_it in 1:length(allMets)) {
    tempMet <- allMets[met_it]
    if (tempMet %in% c("HDMT", "csmGmm")) {next}
    tempDat <- FigS4a_fixq %>% filter(minEff1 == tempEff & Method == tempMet) %>%
      arrange(q)
    minRow <- max(which(tempDat$FDP < 0.1))
    fillRow <- which(FigS4a_correctedPow$minEff1 == tempEff & FigS4a_correctedPow$Method == tempMet)
    FigS4a_correctedPow$qval[fillRow] <- tempDat$q[minRow]
    FigS4a_correctedPow$actualPower[fillRow] <- tempDat$Power[minRow]
    FigS4a_correctedPow$actualfdp[fillRow] <- tempDat$FDP[minRow]
  }
}
FigS4a_correctedPow <- FigS4a_correctedPow %>%
  mutate(actualPower = ifelse(Method %in% c("HDMT", "csmGmm"), Power, actualPower)) %>%
  filter(!(Method %in% c("DACT", "HDMT")))

# plot
FigS4a_plot <- ggplot(data=FigS4a_correctedPow, aes(x=minEff1, y=actualPower, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +  xlab("Min Effect Size") + ylab("Power (3D Pleiotropy)") +
  xlim(c(0.25, 0.45)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))


# Fig S4b
FigS4b_fixq <- fread(paste0(outputDir, "FigS4b_power1_summary.txt"), data.table=F)  %>% mutate(q = 0.01)
for (q_it in 2:9) {
  tempFile <- fread(paste0(outputDir, "FigS4b_power", q_it, "_summary.txt"), data.table=F) %>% mutate(q = q_it / 100)
  FigS4b_fixq <- rbind(FigS4b_fixq, tempFile)
}
FigS4b_fixq <- FigS4b_fixq %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  rbind(., fread(here::here("Fig3/output/Fig3b_summary.txt")) %>% mutate(q = 0.1) %>%
        mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
        mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
        mutate(Method = ifelse(Method == "df50", "locfdr50df", Method))) %>%
    filter(minEff1 <= 0.08)
 
# correct the power
FigS4b_correctedPow <- FigS4b_fixq %>% filter(q == 0.1) %>%
  select(minEff1, Method, Power) %>%
  mutate(qval = NA, actualPower=NA, actualfdp=NA)
allEffs <- unique(FigS4b_fixq$minEff1)
allMets <- unique(FigS4b_fixq$Method)
for (eff_it in 1:length(allEffs)) {
  tempEff <- allEffs[eff_it]
  for (met_it in 1:length(allMets)) {
    tempMet <- allMets[met_it]
    if (tempMet %in% c("DACT", "HDMT", "csmGmm")) {next}
    tempDat <- FigS4b_fixq %>% filter(minEff1 == tempEff & Method == tempMet) %>%
      arrange(q)
    minRow <- which.min(abs(tempDat$FDP - 0.1))
    fillRow <- which(FigS4b_correctedPow$minEff1 == tempEff & FigS4b_correctedPow$Method == tempMet)
    FigS4b_correctedPow$qval[fillRow] <- tempDat$q[minRow]
    FigS4b_correctedPow$actualPower[fillRow] <- tempDat$Power[minRow]
    FigS4b_correctedPow$actualfdp[fillRow] <- tempDat$FDP[minRow]
  }
}
# restrict
FigS4b_correctedPow <- FigS4b_correctedPow %>%
  mutate(actualPower = ifelse(Method %in% c("HDMT", "csmGmm"), Power, actualPower))

# plot
FigS4b_plot <- ggplot(data=FigS4b_correctedPow,
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


# Supp Fig S4c
FigS4c_data <- fread(here::here("Fig3/output/Fig3a_summary.txt"), data.table=F) %>%
  filter(Method != "DACT" & Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method))
FigS4c_data$Method =  factor(FigS4c_data$Method, levels=c("csmGmm", "Kernel", "locfdr7df", "locfdr50df"))

FigS4c_plot <- ggplot(data=FigS4c_data, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  xlab("Min Effect Size") + ylab("Num Incongruous") +
  xlim(c(0.25, 0.45)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))

#Supp Fig S4d
FigS4d_data <- fread(here::here("Fig3/output/Fig3b_summary.txt"), data.table=F)  %>%
  filter(Method != "DACT" & Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method))  %>%
  filter(minEff1 <= 0.08)
FigS4d_data$Method =  factor(FigS4d_data$Method, levels=c("csmGmm", "Kernel", "locfdr7df", "locfdr50df"))

FigS4d_plot <- ggplot(data=FigS4d_data, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  ylab("Num Incongruous") +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))

# put together S4
s4_plot <- plot_grid(FigS4a_plot  + theme(legend.position = "none"),
                      FigS4b_plot + theme(legend.position = "none"),
                     FigS4c_plot  + theme(legend.position = "none"),
                     FigS4d_plot + theme(legend.position = "none"),
                     labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
s4_legend <- get_legend(FigS4d_plot +  theme(legend.direction="horizontal",
                                                                   legend.justification="center",
                                                                   legend.box.just="bottom"))
plot_grid(s4_plot, s4_legend, ncol=1, rel_heights=c(1, 0.15))
#ggsave("s4.pdf",  width=20, height=12)





