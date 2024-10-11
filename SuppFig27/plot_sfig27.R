# Collect results and plot Supp Figure 27

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig27/plot_sfig27.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig27/plot_sfig27.R")

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
outputDir <- here::here("SuppFig27", "output/")
names27f1 <- paste0(outputDir, "SFig27A_aID", 1:500, "_fit1.txt")
names27f2 <- paste0(outputDir, "SFig27A_aID", 1:500, "_fit3.txt")
names27f3 <- paste0(outputDir, "SFig27A_aID", 1:500, "_fit5.txt")
names27f4 <- paste0(outputDir, "SFig27A_aID", 1:500, "_fit7.txt")
names27fr <- paste0(outputDir, "SFig27A_aID", 1:500, "_fitreg.txt")

# read raw output files
res27f1 <- c()
for (file_it in 1:length(names27f1)) {
  tempRes <- tryCatch(fread(names27f1[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res27f1, tempRes)
  res27f1 <- rbindlist(tempList)
}

# Read 27f2
res27f2 <- c()
for (file_it in 1:length(names27f2)) {
  tempRes <- tryCatch(fread(names27f2[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res27f2, tempRes)
  res27f2 <- rbindlist(tempList)
}

# Read 27f3
res27f3 <- c()
for (file_it in 1:length(names27f3)) {
  tempRes <- tryCatch(fread(names27f3[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res27f3, tempRes)
  res27f3 <- rbindlist(tempList)
}

# Read 27f4
res27f4 <- c()
for (file_it in 1:length(names27f4)) {
  tempRes <- tryCatch(fread(names27f4[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res27f4, tempRes)
  res27f4 <- rbindlist(tempList)
}

# Read 27fr
res27fr <- c()
for (file_it in 1:length(names27fr)) {
  tempRes <- tryCatch(fread(names27fr[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res27fr, tempRes %>% mutate(nRejDACT=NA, nRejHDMT=NA, nRejKernel=NA, nRej7df=NA, nRej50df=NA))
  res27fr <- rbindlist(tempList)
}

# summarize
summary27f1 <- summarize_raw(res27f1)
summary27f2 <- summarize_raw(res27f2)
summary27f3 <- summarize_raw(res27f3)
summary27f4 <- summarize_raw(res27f4)
summary27fr <- summarize_raw(res27fr)

# save summaries
setwd(outputDir)
write.table(summary27f1, paste0(outputDir, "med2d_changeeff_mbltrue5_use1_largeProp.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary27f2, paste0(outputDir, "med2d_changeeff_mbltrue5_use3_largeProp.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary27f3, paste0(outputDir, "med2d_changeeff_mbltrue5_use5_largeProp.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary27f4, paste0(outputDir, "med2d_changeeff_mbltrue5_use7_largeProp.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary27fr, paste0(outputDir, "med2d_changeeff_mbltrue5_use1_double_largeProp.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')


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


# read data
true5_use1_largeProp <- fread(paste0(outputDir, "med2d_changeeff_mbltrue5_use1_largeProp.txt")) %>%
  filter(Method == "New") %>%
  mutate(Method = "New1")
true5_double_largeProp <- fread(paste0(outputDir, "med2d_changeeff_mbltrue5_use1_double_largeProp.txt")) %>%
  filter(Method == "New") %>%
  mutate(Method = "NewD")
true5_use3_largeProp <- fread(paste0(outputDir, "med2d_changeeff_mbltrue5_use3_largeProp.txt")) %>%
  filter(Method == "New") %>%
  mutate(Method = "New3")
true5_use5_largeProp <- fread(paste0(outputDir, "med2d_changeeff_mbltrue5_use5_largeProp.txt"))
true5_use7_largeProp <- fread(paste0(outputDir, "med2d_changeeff_mbltrue5_use7_largeProp.txt")) %>%
  filter(Method == "New") %>%
  mutate(Method = "New7")

# put together data
true5_largeProp <- rbind(true5_use1_largeProp, true5_double_largeProp, true5_use3_largeProp,
                         true5_use5_largeProp, true5_use7_largeProp) %>%
  filter(Method != "DACTb") %>%
  #filter(Method != "df7") %>% filter(Method != "DACT") %>% filter(Method != "Kernel") %>% filter(Method != "df50") %>%
  #filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm-5", Method)) %>%
  mutate(Method = ifelse(Method == "NewD", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "New1", "csmGmm-1", Method)) %>%
  mutate(Method = ifelse(Method == "New3", "csmGmm-3", Method)) %>%
  mutate(Method = ifelse(Method == "New7", "csmGmm-7", Method)) %>%
  #mutate(Method = factor(Method, levels=c("csmGmm-3", "DACT", "Kernel", "locfdr7df",
  #                                     "locfdr50df", "HDMT", "csmGmm-1", "csmGmm-2", "csmGmm-4"))) %>%
  filter(minEff1 <= 1)

# only csm
true5_largeProp_csm <- true5_largeProp %>%
  filter(Method != "locfdr7df") %>% filter(Method != "DACT") %>% filter(Method != "Kernel") %>% filter(Method != "locfdr50df") %>%
  filter(Method != "HDMT") 
# others
true5_largeProp_others <- true5_largeProp %>%
  filter(Method != "csmGmm-1") %>% filter(Method != "csmGmm-3") %>% filter(Method != "csmGmm-5") %>%
  filter(Method != "csmGmm-7") %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot s fig 27a
true5_largeProp_csm_fdp_plot <- ggplot(data=true5_largeProp_csm, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("FDP (5 Means, Many Alt)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=c(mycols[1], "darkgreen", "purple", "cyan", "orange")) +
  scale_linetype_manual(values=c(1:5)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true5_largeProp_csm_fdp_plot

# plot s fig 27b
true5_largeProp_csm_power_plot <- ggplot(data=true5_largeProp_csm, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("Power (5 Means, Many Alt)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=c(mycols[1], "darkgreen", "purple", "cyan", "orange")) +
  scale_linetype_manual(values=c(1:5)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true5_largeProp_csm_power_plot


#----------------------------------------------------------#

# plot s fig 27c
true5_largeProp_others_fdp_plot <- ggplot(data=true5_largeProp_others, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("FDP (5 Means, Many Alt)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true5_largeProp_others_fdp_plot

# plot s fig 27d
true5_largeProp_others_power_plot <- ggplot(data=true5_largeProp_others, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("Power (3 Means, Many Alt)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true5_largeProp_others_power_plot

#----------------------------------------------------#
# put it together

mbl_true5_largeProp_plot <- plot_grid(true5_largeProp_csm_fdp_plot + theme(legend.position = "none"),
                                      true5_largeProp_csm_power_plot + theme(legend.position = "none"),
                                      true5_largeProp_others_fdp_plot + theme(legend.position = "none"),
                                      true5_largeProp_others_power_plot + theme(legend.position = "none"),
                                      labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
mbl_true5_largeProp_legend1 <- get_legend(true5_largeProp_csm_fdp_plot +  theme(legend.direction="horizontal",
                                                                                legend.justification="center",
                                                                                legend.box.just="bottom"))
mbl_true5_largeProp_legend2 <- get_legend(true5_largeProp_others_fdp_plot +  theme(legend.direction="horizontal",
                                                                                   legend.justification="center",
                                                                                   legend.box.just="bottom"))
plot_grid(mbl_true5_largeProp_plot, mbl_true5_largeProp_legend1, mbl_true5_largeProp_legend2, ncol=1, rel_heights=c(1, 0.1, 0.1))
#ggsave('mbl_true5_largeprop.pdf', width=18, height=12)











