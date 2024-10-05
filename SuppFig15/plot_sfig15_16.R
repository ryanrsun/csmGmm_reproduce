# Collect results and plot Supp Figure 15 and 16

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig15/plot_sfig15_16.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig15/plot_sfig15_16.R")

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
outputDir <- here::here("SuppFig15", "output/")
names15a <- paste0(outputDir, "SFig15A_aID", 1:160, ".txt")
names15b <- paste0(outputDir, "SFig15B_aID", 1:400, ".txt")
names16a <- paste0(outputDir, "sim_n1k_j100k_med2d_changepi0_asym150diff_aID", 1:160, ".txt")
names16b <- paste0(outputDir, "sim_n1k_j100k_med2d_changepi0_asym200diff_aID", 1:160, ".txt")
#-----------------------------------------#

# read raw output files
res15a <- c()
for (file_it in 1:length(names15a)) {
  tempRes <- tryCatch(fread(names15a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res15a, tempRes)
  res15a <- rbindlist(tempList)
}

# Read 15b
res15b <- c()
for (file_it in 1:length(names15b)) {
  tempRes <- tryCatch(fread(names15b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res15b, tempRes)
  res15b <- rbindlist(tempList)
}

# Read 16a
res16a <- c()
for (file_it in 1:length(names16a)) {
  tempRes <- tryCatch(fread(names16a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res16a, tempRes)
  res16a <- rbindlist(tempList)
}

# Read 16b
res16b <- c()
for (file_it in 1:length(names16b)) {
  tempRes <- tryCatch(fread(names16b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res16b, tempRes)
  res16b <- rbindlist(tempList)
}

# summarize
summary15a <- summarize_raw(res15a)
summary15b <- summarize_raw(res15b)
summary16a <- summarize_raw(res16a)
summary16b <- summarize_raw(res16b)

# save summaries
write.table(summary15a, paste0(outputDir, "ind2d_changepi0_asym120diff_full.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary15b, paste0(outputDir, "ind2d_changeeff_asym120diff_full.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary16a, paste0(outputDir, "ind2d_changepi0_asym150diff_full.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary16b, paste0(outputDir, "ind2d_changepi0_asym200diff_full.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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


# S Fig 15a
changepi0_120 <- fread(paste0(outputDir, "ind2d_changepi0_asym120diff_full.txt")) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) 

# plot S Fig 15a
changepi0_120_fdp_plot <- ggplot(data=changepi0_120, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0)"))) +
  ylab("FDP (20% Swap)") +
  ylim(c(0, 0.4)) +
  xlim(c(0, 0.08)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
changepi0_120_fdp_plot

# S Fig 15b
setwd(outputDir)
changeeff_120 <- fread(paste0(outputDir, "ind2d_changeeff_asym120diff_full.txt")) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) 

# Plot S Fig 15b
changeeff_120_fdp_plot <- ggplot(data=changeeff_120, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (20% Swap)") +
  ylim(c(0, 0.4)) +
  xlim(c(0.25, 0.45)) +  
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
changeeff_120_fdp_plot


# Plot S Fig 15c
changepi0_120_power_plot <- ggplot(data=changepi0_120 %>% filter(Method != "DACT"), aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0)"))) +
  ylab("Power (20% Swap)") +
  xlim(c(0, 0.08)) +
  scale_color_manual(values=mycols[c(1,3,4,5,6)]) +
  scale_linetype_manual(values=c(1,3,4,5,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
changepi0_120_power_plot

# Plot S Fig 15d
changeeff_120_power_plot <- ggplot(data=changeeff_120, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("Power (20% Swap)") +
  xlim(c(0.25, 0.45)) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
changeeff_120_power_plot


# Put together S Fig 15
sfig15_plot <- plot_grid(changepi0_120_fdp_plot + theme(legend.position = "none"),
                                      changeeff_120_fdp_plot + theme(legend.position = "none"),
                                      changepi0_120_power_plot + theme(legend.position = "none"),
                                      changeeff_120_power_plot + theme(legend.position = "none"),
                                      labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
sfig15_legend <- get_legend(changepi0_120_fdp_plot +  theme(legend.direction="horizontal",
                                                                                legend.justification="center",
                                                                                legend.box.just="bottom"))

plot_grid(sfig15_plot, sfig15_legend, ncol=1, rel_heights=c(1, 0.1))
ggsave('asymmetric_ind120.pdf', width=18, height=12)

#---------------------------------------------------------#

# load S Fig 16a
changepi0_150 <- fread(paste0(outputDir, "ind2d_changepi0_asym150diff_full.txt")) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) 

# plot S Fig 16a
changepi0_150_fdp_plot <- ggplot(data=changepi0_150, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0)"))) +
  ylab("FDP (50% Swap)") +
  ylim(c(0, 0.4)) +
  xlim(c(0, 0.08)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.2) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
changepi0_150_fdp_plot

# plot S Fig 16b
changepi0_150_power_plot <- ggplot(data=changepi0_150 %>% filter(Method != "DACT"), aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0)"))) +
  ylab("Power (50% Swap)") +
  xlim(c(0, 0.08)) +
  scale_color_manual(values=mycols[c(1,3,4,5,6)]) +
  scale_linetype_manual(values=c(1,3,4,5,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
changepi0_150_power_plot

# cowplot the 150 and 200 changepi
changepi150_plot <- plot_grid(changepi0_150_fdp_plot + theme(legend.position = "none"),
                         changepi0_150_power_plot + theme(legend.position = "none"),
                         labels=c("A", "B"), nrow=1, label_size=22)
changepi150_legend <- get_legend(changepi0_150_fdp_plot +  theme(legend.direction="horizontal",
                                                            legend.justification="center",
                                                            legend.box.just="bottom"))
#true3_largeProp_legend <- plot_grid(mbl_true3_largeProp_legend1, mbl_true3_largeProp_legend2, ncol=2)
plot_grid(changepi150_plot, changepi150_legend, ncol=1, rel_heights=c(1, 0.1))
ggsave('asymmetric_changepi150.pdf', width=18, height=10)













