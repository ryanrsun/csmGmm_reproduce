# Collect results and plot Supp Figure 1

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig1/plot_suppfig1.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig1/plot_suppfig1.R")

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

outputDir <- here::here("SuppFig1/output/")
snames1a <- paste0(outputDir, "sim_n1k_j100k_med2d_changeeff_flipz_noalt_aID", 1:800, ".txt")
snames1b <- paste0(outputDir, "sim_n1k_j100k_med2d_raisealt_changepi0_flipz_noalt_aID", 1:320, ".txt")
snames1c <- paste0(outputDir, "SuppFig1C_aID", 1:800, ".txt")
snames1d <- paste0(outputDir, "SuppFig1D_aID", 1:320, ".txt")
#-----------------------------------------#

# colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"


# read raw output files
# supp 1a
sres1a <- c()
for (file_it in 1:length(snames1a)) {
  tempRes <- tryCatch(fread(snames1a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(sres1a, tempRes)
  sres1a <- rbindlist(tempList)
}

# Read s1b
sres1b <- c()
for (file_it in 1:length(snames1b)) {
  tempRes <- tryCatch(fread(snames1b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(sres1b, tempRes)
  sres1b <- rbindlist(tempList)
}

# Read s1c
sres1c <- c()
for (file_it in 1:length(snames1c)) {
  tempRes <- tryCatch(fread(snames1c[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(sres1c, tempRes)
  sres1c <- rbindlist(tempList)
}

# Read 1d
sres1d <- c()
for (file_it in 1:length(snames1d)) {
  tempRes <- tryCatch(fread(snames1d[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(sres1d, tempRes)
  sres1d <- rbindlist(tempList)
}


# summarize
summarys1a <- summarize_raw(sres1a)
summarys1b <- summarize_raw(sres1b)
summarys1c <- summarize_raw(sres1c)
summarys1d <- summarize_raw(sres1d)

# save summaries
write.table(summarys1a, paste0(outputDir, "SFig1a_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summarys1b, paste0(outputDir, "SFig1b_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summarys1c, paste0(outputDir, "SFig1c_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summarys1d, paste0(outputDir, "SFig1d_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

#------------------------------------------------#
# plotting starts

# Supp fig - swap alleles
# s1a data
SFig1a_data <- fread(paste0(outputDir, "SFig1a_summary.txt"), data.table=F)  %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot s1a
SFig1a_plot <- ggplot(data=SFig1a_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylab("FDP (Mediation, No Alt)") +
  xlab("Min Effect Size") + ylim(c(0, 0.2)) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# s1b data
setwd(outputDir)
SFig1b_data <- fread(paste0(outputDir, "SFig1b_summary.txt"), data.table=F)  %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 <= 0.08)

# s1b plot
SFig1b_plot <- ggplot(data=SFig1b_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylab("FDP (Mediation, No Alt)") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# s1c data
setwd(outputDir)
SFig1c_data <- fread(paste0(outputDir, "SFig1c_summary.txt"), data.table=F)  %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# s1c plot
SFig1c_plot <- ggplot(data=SFig1c_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylab("FDP (Mediation)") +
  xlab("Min Effect Size") + ylim(c(0, 0.2)) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# s1d data
setwd(outputDir)
SFig1d_data <- fread(paste0(outputDir, "SFig1d_summary.txt"), data.table=F)  %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 <= 0.08)

# plot
SFig1d_plot <- ggplot(data=SFig1d_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylab("FDP (Mediation)") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# put together
S1_plot <- plot_grid(SFig1a_plot + theme(legend.position = "none"),
                       SFig1b_plot + theme(legend.position = "none"),
                       SFig1c_plot + theme(legend.position = "none"),
                       SFig1d_plot + theme(legend.position = "none"),
                        labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
S1_legend <- get_legend(SFig1a_plot +  theme(legend.direction="horizontal",
                                                                     legend.justification="center",
                                                                     legend.box.just="bottom"))
plot_grid(S1_plot, S1_legend, ncol=1, rel_heights=c(1, 0.15))
#ggsave(paste0(outputDir, 'SFig1.pdf'), width=18, height=12)






