# Collect results and plot Supp Figure 17

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig17/plot_sfig17.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig17/plot_sfig17.R")

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
outputDir <- here::here("SuppFig17", "output/")
names17a <- paste0(outputDir, "Fig17A_aID", 1:600, ".txt")
names17b <- paste0(outputDir, "Fig17B_aID", 1:300, ".txt")
#-----------------------------------------#

# read raw output files
res17a <- c()
for (file_it in 1:length(names17a)) {
    tempRes <- tryCatch(fread(names17a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
    tempList <- list(res17a, tempRes)
    res17a <- rbindlist(tempList)
}

# Read 17b
res17b <- c()
for (file_it in 1:length(names17b)) {
    tempRes <- tryCatch(fread(names17b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
    tempList <- list(res17b, tempRes)
    res17b <- rbindlist(tempList)
}

# summarize
summary17a <- summarize_raw(res17a)
summary17b <- summarize_raw(res17b)

# save summaries
write.table(summary17a, paste0(outputDir, "bincor_changeeff_asym120diff.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary17b, paste0(outputDir, "rep2d_changeeff_asym120diff.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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

# read 17a
bincor_asym120diff <- fread(paste0(outputDir, "bincor_changeeff_asym120diff.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 >= 0.3 & minEff1 <= 0.5)

# plot 17a
bincor_asym120diff_fdp_plot <- ggplot(data=bincor_asym120diff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("c-csmGmm FDP (20% Swap)") +
  #ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
bincor_asym120diff_fdp_plot

# plot 17b
bincor_asym120diff_power_plot <- ggplot(data=bincor_asym120diff %>% 
                                          filter(!(Method %in% c("Kernel", "locfdr7df", "locfdr50df", "DACT"))), aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("c-csmGmm Power (20% Swap)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=mycols[c(1,6)]) +
  scale_linetype_manual(values=c(1,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
bincor_asym120diff_power_plot

# read 17c
rep2d_asym120diff <- fread(paste0(outputDir, "rep2d_changeeff_asym120diff.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 >= 0.3 & minEff1 <= 0.5) %>% 
  filter(Method != "HDMT" & Method != "DACT")

# plot 17 c
rep2d_asym120diff_fdp_plot <- ggplot(data=rep2d_asym120diff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("r-csmGmm FDP (20% Swap)") +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[c(1, 3:6)]) +
  scale_linetype_manual(values=c(1, 3:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
rep2d_asym120diff_fdp_plot

# plot 17d
rep2d_asym120diff_power_plot <- ggplot(data=rep2d_asym120diff %>% filter(!(Method %in% c("locfdr50df"))), aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("r-csmGmm Power (20% Swap)") +
  #geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols[c(1, 3, 5, 6)]) +
  scale_linetype_manual(values=c(1, 3, 5, 6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
rep2d_asym120diff_power_plot

# cowplot it
bincor_rep_asym120_plot <- plot_grid(bincor_asym120diff_fdp_plot + theme(legend.position = "none"),
                                     rep2d_asym120diff_fdp_plot + theme(legend.position = "none"),
                                     bincor_asym120diff_power_plot + theme(legend.position = "none"),
                                     rep2d_asym120diff_power_plot + theme(legend.position = "none"),
                                     labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
bincor_rep_asym120_legend <- get_legend(bincor_asym120diff_fdp_plot +  theme(legend.direction="horizontal",
                                                                         legend.justification="center",
                                                                         legend.box.just="bottom"))
plot_grid(bincor_rep_asym120_plot, bincor_rep_asym120_legend, ncol=1, rel_heights=c(1, 0.15))
ggsave('bincor_rep2d_asym.pdf', width=20, height=12)








