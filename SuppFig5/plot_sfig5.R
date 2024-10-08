# Collect results and plot Supp Fig 5

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig5/plot_sfig5.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig5/plot_sfig5.R")

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
outputDir <- here::here("SuppFig5", "output")
names5a <- here::here(outputDir, "SFig5A_aID", 1:880, ".txt")
names5b <- here::here(outputDir, "SFig5B_aID", 1:320, ".txt")

# read raw output files
res5a <- c()
for (file_it in 1:length(names5a)) {
  tempRes <- tryCatch(fread(names5a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res5a, tempRes)
  res5a <- rbindlist(tempList)
}

# Read 5b
res5b <- c()
for (file_it in 1:length(names5b)) {
  tempRes <- tryCatch(fread(names5b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res5b, tempRes)
  res5b <- rbindlist(tempList)
}

# summarize
summary5a <- summarize_raw(res5a, FDP2=TRUE)
summary5b <- summarize_raw(res5b, FDP2=TRUE)

# save summaries
write.table(summary5a, here::here(outputDir, "SFig5a_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary5b, here::here(outputDir, "SFig5b_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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

# Supp Fig 5a
SFig5a_dat <- fread(here::here(outputDir, "SFig5a_summary.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 >= 0.1 & minEff1 <= 0.3)

# Plot S Fig 5a
SFig5a_plot <- ggplot(data=SFig5a_dat, aes(x=minEff1, y=sig2FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylim(c(0, 0.7)) + 
  #xlim(c(0.2, 0.4)) +
  xlab("Min Effect Size") + ylab("FDP (NCP Threshold=2)") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# Supp Fig 5b
SFig5b_dat <- fread(here::here(outputDir, "SFig5b_summary.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot S Fig 5b
SFig5b_plot <- ggplot(data=SFig5b_dat, aes(x=minEff1, y=sig2FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylim(c(0, 0.7)) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0)"))) +
  ylab("FDP (NCP Threshold=2)") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# S Fig 5c
SFig5c_dat <- fread(here::here(outputDir, "SFig5a_summary.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 >= 0.1 & minEff1 <= 0.3)

# plot S Fig 5c
SFig5c_plot <- ggplot(data=SFig5c_dat, aes(x=minEff1, y=sig2Pow, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  xlab("Min Effect Size") + ylab("Power (NCP Threshold=2)") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# plot S Fig 5d
SFig5d_plot <- ggplot(data=SFig5b_dat, aes(x=minEff1, y=sig2Pow, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0)"))) +
  ylab("Power (NCP Threshold=2)") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# put together S5
s5_plot <- plot_grid(SFig5a_plot + theme(legend.position = "none"),
                     SFig5b_plot + theme(legend.position = "none"),
                     SFig5c_plot + theme(legend.position = "none"),
                     SFig5d_plot + theme(legend.position = "none"),
                     labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
s5_plot_legend <- get_legend(SFig5a_plot +  theme(legend.direction="horizontal",
                                                                       legend.justification="center",
                                                                       legend.box.just="bottom"))
plot_grid(s5_plot, s5_plot_legend, ncol=1, rel_heights=c(1, 0.15))
#ggsave('sfig5.pdf', width=18, height=12)







