# Collect results and plot Supp Fig 3

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig3/plot_suppfig3.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig3/plot_suppfig3.R")

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
outputDir <- here::here("SuppFig3", "output")
fnames3a <- paste0(outputDir, "/SuppFig3A_aID", 221:620, ".txt")
fnames3b <- paste0(outputDir, "/SuppFig3B_aID", 1:160, ".txt")

# colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"

# read raw output files
res3a <- c()
for (file_it in 1:length(fnames3a)) {
  tempRes <- tryCatch(fread(fnames3a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res3a, tempRes)
  res3a <- rbindlist(tempList)
}

res3b <- c()
for (file_it in 1:length(fnames3b)) {
  tempRes <- tryCatch(fread(fnames3b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res3b, tempRes)
  res3b <- rbindlist(tempList)
}

# summarize
summarys3a <- summarize_raw(res3a)
summarys3b <- summarize_raw(res3b)

# save summaries
write.table(summarys3a, here::here(outputDir, "SFig3a_summary_final.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summarys3b, here::here(outputDir, "SFig3b_summary_final.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

#------------------------------------------------#
# plotting starts


# S3a 
ind3d_changeeff_noalt <- fread(here::here(outputDir, "SFig3a_summary_final.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  filter(Method != "DACT") %>%
  filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 >= 0.24 & minEff1 < 0.44)

# plot S3a
ind3d_changeeff_noalt_plot <- ggplot(data=ind3d_changeeff_noalt, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  ylab("FDP (3D Pleiotropy, no alt)") +
  xlab("Min Effect Size") +  ylim(c(0, 0.3)) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind3d_changeeff_noalt_plot


# S3b
ind3d_changepi0_noalt <- fread(here::here(outputDir, "SFig3b_summary_final.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  filter(Method != "DACT") %>%
  filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 <= 0.08)

# plot S3b
ind3d_changepi0_noalt_plot <- ggplot(data=ind3d_changepi0_noalt, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  ylab("FDP (3D Pleiotropy, no alt)") +
  ylim(c(0, 0.3)) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind3d_changepi0_noalt_plot

# put together S3 - just A and B
s2_plot <- plot_grid(ind3d_changeeff_noalt_plot  + theme(legend.position = "none"),
                     ind3d_changepi0_noalt_plot  + theme(legend.position = "none"),
                     labels=c("A", "B"), nrow=1, label_size=22)
s2_plot_legend <- get_legend(ind3d_changeeff_noalt_plot +  theme(legend.direction="horizontal",
                                                                 legend.justification="center",
                                                                 legend.box.just="bottom"))
plot_grid(s2_plot, s2_plot_legend, ncol=1, rel_heights=c(1, 0.15))
#ggsave('supp_3dnoalt.pdf', width=18, height=6)




