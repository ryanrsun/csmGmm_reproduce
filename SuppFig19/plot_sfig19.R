# Collect results and plot Supp Figure 19

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig19/plot_sfig19.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig19/plot_sfig19.R")

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
outputDir <- here::here("SuppFig19", "output/")
names19a <- paste0(outputDir, "SFig19A_aID101", 1:160, ".txt")
names19b <- paste0(outputDir, "SFig19B_aID101", 1:160, ".txt")
#-----------------------------------------#

# read raw output files
res19a <- c()
for (file_it in 1:length(names19a)) {
  tempRes <- tryCatch(fread(names19a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res19a, tempRes)
  res19a <- rbindlist(tempList)
}

# Read 19b
res19b <- c()
for (file_it in 1:length(names19b)) {
  tempRes <- tryCatch(fread(names19b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res19b, tempRes)
  res19b <- rbindlist(tempList)
}

# summarize
summary19a <- summarize_raw(res19a, cor=T)
summary19b <- summarize_raw(res19b, cor=T)
#summary19a <- summarize_raw(res19a)
#summary19b <- summarize_raw(res19b)

# save summaries
write.table(summary19a, paste0(outputDir, "bincor_changepi0_cor5_large.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary19b, paste0(outputDir, "bincor_changepi0_cor8_large.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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

# read S fig 19a
bincor_changepi0_5 <- fread(paste0(outputDir, "bincor_changepi0_cor5_large.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>% mutate(Method = ifelse(Method == "New", "c-csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("c-csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot s fig 19a
bincor_changepi0_5_fdp_plot <- ggplot(data=bincor_changepi0_5,
                                      aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[1:6]) +
  scale_linetype_manual(values=1:6) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP (2D Correlation = 0.5)") +
  #ylim(c(0, 0.45)) + 
  #xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changepi0_5_fdp_plot

# plot s fig 19c
bincor_changepi0_5_power_plot <- ggplot(data=bincor_changepi0_5 %>% filter(Method %in% c("c-csmGmm")),
                                        aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  scale_color_manual(values=mycols[1]) +
  scale_linetype_manual(values=1) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("Power (2D Correlation = 0.5)") +
  ylim(c(0, 0.5)) + 
  #xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changepi0_5_power_plot

# read s fig 19b
setwd(outputDir)
bincor_changepi0_8 <- fread(paste0(outputDir, "bincor_changepi0_cor8_large.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "c-csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("c-csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot bincor changepi0 8 fdp
bincor_changepi0_8_fdp_plot <- ggplot(data=bincor_changepi0_8,
                                     aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[1:6]) +
  scale_linetype_manual(values=1:6) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP (2D Correlation = 0.8)") +
  #ylim(c(0, 0.45)) + 
  #xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changepi0_8_fdp_plot

# plot bincor changepi0 8 power
bincor_changepi0_8_power_plot <- ggplot(data=bincor_changepi0_8 %>% filter(Method %in% c("c-csmGmm")),
                                       aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  scale_color_manual(values=mycols[1]) +
  scale_linetype_manual(values=1) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("Power (2D Correlation = 0.8)") +
  ylim(c(0, 0.5)) + 
  #xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changepi0_8_power_plot


# save
morecore_changeepi0_fdp_plot <- plot_grid(bincor_changepi0_5_fdp_plot + theme(legend.position = "none"),
                                          bincor_changepi0_8_fdp_plot + theme(legend.position = "none"),
                                          bincor_changepi0_5_power_plot + theme(legend.position = "none"),
                                          bincor_changepi0_8_power_plot + theme(legend.position = "none"),
                                          labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
morecore_changeepi0_fdp_legend <- get_legend(bincor_changepi0_5_fdp_plot +  theme(legend.direction="horizontal",
                                                                                 legend.justification="center",
                                                                                 legend.box.just="bottom"))
plot_grid(morecore_changeepi0_fdp_plot, morecore_changeepi0_fdp_legend, ncol=1, rel_heights=c(1, 0.15))
ggsave('morecor_changepi0.pdf', width=20, height=12)











