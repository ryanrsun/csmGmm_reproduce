# Collect results and plot Supp Fig 6
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig6/origOutput"
#snames6b <- paste0("sim_n1k_j100k_bincor_changepi0_raiseAlt_fixedhdmt_aID", 1:320, ".txt")
snames6b <- paste0("SFig6C_aID", 1:320, ".txt")
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
setwd(outputDir)
sres6b <- c()
for (file_it in 1:length(snames6b)) {
  tempRes <- tryCatch(fread(snames6b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(sres6b, tempRes)
  sres6b <- rbindlist(tempList)
}

# summarize
summarys6b <- summarize_raw(sres6b)

# save summaries
setwd(outputDir)
write.table(summarys6b, "SFig6b_summary_final.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

#------------------------------------------------#
# plotting starts

# SFig 6A is repeated from Fig 3C
setwd(outputDir)
SFig6a_data <- fread("Fig3c_summary_final.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "c-csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("c-csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 >= 0.2 & minEff1 <= 0.39) %>%
  filter(!is.na(Method))

SFig6a_plot <- ggplot(data=SFig6a_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Correlated Test Stat.)") +
  ylim(c(0, 0.45)) + xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))


# S6B data
setwd(outputDir)
SFig6b_data <- fread("SFig6b_summary_final.txt", data.table=F)  %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 <= 0.08)

# s6b plot
SFig6b_plot <- ggplot(data=SFig6b_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# s6c data
setwd(outputDir)
SFig6c_data <- fread("Fig3c_summary_final.txt", data.table=F)  %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# s6c plot
SFig6c_plot <- ggplot(data=SFig6c_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  xlab("Min Effect Size") + ylab("Power") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# plot supp fig 6d
SFig6d_plot <- ggplot(data=SFig6b_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("Power") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))



# put together s6
s6_plot <- plot_grid(SFig6a_plot + theme(legend.position = "none"),
                            SFig6b_plot + theme(legend.position = "none"),
                            SFig6c_plot + theme(legend.position = "none"),
                            SFig6d_plot + theme(legend.position = "none"),
                               labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
s6_plot_legend <- get_legend(SFig6a_plot +  theme(legend.direction="horizontal",
                                                                             legend.justification="center",
                                                                             legend.box.just="bottom"))
plot_grid(s6_plot, s6_plot_legend, ncol=1, rel_heights=c(1, 0.15))
setwd(outputDir)
ggsave('s6_bincor.pdf', width=18, height=12)




