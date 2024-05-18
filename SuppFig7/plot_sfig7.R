# Collect results and plot Supp Fig 7
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig7/output"
snames7b <- paste0("SFig7B_aID", 1:160, ".txt")
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
sres7b <- c()
for (file_it in 1:length(snames7b)) {
  tempRes <- tryCatch(fread(snames7b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(sres7b, tempRes)
  sres7b <- rbindlist(tempList)
}

# summarize
summarys7b <- summarize_raw(sres7b)

# save summaries
setwd(outputDir)
write.table(summarys7b, "SFig7b_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

#------------------------------------------------#
# plotting starts


# S7a 
SFig7a_data <- fread("Fig3d_summary.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  filter(Method != "DACT") %>%
  filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))   %>%
  #filter(Method %in% c("locfdr7df", "locfdr50df", "Kernel", "csmGmm")) %>%
  filter(minEff1 >= 0.3 & minEff1 <= 0.49)

# plot S7a
SFig7a_plot <- ggplot(data=SFig7a_data %>% mutate(FDP = ifelse(Method %in% c("HDMT", "DACT"), 1, FDP)), aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylim(c(0, 0.2)) + xlim(c(0.3, 0.5)) +
  xlab("Min Effect Size") + ylab("FDP") +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))


# S7c 
# plot
SFig7c_plot <- ggplot(data=SFig7a_data %>% mutate(Power = ifelse(Method %in% c("HDMT", "DACT"), 1, Power)), aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  ylim(c(0, 1)) + xlim(c(0.3, 0.5)) +
  xlab("Min Effect Size") + ylab("Power") +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))


# S7b
setwd(outputDir)
SFig7b_data <- fread("SFig7b_summary.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  filter(Method != "DACT") %>%
  filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 <= 0.08)

# plot
SFig7b_plot <- ggplot(data=SFig7b_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP") +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# S7d 
# plot
SFig7d_plot <- ggplot(data=SFig7b_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("Power") +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# put together s7
s7_plot <- plot_grid(SFig7a_plot + theme(legend.position = "none"),
                           SFig7b_plot + theme(legend.position = "none"),
                           SFig7c_plot + theme(legend.position = "none"),
                           SFig7d_plot + theme(legend.position = "none"),
                            labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
s7_legend <- get_legend(SFig7a_plot +  theme(legend.direction="horizontal",
                                                                       legend.justification="center",
                                                                       legend.box.just="bottom"))
plot_grid(s7_plot, s7_legend, ncol=1, rel_heights=c(1, 0.15))
setwd(outputDir)
ggsave('s7.pdf', width=18, height=12)





