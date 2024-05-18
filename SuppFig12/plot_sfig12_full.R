# Collect results and plot Supp Figure 12
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig12/origOutput"
# MAF 1%
names12a <- paste0("sim_n1k_j100k_med2d_changeeff_noalt_maf01_aID", 801:1400, ".txt")
names12b <- paste0("SFig12B_aID", 1:160, ".txt")
names12c <- paste0("sim_n1k_j100k_med2d_changeeff_raisealt_maf01_aID", 801:1400, ".txt")
names12d <- paste0("sim_n1k_j100k_med2d_raisealt_changepi0_maf01_aID", 1:160, ".txt")
# MAF 10%
names12a_maf10 <- paste0("SFig12A1_aID", 1:500, ".txt")
names12b_maf10 <- paste0("SFig12B1_aID", 1:160, ".txt")
names12c_maf10 <- paste0("SFig12C1_aID", 1:500, ".txt")
names12d_maf10 <- paste0("SFig12D1_aID", 1:160, ".txt")

#-----------------------------------------#

# read raw output files
setwd(outputDir)
res12a <- c()
for (file_it in 1:length(names12a)) {
  tempRes <- tryCatch(fread(names12a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res12a, tempRes)
  res12a <- rbindlist(tempList)
}

# Read 12b
res12b <- c()
for (file_it in 1:length(names12b)) {
  tempRes <- tryCatch(fread(names12b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res12b, tempRes)
  res12b <- rbindlist(tempList)
}



# Read 12c
res12c <- c()
for (file_it in 1:length(names12c)) {
  tempRes <- tryCatch(fread(names12c[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res12c, tempRes)
  res12c <- rbindlist(tempList)
}

# Read 12d
res12d <- c()
for (file_it in 1:length(names12d)) {
  tempRes <- tryCatch(fread(names12d[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res12d, tempRes)
  res12d <- rbindlist(tempList)
}

# summarize
summary12a <- summarize_raw(res12a)
summary12b <- summarize_raw(res12b)
summary12c <- summarize_raw(res12c)
summary12d <- summarize_raw(res12d)

# save summaries
setwd(outputDir)
write.table(summary12a, "med2d_changeeff_noAlt_maf01.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary12b, "med2d_changepi0_noAlt_maf01.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary12c, "med2d_changeeff_maf01.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary12d, "med2d_changepi0_maf01.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')


#-------------------------------------------------#
# now for MAF=0.1

# read raw output files
setwd(outputDir)
res12a <- c()
for (file_it in 1:length(names12a_maf10)) {
  tempRes <- tryCatch(fread(names12a_maf10[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  if (ncol(tempRes) == 26) {tempRes <- tempRes %>% mutate(nRej50df = NA)}
  tempList <- list(res12a, tempRes)
  res12a <- rbindlist(tempList, use.names=T)
}

# Read 12b
res12b <- c()
for (file_it in 1:length(names12b_maf10)) {
  tempRes <- tryCatch(fread(names12b_maf10[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res12b, tempRes)
  res12b <- rbindlist(tempList)
}

# Read 12c
res12c <- c()
for (file_it in 1:length(names12c_maf10)) {
  tempRes <- tryCatch(fread(names12c_maf10[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  if (ncol(tempRes) == 26) {tempRes <- tempRes %>% mutate(nRej50df = NA)}
  tempList <- list(res12c, tempRes)
  res12c <- rbindlist(tempList, use.names=T)
}

# Read 12d
res12d <- c()
for (file_it in 1:length(names12d_maf10)) {
  tempRes <- tryCatch(fread(names12d_maf10[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res12d, tempRes)
  res12d <- rbindlist(tempList)
}

# summarize
summary12a <- summarize_raw(res12a)
summary12b <- summarize_raw(res12b)
summary12c <- summarize_raw(res12c)
summary12d <- summarize_raw(res12d)

# save summaries
setwd(outputDir)
write.table(summary12a, "SFig12a_maf10_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary12b, "SFig12b_maf10_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary12c, "SFig12c_maf10_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary12d, "SFig12d_maf10_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

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


# S Fig 12A
setwd(outputDir)
ind2d_noalt_changeeff_maf01 <- fread("med2d_changeeff_noAlt_maf01.txt", data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot 12a
ind2d_changeeff_noalt_maf01_plot <- ggplot(data=ind2d_noalt_changeeff_maf01, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (Mediation, No Alt)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.2)) +
  xlim(c(0.4, 0.6)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_noalt_maf01_plot


# S Fig 12B
setwd(outputDir)
ind2d_noalt_changepi0_maf01 <- fread("med2d_changepi0_noAlt_maf01.txt", data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot S Fig 12B
ind2d_changepi0_noalt_maf01_plot <-ggplot(data=ind2d_noalt_changepi0_maf01, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP (Mediation, No Alt)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changepi0_noalt_maf01_plot


# S Fig 12C
setwd(outputDir)
ind2d_changeeff_maf01 <- fread("med2d_changeeff_maf01.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot s fig 12c
ind2d_changeeff_fdp_maf01_plot <- ggplot(data=ind2d_changeeff_maf01, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (Mediation, With Alt)") +
  ylim(c(0, 0.2)) +
  xlim(c(0.4, 0.6)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_fdp_maf01_plot

# s fig 12 d
setwd(outputDir)
ind2d_changepi0_maf01 <- fread("med2d_changepi0_maf01.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot s fig 12 d
ind2d_changepi0_fdp_maf01_plot <- ggplot(data=ind2d_changepi0_maf01, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP (Mediation, With Alt)") +
  ylim(c(0, 0.45)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changepi0_fdp_maf01_plot

#--------------------------------------------------#
# MAF10%

# S Fig 12A
setwd(outputDir)
ind2d_noalt_changeeff_maf10 <- fread("SFig12a_maf10_summary.txt", data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot 12a
ind2d_changeeff_noalt_maf10_plot <- ggplot(data=ind2d_noalt_changeeff_maf10, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (Mediation, No Alt)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.2)) +
  xlim(c(0.2, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# S Fig 12B
setwd(outputDir)
ind2d_noalt_changepi0_maf10 <- fread("SFig12b_maf10_summary.txt", data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot S Fig 12B
ind2d_changepi0_noalt_maf10_plot <-ggplot(data=ind2d_noalt_changepi0_maf10, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP (Mediation, No Alt)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# S Fig 12C
setwd(outputDir)
ind2d_changeeff_maf10 <- fread("SFig12c_maf10_summary.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot s fig 12c
ind2d_changeeff_fdp_maf10_plot <- ggplot(data=ind2d_changeeff_maf10, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (Mediation, With Alt)") +
  ylim(c(0, 0.2)) +
  xlim(c(0.2, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# s fig 12 d
setwd(outputDir)
ind2d_changepi0_maf10 <- fread("SFig12d_maf10_summary.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot s fig 12 d
ind2d_changepi0_fdp_maf10_plot <- ggplot(data=ind2d_changepi0_maf10, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP (Mediation, With Alt)") +
  ylim(c(0, 0.45)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# save
ind2d_fdp_maf01_plot <- plot_grid(ind2d_changeeff_noalt_maf01_plot + theme(legend.position = "none"),
                            ind2d_changepi0_noalt_maf01_plot + theme(legend.position = "none"),
                            ind2d_changeeff_fdp_maf01_plot + theme(legend.position = "none"),
                            ind2d_changepi0_fdp_maf01_plot + theme(legend.position = "none"),
														ind2d_changeeff_noalt_maf10_plot + theme(legend.position = "none"),
                            ind2d_changepi0_noalt_maf10_plot + theme(legend.position = "none"),
                            ind2d_changeeff_fdp_maf10_plot + theme(legend.position = "none"),
                            ind2d_changepi0_fdp_maf10_plot + theme(legend.position = "none"),
                            labels=c("A", "B", "C", "D", "E", "F", "G", "H"), ncol=2, label_size=22)
ind2d_fdp_maf01_legend <- get_legend(ind2d_changeeff_noalt_maf01_plot +  theme(legend.direction="horizontal",
                                                                   legend.justification="center",
                                                                   legend.box.just="bottom"))
plot_grid(ind2d_fdp_maf01_plot, ind2d_fdp_maf01_legend, ncol=1, rel_heights=c(1, 0.15))

setwd(outputDir)
ggsave("sfig12_full.pdf", width=20, height=24)


