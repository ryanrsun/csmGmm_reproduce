# Collect results and plot Supp Fig 8
library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(devtools)
devtools::install_github("ryanrsun/csmGmm")
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#-----------------------------------------#
# change to where the output files are stored
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig8/origOutput"
names8a <- paste0("sim_n1k_j100k_med2d_noalt_changepi0_randomeff_pi0_aID", 1:160, ".txt")
names8b <- paste0("sim_n1k_j100k_med2d_raisealt_changepi0_randomeff_pi0_aID", 1:160, ".txt")
names8c <- paste0("sim_n1k_j100k_ind3d_changepi0_raiseAlt_randomeff_aID", 1:160, ".txt")
#-----------------------------------------#

# read raw output files
setwd(outputDir)
res8a <- c()
for (file_it in 1:length(names8a)) {
  tempRes <- tryCatch(fread(names8a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res8a, tempRes)
  res8a <- rbindlist(tempList)
}

# Read 8b
res8b <- c()
for (file_it in 1:length(names8b)) {
  tempRes <- tryCatch(fread(names8b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res8b, tempRes)
  res8b <- rbindlist(tempList)
}

# Read 8c
res8c <- c()
for (file_it in 1:length(names8c)) {
  tempRes <- tryCatch(fread(names8c[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res8c, tempRes)
  res8c <- rbindlist(tempList)
}

# summarize
summary8a <- summarize_raw(res8a)
summary8b <- summarize_raw(res8b)
summary8c <- summarize_raw(res8c)

# save summaries
setwd(outputDir)
write.table(summary8a, "med2d_noalt_changepi0_randomeff.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary8b, "med2d_changepi0_randomeff.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary8c, "ind3d_changepi0_randomeff.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

#------------------------------------------------#
# plotting starts

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"


# no effect change pi0
setwd(outputDir)
ind2d_noalt_changepi0_randomeff <- fread("med2d_noalt_changepi0_randomeff.txt", data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot no effect change pi0
ind2d_changepi0_noalt_fdp_randomeff_plot <-ggplot(data=ind2d_noalt_changepi0_randomeff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
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
ind2d_changepi0_noalt_fdp_randomeff_plot

# change pi0 data
setwd(outputDir)
ind2d_changepi0_randomeff <- fread("med2d_changepi0_randomeff.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot changepi0 fdp
ind2d_changepi0_fdp_randomeff_plot <- ggplot(data=ind2d_changepi0_randomeff, aes(x=minEff1, y=FDP, group=Method)) +
  #geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
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
ind2d_changepi0_fdp_randomeff_plot


# plot changepi0 power
ind2d_changepi0_power_randomeff_plot <- ggplot(data=ind2d_changepi0_randomeff %>% filter(Method != "DACT"), aes(x=minEff1, y=Power, group=Method)) +
  #geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("Power (Mediation, With Alt)") +
  ylim(c(0, 0.45)) +
 # geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[c(1, 3:6)]) +
  scale_linetype_manual(values=c(1, 3:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changepi0_power_randomeff_plot


# 3d change pi0
setwd(outputDir)
ind3d_changepi0_randomeff <- fread("ind3d_changepi0_randomeff.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  filter(Method != "DACT") %>%
  filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  

# plot 3d fdp
ind3d_changepi0_fdp_randomeff_plot <- ggplot(data=ind3d_changepi0_randomeff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  ylab("FDP (3D Pleiotropy)") +
  ylim(c(0, 0.3)) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))
ind3d_changepi0_fdp_randomeff_plot



# save
changepi0_randomeff_plot <- plot_grid(ind2d_changepi0_noalt_fdp_randomeff_plot + theme(legend.position = "none"),
                                      ind2d_changepi0_fdp_randomeff_plot + theme(legend.position = "none"),
                                  ind3d_changepi0_fdp_randomeff_plot + theme(legend.position = "none"),
                                  ind2d_changepi0_power_randomeff_plot + theme(legend.position = "none"),
                            labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
changepi0_randomeff_legend <- get_legend(ind2d_changepi0_noalt_fdp_randomeff_plot +  theme(legend.direction="horizontal",
                                                                   legend.justification="center",
                                                                   legend.box.just="bottom"))
plot_grid(changepi0_randomeff_plot, changepi0_randomeff_legend, ncol=1, rel_heights=c(1, 0.15))
ggsave('SFig8.pdf', width=20, height=12)



