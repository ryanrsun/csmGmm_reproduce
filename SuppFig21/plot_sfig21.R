# Collect results and plot Supp Figure 21

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig21/plot_sfig21.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig21/plot_sfig21.R")

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
outputDir <- here::here("SuppFig21", "output/")
names21a <- paste0(outputDir, "SFig21B_aID", 1:400, ".txt")
names21b <- paste0(outputDir, "SFig21B_aID", 1:160, ".txt")
names21c <- paste0(outputDir, "SFig21B_aID", 1:400, ".txt")
names21d <- paste0(outputDir, "SFig21B_aID", 1:160, ".txt")
#-----------------------------------------#

# read raw output files
res21a <- c()
for (file_it in 1:length(names21a)) {
  tempRes <- tryCatch(fread(names21a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res21a, tempRes)
  res21a <- rbindlist(tempList)
}

# Read 21b
res21b <- c()
for (file_it in 1:length(names21b)) {
  tempRes <- tryCatch(fread(names21b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res21b, tempRes)
  res21b <- rbindlist(tempList)
}

# Read 21c
res21c <- c()
for (file_it in 1:length(names21c)) {
  tempRes <- tryCatch(fread(names21c[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res21c, tempRes)
  res21c <- rbindlist(tempList)
}

# Read 21d
res21d <- c()
for (file_it in 1:length(names21d)) {
  tempRes <- tryCatch(fread(names21d[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res21d, tempRes)
  res21d <- rbindlist(tempList)
}

# summarize
summary21a <- summarize_raw(res21a)
summary21b <- summarize_raw(res21b)
summary21c <- summarize_raw(res21c)
summary21d <- summarize_raw(res21d)

# save summaries
write.table(summary21a, paste0(outputDir, "med2d_changeeff_noalt_noAssumption.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary21b, paste0(outputDir, "med2d_changepi0_noalt_noAssumption.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary21c, paste0(outputDir, "med2d_changeeff_noAssumption.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary21d, paste0(outputDir, "med2d_changepi0_noAssumption.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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

# SFig 21a
# no effect change eff
ind2d_noalt_changeeff_noAss <- fread(paste0(outputDir, "med2d_changeeff_noalt_noAssumption.txt"), data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  

# plot 21a
ind2d_changeeff_noalt_fdp_noAss_plot <- ggplot(data=ind2d_noalt_changeeff_noAss, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (Mediation, No Alt)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.2)) + 
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_noalt_fdp_noAss_plot


# S Fig 21b
ind2d_noalt_changepi0_noAss <- fread(paste0(outputDir, "med2d_changepi0_noalt_noAssumption.txt"), data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot 21b
ind2d_changepi0_noalt_fdp_noAss_plot <-ggplot(data=ind2d_noalt_changepi0_noAss, aes(x=minEff1, y=FDP, group=Method)) +
  #geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
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
ind2d_changepi0_noalt_fdp_noAss_plot


# S Fig 21 c
ind2d_changeeff_noAss <- fread(paste0(outputDir, "med2d_changeeff_noAssumption.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) 

# plot 21c
ind2d_changeeff_fdp_noAss_plot <- ggplot(data=ind2d_changeeff_noAss, aes(x=minEff1, y=FDP, group=Method)) +
  #geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (Mediation, With Alt)") +
  ylim(c(0, 0.2)) +
  #xlim(c(0.4, 0.6)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_fdp_noAss_plot


# S Fig 21 d
ind2d_changepi0_noAss <- fread(paste0(outputDir, "med2d_changepi0_noAssumption.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# plot 1d
ind2d_changepi0_fdp_noAss_plot <- ggplot(data=ind2d_changepi0_noAss, aes(x=minEff1, y=FDP, group=Method)) +
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
ind2d_changepi0_fdp_noAss_plot


# save
ind2d_fdp_noAss_plot <- plot_grid(ind2d_changeeff_noalt_fdp_noAss_plot + theme(legend.position = "none"),
                                  ind2d_changepi0_noalt_fdp_noAss_plot + theme(legend.position = "none"),
                                  ind2d_changeeff_fdp_noAss_plot + theme(legend.position = "none"),
                                  ind2d_changepi0_fdp_noAss_plot + theme(legend.position = "none"),
                                  labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
ind2d_fdp_noAss_legend <- get_legend(ind2d_changeeff_noalt_fdp_noAss_plot +  theme(legend.direction="horizontal",
                                                                               legend.justification="center",
                                                                               legend.box.just="bottom"))
plot_grid(ind2d_fdp_noAss_plot, ind2d_fdp_noAss_legend, ncol=1, rel_heights=c(1, 0.15))
#ggsave('ind2d_noalt_noAssumption.pdf', width=20, height=12)


#--------------------------------------------------------------------------------------------#
# SFig 22

# S Fig 22a
ind2d_changeeff_pow_noAss_plot <- ggplot(data=ind2d_changeeff_noAss,
                                         aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") + ylab("Power (2D Mediation)") +
  ylim(c(0, 0.8)) +
  #xlim(c(0.4, 0.6)) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_pow_noAss_plot


# S fig 22b
ind2d_changepi0_pow_noAss_plot <- ggplot(data=ind2d_changepi0_noAss %>% filter(Method != "DACT"),
                                         aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("Power (2D Mediation)") +
  scale_color_manual(values=mycols[c(1, 3:6)]) +
  scale_linetype_manual(values=c(1,3:6)) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))
ind2d_changepi0_pow_noAss_plot


# s fig 22 c
ind2d_changeeff_incon_noAss_plot <- ggplot(data=ind2d_changeeff_noAss, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") +
  ylab("Num Incongruous") +
  #xlim(c(0.4, 0.6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_incon_noAss_plot

# s fig 22 d
ind2d_changepi0_incon_noAss_plot <- ggplot(data=ind2d_changepi0_noAss, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("Num Incongruous") +
  #xlim(c(0, 0.08)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changepi0_incon_noAss_plot

# put together supp figure 22
ind2d_power_incon_noAss_plot <- plot_grid(ind2d_changeeff_pow_noAss_plot + theme(legend.position = "none"),
                                          ind2d_changepi0_pow_noAss_plot + theme(legend.position = "none"),
                                          ind2d_changeeff_incon_noAss_plot + theme(legend.position = "none"),
                                          ind2d_changepi0_incon_noAss_plot + theme(legend.position = "none"),
                                          labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
ind_changeeff_power_incon_noAss_legend <- get_legend(ind2d_changeeff_pow_noAss_plot +  theme(legend.direction="horizontal",
                                                                                             legend.justification="center",
                                                                                             legend.box.just="bottom"))
plot_grid(ind2d_power_incon_noAss_plot, ind_changeeff_power_incon_noAss_legend, ncol=1, rel_heights=c(1, 0.1))
#ggsave('ind2d_power_incon_noAss.pdf', width=18, height=12)


