# Collect results and plot Supp Figure 13

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig13/plot_sfig13_full.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig13/plot_sfig13_full.R")

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
outputDir <- here::here("SuppFig12", "output/")

# first panel
ind2d_changeeff_maf01 <- fread(paste0(outputDir, "med2d_changeeff_maf01.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# read data
ind2d_changepi0_maf01 <- fread(paste0(outputDir, "med2d_changepi0_maf01.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# read data
ind2d_changeeff_maf10 <- fread(paste0(outputDir, "SFig12c_maf10_summary.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# read data
ind2d_changepi0_maf10 <- fread(paste0(outputDir, "SFig12d_maf10_summary.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

# define colors
gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"


# plot ind2d changeeff uncorrected power
ind2d_changeeff_pow_maf01_plot <- ggplot(data=ind2d_changeeff_maf01,
                                           aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") + ylab("Power (2D Mediation)") +
  ylim(c(0, 0.8)) +
  xlim(c(0.4, 0.6)) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_pow_maf01_plot


# plot ind2d changepi0 uncorrected power
ind2d_changepi0_pow_maf01_plot <- ggplot(data=ind2d_changepi0_maf01 %>% filter(Method != "DACT"),
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
ind2d_changepi0_pow_maf01_plot


# ind2d change eff number incongruous
ind2d_changeeff_incon_maf01_plot <- ggplot(data=ind2d_changeeff_maf01, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") +
  ylab("Num Incongruous") +
  xlim(c(0.4, 0.6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_incon_maf01_plot

# ind2d change pi0 number incongruous
ind2d_changepi0_incon_maf01_plot <- ggplot(data=ind2d_changepi0_maf01, aes(x=minEff1, y=Incongruous, group=Method)) +
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
ind2d_changepi0_incon_maf01_plot

#---------------------------------------------------------#
# MAF 10%

# plot ind2d changeeff uncorrected power
ind2d_changeeff_pow_maf10_plot <- ggplot(data=ind2d_changeeff_maf10,
                                           aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") + ylab("Power (2D Mediation)") +
  xlim(c(0.2, 0.4)) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# plot ind2d changepi0 uncorrected power
ind2d_changepi0_pow_maf10_plot <- ggplot(data=ind2d_changepi0_maf10 %>% filter(Method != "DACT"),
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
ind2d_changepi0_pow_maf01_plot


# ind2d change eff number incongruous
ind2d_changeeff_incon_maf10_plot <- ggplot(data=ind2d_changeeff_maf10, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") +
  ylab("Num Incongruous") +
  xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

# ind2d change pi0 number incongruous
ind2d_changepi0_incon_maf10_plot <- ggplot(data=ind2d_changepi0_maf10, aes(x=minEff1, y=Incongruous, group=Method)) +
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


# put together figure 
ind2d_power_incon_maf01_plot <- plot_grid(ind2d_changeeff_pow_maf01_plot + theme(legend.position = "none"),
                                    ind2d_changepi0_pow_maf01_plot + theme(legend.position = "none"),
                                    ind2d_changeeff_incon_maf01_plot + theme(legend.position = "none"),
                                    ind2d_changepi0_incon_maf01_plot + theme(legend.position = "none"),
                                    ind2d_changeeff_pow_maf10_plot + theme(legend.position = "none"),
                                    ind2d_changepi0_pow_maf10_plot + theme(legend.position = "none"),
                                    ind2d_changeeff_incon_maf10_plot + theme(legend.position = "none"),
                                    ind2d_changepi0_incon_maf10_plot + theme(legend.position = "none"),
                                    labels=c("A", "B", "C", "D", "E", "F", "G", "H"), ncol=2, label_size=22)
ind_changeeff_power_incon_maf01_legend <- get_legend(ind2d_changeeff_pow_maf01_plot +  theme(legend.direction="horizontal",
                                                                                  legend.justification="center",
                                                                                  legend.box.just="bottom"))
plot_grid(ind2d_power_incon_maf01_plot, ind_changeeff_power_incon_maf01_legend, ncol=1, rel_heights=c(1, 0.1))

ggsave("sfig13_full.pdf", width=20, height=24)



