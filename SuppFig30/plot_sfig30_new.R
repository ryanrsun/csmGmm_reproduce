# Collect results and plot Supp Figure 30

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig30/plot_sfig30_new.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig30/plot_sfig30_new.R")

library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("SuppFig29", "output/")

# plotting starts

# define colors
gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"

# read data
med2d_mbl5_use1_high <- fread(paste0(outputDir, "med_mbl5_onemix_high.txt"), data.table=F) %>%
  filter(Method == "New") %>%
  mutate(Method = "csmGmm-1-high") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_use1_low <- fread(paste0(outputDir, "med_mbl5_onemix_low.txt"), data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-low") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_use1_orig <- fread(here::here("SuppFig26/output/med2d_changeeff_mbltrue5_use1_2pct_init.txt"), data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-med") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_standard <- fread(here::here("SuppFig26/output/med2d_changeeff_mbltrue5_use1_2pct_double_init.txt"), data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-standard") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_OM <- rbind(med2d_mbl5_use1_high,
                            med2d_mbl5_use1_low,
                            med2d_mbl5_use1_orig,
                            med2d_mbl5_standard)

# plot 
med2d_mbl5_OM_power_plot <- ggplot(data=med2d_mbl5_OM, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (Mediation, 5 Signals)") +
  xlab("Min Effect Magnitude") +
  xlim(c(0, 0.15)) + 
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=c(4,3,2,1)) +
  scale_color_manual(values=c("purple", "orange", "darkgreen", mycols[1])) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
med2d_mbl5_OM_power_plot

# save
#ggsave("med2d_mbl5_onemix_power.pdf", width=12, height=6)











