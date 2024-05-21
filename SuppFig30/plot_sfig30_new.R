# Collect results and plot Supp Figure 30
library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
setwd('/rsrch3/home/biostatistics/rsun3/github/ancillaryFunctions')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#-----------------------------------------#
# change to where the output files are stored
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig29/output"
names29f1 <- paste0("SFig29E_aID", 1:1000, "_fithigh.txt")
names29f3 <- paste0("SFig29E_aID", 1:1000, "_fitlow.txt")
#-----------------------------------------#

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

# read data
setwd(outputDir)
med2d_mbl5_use1_high <- fread("med_mbl5_onemix_high.txt", data.table=F) %>%
  filter(Method == "New") %>%
  mutate(Method = "csmGmm-1-high") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_use1_low <- fread("med_mbl5_onemix_low.txt", data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-low") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_use1_orig <- fread("med2d_changeeff_mbltrue5_use1_2pct_init.txt", data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-med") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_standard <- fread("med2d_changeeff_mbltrue5_use1_2pct_double_init.txt", data.table=F) %>%
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
setwd(outputDir)
ggsave("med2d_mbl5_onemix_power.pdf", width=12, height=6)











