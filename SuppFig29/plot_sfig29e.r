# Collect results and plot Supp Figure 29

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig29/plot_sfig29e.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig29/plot_sfig29e.R")

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
names29f1 <- paste0(outputDir, "SFig29E_aID", 1:1000, "_fithigh.txt")
names29f3 <- paste0(outputDir, "SFig29E_aID", 1:1000, "_fitlow.txt")
#-----------------------------------------#

# read raw output files
res29f1 <- c()
for (file_it in 1:length(names29f1)) {
  tempRes <- tryCatch(fread(names29f1[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res29f1, tempRes %>% mutate(nRejDACT=0, fdpDACT=0, nRejHDMT=0, fdpHDMT=0, nRejKernel=0, fdpKernel=0,
                                               nRej7df=0, fdp7df=0, nRej50df=0, fdp50df=0))
  res29f1 <- rbindlist(tempList)
}

# Read 29f3
res29f3 <- c()
for (file_it in 1:length(names29f3)) {
  tempRes <- tryCatch(fread(names29f3[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res29f3, tempRes %>% mutate(nRejDACT=0, fdpDACT=0, nRejHDMT=0, fdpHDMT=0, nRejKernel=0, fdpKernel=0,
                                               nRej7df=0, fdp7df=0, nRej50df=0, fdp50df=0))
  res29f3 <- rbindlist(tempList)
}

# summarize
summary29f1 <- summarize_raw(res29f1)
summary29f3 <- summarize_raw(res29f3)

# save summaries
write.table(summary29f1, paste0(outputDir, "med_mbl5_onemix_high.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary29f3, paste0(outputDir, "med_mbl5_onemix_low.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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
med2d_mbl5_use1_high <- fread(paste0(outputDir, "med_mbl5_onemix_high.txt"), data.table=F) %>%
  filter(Method == "New") %>%
  mutate(Method = "csmGmm-1-high") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_use1_low <- fread(paste0(outputDir, "med_mbl5_onemix_low.txt"), data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-low") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_use1_orig <- fread(here::here("SuppFig24/output/med2d_changeeff_mbltrue5_use1_2pct_init.txt"), data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-med") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_standard <- fread(here::here("SuppFig24/output/med2d_changeeff_mbltrue5_use1_2pct_double_init.txt"), data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-standard") %>%
  select(minEff1, Method, Power, FDP)
med2d_mbl5_OM <- rbind(med2d_mbl5_use1_high,
                            med2d_mbl5_use1_low,
                            med2d_mbl5_use1_orig,
                            med2d_mbl5_standard)

# plot 
med2d_mbl5_OM_plot <- ggplot(data=med2d_mbl5_OM, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (Mediation, 5 Signals)") +
  xlab("Min Effect Magnitude") +
  xlim(c(0, 0.15)) + 
  ylim(c(0, 0.2)) + 
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=c(4,3,2,1)) +
  scale_color_manual(values=c("purple", "orange", "darkgreen", mycols[1])) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
med2d_mbl5_OM_plot

# save
ggsave("med2d_mbl5_onemix.pdf", width=12, height=6)











