# Collect results and plot Supp Figure 28

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig28/plot_sfig28.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig28/plot_sfig28.R")

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
outputDir <- here::here("SuppFig28", "output/")
names28f1 <- paste0("SFig28A_", 1:500, "_fitlow.txt")
names28f2 <- paste0("SFig28A_", 1:500, "_fitmed.txt")
names28f3 <- paste0("SFig28A_", 1:500, "_fithigh.txt")

# read raw output files
res28f1 <- c()
for (file_it in 1:length(names28f1)) {
  tempRes <- tryCatch(fread(names28f1[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res28f1, tempRes)
  res28f1 <- rbindlist(tempList)
}

# Read 28f2
res28f2 <- c()
for (file_it in 1:length(names28f2)) {
  tempRes <- tryCatch(fread(names28f2[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res28f2, tempRes)
  res28f2 <- rbindlist(tempList)
}

# Read 28f3
res28f3 <- c()
for (file_it in 1:length(names28f3)) {
  tempRes <- tryCatch(fread(names28f3[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res28f3, tempRes)
  res28f3 <- rbindlist(tempList)
}

# summarize
summary28f1 <- summarize_raw(res28f1)
summary28f2 <- summarize_raw(res28f2)
summary28f3 <- summarize_raw(res28f3)

# save summaries
write.table(summary28f1, paste0(outputDir, "med2d_changeeff_true3_use1_low.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary28f2, paste0(outputDir, "med2d_changeeff_true3_use1_orig.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary28f3, paste0(outputDir, "med2d_changeeff_true3_use1_high.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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
med2d_changeeff_true3_use1_high <- fread(paste0(outputDir, "med2d_changeeff_true3_use1_high.txt"), data.table=F) %>%
  filter(Method == "New") %>%
  mutate(Method = "csmGmm-1-high") %>%
  select(minEff1, Method, Power, FDP)
med2d_changeeff_true3_use1_low <- fread(paste0(outputDir, "med2d_changeeff_true3_use1_low.txt"), data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-low") %>%
  select(minEff1, Method, Power, FDP)
med2d_changeeff_true3_use1_orig <- fread(paste0(outputDir, "med2d_changeeff_true3_use1_orig.txt"), data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-med") %>%
  select(minEff1, Method, Power, FDP)
med2d_changeeff_true3_use1_standard <- fread(paste0(here::here("SuppFig24", "output"), "med2d_changeeff_mbltrue3_use1_2pct_double_init.txt")) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-Standard") %>%
  select(minEff1, Method, Power, FDP)

med2d_changeeff_true3_OM <- rbind(med2d_changeeff_true3_use1_high,
                                  med2d_changeeff_true3_use1_low,
                                  med2d_changeeff_true3_use1_orig,
                                  med2d_changeeff_true3_use1_standard)

# plot med
ind2d_changeeff_true3_OM_plot <- ggplot(data=med2d_changeeff_true3_OM, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (Mediation, 3 Signals)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.2)) + 
  xlim(c(0, 0.15)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=c(4,3,2,1)) +
  scale_color_manual(values=c("purple", "orange", "darkgreen", mycols[1])) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_true3_OM_plot

# save
#ggsave("med2d_changeeff_onemix.pdf", width=12, height=6)

















