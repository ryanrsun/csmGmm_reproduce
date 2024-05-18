# Collect results and plot Supp Figure 28
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig28/origOutput"
names28f1 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use1_low_aID", 1:500, ".txt")
names28f2 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use1_orig_aID", 1:500, ".txt")
names28f3 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use1_high_aID", 1:500, ".txt")
names28fr <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use1_2pct_double_init_aID", 1:500, ".txt")
#-----------------------------------------#

# read raw output files
setwd(outputDir)
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

# Read 28fr
res28fr <- c()
for (file_it in 1:length(names28fr)) {
  tempRes <- tryCatch(fread(names28fr[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res28fr, tempRes %>% mutate(nRejDACT=NA, nRejHDMT=NA, nRejKernel=NA, nRej7df=NA, nRej50df=NA))
  res28fr <- rbindlist(tempList)
}

# summarize
summary28f1 <- summarize_raw(res28f1)
summary28f2 <- summarize_raw(res28f2)
summary28f3 <- summarize_raw(res28f3)
summary28fr <- summarize_raw(res28fr)

# save summaries
setwd(outputDir)
write.table(summary28f1, "med2d_changeeff_true3_use1_low.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary28f2, "med2d_changeeff_true3_use1_orig.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary28f3, "med2d_changeeff_true3_use1_high.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary28fr, "med2d_changeeff_mbltrue3_use1_2pct_double_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')


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
med2d_changeeff_true3_use1_high <- fread("med2d_changeeff_true3_use1_high.txt", data.table=F) %>%
  filter(Method == "New") %>%
  mutate(Method = "csmGmm-1-high") %>%
  select(minEff1, Method, Power, FDP)
med2d_changeeff_true3_use1_low <- fread("med2d_changeeff_true3_use1_low.txt", data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-low") %>%
  select(minEff1, Method, Power, FDP)
med2d_changeeff_true3_use1_orig <- fread("med2d_changeeff_true3_use1_orig.txt", data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "csmGmm-1-med") %>%
  select(minEff1, Method, Power, FDP)
med2d_changeeff_true3_use1_standard <- fread("med2d_changeeff_mbltrue3_use1_2pct_double_init.txt") %>%
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
setwd(outputDir)
ggsave("med2d_changeeff_onemix.pdf", width=12, height=6)

















