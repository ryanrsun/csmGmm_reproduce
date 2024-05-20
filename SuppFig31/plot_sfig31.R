# Collect results and plot Supp Figure 31
library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
devtools::install_github("ryanrsun/csmGmm")
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#-----------------------------------------#
# change to where the output files are stored
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig31/origOutput"
names31a <- paste0("sim_n1k_j100k_med2d_changeeff_maxP_aID", 1:400, ".txt")
names31b <- paste0("sim_n1k_j100k_med2d_changepi0_maxP_aID", 1:160, ".txt")
#-----------------------------------------#

# read raw output files
setwd(outputDir)
res31a <- c()
for (file_it in 1:length(names31a)) {
  tempRes <- tryCatch(fread(names31a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res31a, tempRes)
  res31a <- rbindlist(tempList)
}

# Read 31b
res31b <- c()
for (file_it in 1:length(names31b)) {
  tempRes <- tryCatch(fread(names31b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res31b, tempRes)
  res31b <- rbindlist(tempList)
}

# summarize
summary31a <- summarize_raw(res31a, maxP=T)
summary31b <- summarize_raw(res31b, maxP=T)

# save summaries
setwd(outputDir)
write.table(summary31a, "med2d_changeeff_maxP.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary31b, "med2d_changepi0_maxP.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')


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


# maxP change eff
setwd(outputDir)
maxP_med2d_changeeff <- fread("med2d_changeeff_maxP.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "MaxP", "maxP", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT", "maxP"))) %>%
  filter(minEff1 <= 1) %>%
  mutate(FDP = ifelse(is.na(FDP), 0, FDP))

# plot maxP change eff fdp
mycols2 <- c(mycols, "purple")
maxP_med2d_changeeff_fdp_plot <- ggplot(data=maxP_med2d_changeeff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (Mediation)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols2) +
  scale_linetype_manual(values=c(1:7)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
maxP_med2d_changeeff_fdp_plot

# plot maxP change eff power
maxP_med2d_changeeff_power_plot <- ggplot(data=maxP_med2d_changeeff %>% filter(!(Method %in% c("DACT", "locfdr7df"))),
                                          aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("Power (Mediation)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=mycols2[c(1, 3, 5, 6, 7)]) +
  scale_linetype_manual(values=c(1, 3, 5, 6, 7)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
maxP_med2d_changeeff_power_plot


# maxP changepi0
maxP_med2d_changepi0 <- fread("med2d_changepi0_maxP.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "MaxP", "maxP", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT", "maxP"))) %>%
  filter(minEff1 <= 1) %>%
  mutate(FDP = ifelse(is.na(FDP), 0, FDP))

# plot maxP change pi0
maxP_med2d_changepi0_fdp_plot <- ggplot(data=maxP_med2d_changepi0, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (Mediation)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols2) +
  scale_linetype_manual(values=1:7) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
maxP_med2d_changepi0_fdp_plot



# plot maxP change pi0 power
maxP_med2d_changepi0_power_plot <- ggplot(data=maxP_med2d_changepi0 %>% filter(!(Method %in% c("DACT"))),
                                          aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("Power (Mediation)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=mycols2[c(1, 3, 4, 5, 6, 7)]) +
  scale_linetype_manual(values=c(1, 3, 4, 5, 6, 7)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
maxP_med2d_changepi0_power_plot


# save ind 2d
maxP_plot <- plot_grid(maxP_med2d_changeeff_fdp_plot + theme(legend.position = "none"),
                       maxP_med2d_changepi0_fdp_plot + theme(legend.position = "none"),
                       maxP_med2d_changeeff_power_plot + theme(legend.position = "none"),
                       maxP_med2d_changepi0_power_plot + theme(legend.position = "none"),
                                labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
maxP_legend <- get_legend(maxP_med2d_changeeff_fdp_plot +  theme(legend.direction="horizontal",
                                                               legend.justification="center",
                                                               legend.box.just="bottom"))
plot_grid(maxP_plot, maxP_legend, ncol=1, rel_heights=c(1, 0.15))
setwd(outputDir)
ggsave('maxp.pdf', width=20, height=12)












