# Collect results and plot Supp Figure 31

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig31/plot_sfig31.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig31/plot_sfig31.R")

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
outputDir <- here::here("SuppFig14", "output/")
names31a <- paste0(outputDir, "SFig31A_aID", 1:400, ".txt")
names31b <- paste0(outputDir, "SFig31B_aID", 1:160, ".txt")
#-----------------------------------------#

# read raw output files
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
write.table(summary31a, paste0(outputDir, "med2d_changeeff_maxP.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary31b, paste0(outputDir, "med2d_changepi0_maxP.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')


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
maxP_med2d_changeeff <- fread(paste0(outputDir, "med2d_changeeff_maxP.txt"), data.table=F) %>%
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
maxP_med2d_changepi0 <- fread(paste0(outputDir, "med2d_changepi0_maxP.txt"), data.table=F) %>%
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
#ggsave('maxp.pdf', width=20, height=12)












