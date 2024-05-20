# Collect results and plot Supp Figure 26
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig24/origOutput"
names26f1 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue5_use1_2pct_init_aID", 1:300, ".txt")
names26f2 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue5_use3_2pct_init_aID", 1:300, ".txt")
names26f3 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue5_use5_2pct_init_aID", 1:300, ".txt")
names26f4 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue5_use7_2pct_init_aID", 1:300, ".txt")
names26fr <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue5_use1_2pct_double_init_aID", 1:300, ".txt")

#-----------------------------------------#

# read raw output files
setwd(outputDir)
res26f1 <- c()
for (file_it in 1:length(names26f1)) {
  tempRes <- tryCatch(fread(names26f1[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res26f1, tempRes)
  res26f1 <- rbindlist(tempList)
}

# Read 26f2
res26f2 <- c()
for (file_it in 1:length(names26f2)) {
  tempRes <- tryCatch(fread(names26f2[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res26f2, tempRes)
  res26f2 <- rbindlist(tempList)
}

# Read 26f3
res26f3 <- c()
for (file_it in 1:length(names26f3)) {
  tempRes <- tryCatch(fread(names26f3[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res26f3, tempRes)
  res26f3 <- rbindlist(tempList)
}

# Read 26f4
res26f4 <- c()
for (file_it in 1:length(names26f4)) {
  tempRes <- tryCatch(fread(names26f4[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res26f4, tempRes)
  res26f4 <- rbindlist(tempList)
}

# Read 26fr
res26fr <- c()
for (file_it in 1:length(names26fr)) {
  tempRes <- tryCatch(fread(names26fr[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res26fr, tempRes)
  res26fr <- rbindlist(tempList)
}


# summarize
summary26f1 <- summarize_raw(res26f1)
summary26f2 <- summarize_raw(res26f2)
summary26f3 <- summarize_raw(res26f3)
summary26f4 <- summarize_raw(res26f4)
summary26fr <- summarize_raw(res26fr)

# save summaries
setwd(outputDir)
write.table(summary26f1, "med2d_changeeff_mbltrue5_use1_2pct_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary26f2, "med2d_changeeff_mbltrue5_use3_2pct_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary26f3, "med2d_changeeff_mbltrue5_use5_2pct_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary26f4, "med2d_changeeff_mbltrue5_use7_2pct_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary26fr, "med2d_changeeff_mbltrue5_use1_2pct_double_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')


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


# read files
setwd(outputDir)
true5_use1_2pct <- fread("med2d_changeeff_mbltrue5_use1_2pct_init.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New1")
true5_double_2pct <- fread("med2d_changeeff_mbltrue5_use1_2pct_double_init.txt") 
true5_use3_2pct <- fread("med2d_changeeff_mbltrue5_use3_2pct_init.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New3")
true5_use5_2pct <- fread("med2d_changeeff_mbltrue5_use5_2pct_init.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New5")
true5_use7_2pct <- fread("med2d_changeeff_mbltrue5_use7_2pct_init.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New7")

# put together data
true5_2pct <- rbind(true5_use1_2pct, true5_double_2pct, true5_use3_2pct,
                    true5_use5_2pct, true5_use7_2pct) %>%
  filter(Method != "DACTb") %>%
  #filter(Method != "df7") %>% filter(Method != "DACT") %>% filter(Method != "Kernel") %>% filter(Method != "df50") %>%
  #filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "New1", "csmGmm-1", Method)) %>%
  mutate(Method = ifelse(Method == "New3", "csmGmm-3", Method)) %>%
  mutate(Method = ifelse(Method == "New5", "csmGmm-5", Method)) %>%
  mutate(Method = ifelse(Method == "New7", "csmGmm-7", Method)) %>%
  #mutate(Method = factor(Method, levels=c("csmGmm-3", "DACT", "Kernel", "locfdr7df",
  #                                     "locfdr50df", "HDMT", "csmGmm-1", "csmGmm-2", "csmGmm-4"))) %>%
  filter(minEff1 <= 1)

# only csm
true5_2pct_csm <- true5_2pct %>%
  filter(Method != "locfdr7df") %>% filter(Method != "DACT") %>% filter(Method != "Kernel") %>% filter(Method != "locfdr50df") %>%
  filter(Method != "HDMT") 
# others
true5_2pct_others <- true5_2pct %>%
  filter(Method != "csmGmm-1") %>% filter(Method != "csmGmm-3") %>% 
  filter(Method != "csmGmm-5") %>% filter(Method != "csmGmm-7") %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

#plot sfig 26a
true5_2pct_csm_fdp_plot <- ggplot(data=true5_2pct_csm, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("FDP (5 Means)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=c(mycols[1], "darkgreen", "purple", "cyan", "orange")) +
  scale_linetype_manual(values=c(1:5)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true5_2pct_csm_fdp_plot

# plot s fig 26b
true5_2pct_csm_power_plot <- ggplot(data=true5_2pct_csm, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("Power (5 Means)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=c(mycols[1], "darkgreen", "purple", "cyan", "orange")) +
  scale_linetype_manual(values=c(1:5)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true5_2pct_csm_power_plot

#---------------------------------------------------------------------------------------#

# plot s fig 26c
true5_2pct_others_fdp_plot <- ggplot(data=true5_2pct_others, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("FDP (5 Means)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true5_2pct_others_fdp_plot

# plot s fig 26d
true5_2pct_others_power_plot <- ggplot(data=true5_2pct_others, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("Power (5 Means)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true5_2pct_others_power_plot

#---------------------------------------------------------------------------------------#

# cowplot together
mbl_true5_2pct_plot <- plot_grid(true5_2pct_csm_fdp_plot + theme(legend.position = "none"),
                                 true5_2pct_csm_power_plot + theme(legend.position = "none"),
                                 true5_2pct_others_fdp_plot + theme(legend.position = "none"),
                                 true5_2pct_others_power_plot + theme(legend.position = "none"),
                                 labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
mbl_true5_2pct_legend1 <- get_legend(true5_2pct_csm_fdp_plot +  theme(legend.direction="horizontal",
                                                                      legend.justification="center",
                                                                      legend.box.just="bottom"))
mbl_true5_2pct_legend2 <- get_legend(true5_2pct_others_fdp_plot +  theme(legend.direction="horizontal",
                                                                         legend.justification="center",
                                                                         legend.box.just="bottom"))
plot_grid(mbl_true5_2pct_plot, mbl_true5_2pct_legend1, mbl_true5_2pct_legend2, ncol=1, rel_heights=c(1, 0.1, 0.1))
setwd(outputDir)
ggsave('mbl_true5_2pct.pdf', width=18, height=12)









