# Collect results and plot Supp Figure 25
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig25/origOutput"
names25f1 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use1_changeeff_aID", 1:500, ".txt")
names25f2 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use2_changeeff_aID", 1:500, ".txt")
names25f3 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use3_changeeff_aID", 1:500, ".txt")
names25f4 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use4_changeeff_aID", 1:500, ".txt")
names25fr <- paste0("SFig25A_aID", 1:500, "_fitreg.txt")

#-----------------------------------------#

# read raw output files
setwd(outputDir)
res25f1 <- c()
for (file_it in 1:length(names25f1)) {
  tempRes <- tryCatch(fread(names25f1[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res25f1, tempRes)
  res25f1 <- rbindlist(tempList)
}

# Read 25f2
res25f2 <- c()
for (file_it in 1:length(names25f2)) {
  tempRes <- tryCatch(fread(names25f2[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res25f2, tempRes)
  res25f2 <- rbindlist(tempList)
}

# Read 25f3
res25f3 <- c()
for (file_it in 1:length(names25f3)) {
  tempRes <- tryCatch(fread(names25f3[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res25f3, tempRes)
  res25f3 <- rbindlist(tempList)
}

# Read 25f4
res25f4 <- c()
for (file_it in 1:length(names25f4)) {
  tempRes <- tryCatch(fread(names25f4[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res25f4, tempRes)
  res25f4 <- rbindlist(tempList)
}

# Read 25fr
res25fr <- c()
for (file_it in 1:length(names25fr)) {
  tempRes <- tryCatch(fread(names25fr[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res25fr, tempRes)
  res25fr <- rbindlist(tempList)
}

# summarize
summary25f1 <- summarize_raw(res25f1)
summary25f2 <- summarize_raw(res25f2)
summary25f3 <- summarize_raw(res25f3)
summary25f4 <- summarize_raw(res25f4)
summary25fr <- summarize_raw(res25fr)

# save summaries
setwd(outputDir)
write.table(summary25f1, "med2d_changeeff_mbltrue3_use1_largeProp.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary25f2, "med2d_changeeff_mbltrue3_use2_largeProp.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary25f3, "med2d_changeeff_mbltrue3_use3_largeProp.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary25f4, "med2d_changeeff_mbltrue3_use4_largeProp.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary25fr, "med2d_changeeff_mbltrue3_use1_double_largeProp.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')


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
true3_use1_largeProp <- fread("med2d_changeeff_mbltrue3_use1_largeProp.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New1")
true3_double_largeProp <- fread("med2d_changeeff_mbltrue3_use1_double_largeProp.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "NewD")
true3_use2_largeProp <- fread("med2d_changeeff_mbltrue3_use2_largeProp.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New2")
true3_use3_largeProp <- fread("med2d_changeeff_mbltrue3_use3_largeProp.txt")
true3_use4_largeProp <- fread("med2d_changeeff_mbltrue3_use4_largeProp.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New4")

# put together data
true3_largeProp <- rbind(true3_use1_largeProp, true3_double_largeProp, true3_use2_largeProp,
                         true3_use3_largeProp, true3_use4_largeProp) %>%
  filter(Method != "DACTb") %>%
  #filter(Method != "df7") %>% filter(Method != "DACT") %>% filter(Method != "Kernel") %>% filter(Method != "df50") %>%
  #filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm-3", Method)) %>%
  mutate(Method = ifelse(Method == "NewD", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "New1", "csmGmm-1", Method)) %>%
  mutate(Method = ifelse(Method == "New2", "csmGmm-2", Method)) %>%
  mutate(Method = ifelse(Method == "New4", "csmGmm-4", Method)) %>%
  #mutate(Method = factor(Method, levels=c("csmGmm-3", "DACT", "Kernel", "locfdr7df",
     #                                     "locfdr50df", "HDMT", "csmGmm-1", "csmGmm-2", "csmGmm-4"))) %>%
  filter(minEff1 <= 1)

# only csm
true3_largeProp_csm <- true3_largeProp %>%
  filter(Method != "locfdr7df") %>% filter(Method != "DACT") %>% filter(Method != "Kernel") %>% filter(Method != "locfdr50df") %>%
  filter(Method != "HDMT") 
# others
true3_largeProp_others <- true3_largeProp %>%
  filter(Method != "csmGmm-1") %>% filter(Method != "csmGmm-2") %>% filter(Method != "csmGmm-3") %>%
  filter(Method != "csmGmm-4") %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                      "locfdr50df", "HDMT")))

# plot sfig 25 a
true3_largeProp_csm_fdp_plot <- ggplot(data=true3_largeProp_csm, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("FDP (3 Means, Many Alt)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=c(mycols[1], "darkgreen", "purple", "cyan", "orange")) +
  scale_linetype_manual(values=c(1:5)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true3_largeProp_csm_fdp_plot

# s fig 25 b
true3_largeProp_csm_power_plot <- ggplot(data=true3_largeProp_csm, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("Power (3 Means, Many Alt)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=c(mycols[1], "darkgreen", "purple", "cyan", "orange")) +
  scale_linetype_manual(values=c(1:5)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true3_largeProp_csm_power_plot


#----------------------------------------------------------#

# s fig 25 c
true3_largeProp_others_fdp_plot <- ggplot(data=true3_largeProp_others, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("FDP (3 Means, Many Alt)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true3_largeProp_others_fdp_plot

# s fig 25 d
true3_largeProp_others_power_plot <- ggplot(data=true3_largeProp_others, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("Power (3 Means, Many Alt)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true3_largeProp_others_power_plot

#----------------------------------------------------#
# put it together

# cowplot the large prop
mbl_true3_largeProp_plot <- plot_grid(true3_largeProp_csm_fdp_plot + theme(legend.position = "none"),
                           true3_largeProp_csm_power_plot + theme(legend.position = "none"),
                           true3_largeProp_others_fdp_plot + theme(legend.position = "none"),
                           true3_largeProp_others_power_plot + theme(legend.position = "none"),
                                    labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
mbl_true3_largeProp_legend1 <- get_legend(true3_largeProp_csm_fdp_plot +  theme(legend.direction="horizontal",
                                                                                  legend.justification="center",
                                                                                  legend.box.just="bottom"))
mbl_true3_largeProp_legend2 <- get_legend(true3_largeProp_others_fdp_plot +  theme(legend.direction="horizontal",
                                                                                legend.justification="center",
                                                                                legend.box.just="bottom"))
plot_grid(mbl_true3_largeProp_plot, mbl_true3_largeProp_legend1, mbl_true3_largeProp_legend2, ncol=1, rel_heights=c(1, 0.1, 0.1))
setwd(outputDir)
ggsave('mbl_true3_largeprop.pdf', width=18, height=12)











