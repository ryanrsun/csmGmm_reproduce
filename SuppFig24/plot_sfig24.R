# Collect results and plot Supp Figure 24
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig24/origOutput"
names24f1 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use1_2pct_init_aID", 1:300, ".txt")
names24f2 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use2_2pct_init_aID", 1:300, ".txt")
names24f3 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use3_2pct_init_aID", 1:300, ".txt")
names24f4 <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use4_2pct_init_aID", 1:300, ".txt")
names24fr <- paste0("sim_n1k_j100k_med2d_changeeff_mbltrue3_use1_2pct_double_init_aID", 1:300, ".txt")

#-----------------------------------------#

# read raw output files
setwd(outputDir)
res24f1 <- c()
for (file_it in 1:length(names24f1)) {
  tempRes <- tryCatch(fread(names24f1[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res24f1, tempRes)
  res24f1 <- rbindlist(tempList)
}

# Read 24f2
res24f2 <- c()
for (file_it in 1:length(names24f2)) {
  tempRes <- tryCatch(fread(names24f2[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res24f2, tempRes)
  res24f2 <- rbindlist(tempList)
}

# Read 24f3
res24f3 <- c()
for (file_it in 1:length(names24f3)) {
  tempRes <- tryCatch(fread(names24f3[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res24f3, tempRes)
  res24f3 <- rbindlist(tempList)
}

# Read 24f4
res24f4 <- c()
for (file_it in 1:length(names24f4)) {
  tempRes <- tryCatch(fread(names24f4[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res24f4, tempRes)
  res24f4 <- rbindlist(tempList)
}

# Read 24fr
res24fr <- c()
for (file_it in 1:length(names24fr)) {
  tempRes <- tryCatch(fread(names24fr[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res24fr, tempRes)
  res24fr <- rbindlist(tempList)
}


# summarize
summary24f1 <- summarize_raw(res24f1)
summary24f2 <- summarize_raw(res24f2)
summary24f3 <- summarize_raw(res24f3)
summary24f4 <- summarize_raw(res24f4)
summary24fr <- summarize_raw(res24fr)

# save summaries
setwd(outputDir)
write.table(summary24f1, "med2d_changeeff_mbltrue3_use1_2pct_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary24f2, "med2d_changeeff_mbltrue3_use2_2pct_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary24f3, "med2d_changeeff_mbltrue3_use3_2pct_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary24f4, "med2d_changeeff_mbltrue3_use4_2pct_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary24fr, "med2d_changeeff_mbltrue3_use1_2pct_double_init.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')


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
true3_use1_2pct <- fread("med2d_changeeff_mbltrue3_use1_2pct_init.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New1")
true3_double_2pct <- fread("med2d_changeeff_mbltrue3_use1_2pct_double_init.txt") 
true3_use2_2pct <- fread("med2d_changeeff_mbltrue3_use2_2pct_init.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New2")
true3_use3_2pct <- fread("med2d_changeeff_mbltrue3_use3_2pct_init.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New3")
true3_use4_2pct <- fread("med2d_changeeff_mbltrue3_use4_2pct_init.txt") %>%
  filter(Method == "New") %>%
  mutate(Method = "New4")

# put together data
true3_2pct <- rbind(true3_use1_2pct, true3_double_2pct, true3_use2_2pct,
                    true3_use3_2pct, true3_use4_2pct) %>%
  filter(Method != "DACTb") %>%
  #filter(Method != "df7") %>% filter(Method != "DACT") %>% filter(Method != "Kernel") %>% filter(Method != "df50") %>%
  #filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "New1", "csmGmm-1", Method)) %>%
  mutate(Method = ifelse(Method == "New2", "csmGmm-2", Method)) %>%
  mutate(Method = ifelse(Method == "New3", "csmGmm-3", Method)) %>%
  mutate(Method = ifelse(Method == "New4", "csmGmm-4", Method)) %>%
  #mutate(Method = factor(Method, levels=c("csmGmm-3", "DACT", "Kernel", "locfdr7df",
  #                                     "locfdr50df", "HDMT", "csmGmm-1", "csmGmm-2", "csmGmm-4"))) %>%
  filter(minEff1 <= 1)

# only csm
true3_2pct_csm <- true3_2pct %>%
  filter(Method != "locfdr7df") %>% filter(Method != "DACT") %>% filter(Method != "Kernel") %>% filter(Method != "locfdr50df") %>%
  filter(Method != "HDMT") 
# others
true3_2pct_others <- true3_2pct %>%
  filter(Method != "csmGmm-1") %>% filter(Method != "csmGmm-2") %>% 
  filter(Method != "csmGmm-3") %>% filter(Method != "csmGmm-4") %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot mbl true 3 change eff fdp csm
true3_2pct_csm_fdp_plot <- ggplot(data=true3_2pct_csm, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("FDP (3 Means)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=c(mycols[1], "darkgreen", "purple", "cyan", "orange")) +
  scale_linetype_manual(values=c(1:5)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true3_2pct_csm_fdp_plot

# plot mbl true 3 change eff power csm
true3_2pct_csm_power_plot <- ggplot(data=true3_2pct_csm, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("Power (3 Means)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=c(mycols[1], "darkgreen", "purple", "cyan", "orange")) +
  scale_linetype_manual(values=c(1:5)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true3_2pct_csm_power_plot

#---------------------------------------------------------------------------------------#

# plot mbl true 3 change eff fdp others
true3_2pct_others_fdp_plot <- ggplot(data=true3_2pct_others, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("FDP (3 Means)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true3_2pct_others_fdp_plot

# plot mbl true 3 change eff power others
true3_2pct_others_power_plot <- ggplot(data=true3_2pct_others, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Effect Magnitude") +
  ylab("Power (3 Means)") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
true3_2pct_others_power_plot

#---------------------------------------------------------------------------------------#

# cowplot all together
mbl_true3_2pct_plot <- plot_grid(true3_2pct_csm_fdp_plot + theme(legend.position = "none"),
                                      true3_2pct_csm_power_plot + theme(legend.position = "none"),
                                      true3_2pct_others_fdp_plot + theme(legend.position = "none"),
                                      true3_2pct_others_power_plot + theme(legend.position = "none"),
                                      labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
mbl_true3_2pct_legend1 <- get_legend(true3_2pct_csm_fdp_plot +  theme(legend.direction="horizontal",
                                                                                legend.justification="center",
                                                                                legend.box.just="bottom"))
mbl_true3_2pct_legend2 <- get_legend(true3_2pct_others_fdp_plot +  theme(legend.direction="horizontal",
                                                                                   legend.justification="center",
                                                                                   legend.box.just="bottom"))
#true3_largeProp_legend <- plot_grid(mbl_true3_largeProp_legend1, mbl_true3_largeProp_legend2, ncol=2)
plot_grid(mbl_true3_2pct_plot, mbl_true3_2pct_legend1, mbl_true3_2pct_legend2, ncol=1, rel_heights=c(1, 0.1, 0.1))
setwd(outputDir)
ggsave('mbl_true3_2pct.pdf', width=18, height=12)












