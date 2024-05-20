# Plot supp fig 9
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggformula)
library(devtools)
devtools::install_github("ryanrsun/csmGmm")
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#-----------------------------------------#
# change to where the output files are stored
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig9/origOutput"
names9a7 <- paste0("sim_n1k_j100k_ind5d_changeeff2_7df_aID", 1:2000, ".txt")
names9a50 <- paste0("sim_n1k_j100k_ind5d_changeeff2_50df_aID", 1:2000, ".txt")
names9ak <- paste0("sim_n1k_j100k_ind5d_changeeff2_kernel_aID", 1:2000, ".txt")
names9ac <- paste0("sim_n1k_j100k_ind5d_changeeff2_new_aID", 1:2000, ".txt")
names9b <- paste0("sim_n1k_j100k_ind5d_changepi0_combined_aID", 1:800, ".txt")
names9c <- paste0("sim_n1k_j100k_ind6d_changeeff2_combined_aID", 1:2000, ".txt")
names9d <- paste0("sim_n1k_j100k_ind6d_changepi0_combined_aID", 1:800, ".txt")
#-----------------------------------------#


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"

# start reading for s fig 9a
setwd(outputDir)
fullResNew <- c()
fullRes7df <- c()
fullRes50df <- c()
fullResKernel <- c()

for (file_it in 1:length(names9a7)) {
  tempResNew <- fread(names9ac[file_it])
  tempListNew <- list(fullResNew, tempResNew)
  fullResNew <- rbindlist(tempListNew)
    
  tempRes7df <- fread(names9a7[file_it])
  tempList7df <- list(fullRes7df, tempRes7df)
  fullRes7df <- rbindlist(tempList7df)

  tempRes50df <- fread(names9a50[file_it])
  tempList50df <- list(fullRes50df, tempRes50df)
  fullRes50df <- rbindlist(tempList50df)

  tempResKernel <- fread(names9ak[file_it])
  tempListKernel <- list(fullResKernel, tempResKernel)
  fullResKernel <- rbindlist(tempListKernel)
}

outDF <- c()
effSizes <- sort(unique(fullResNew$minEff1))
for (eff_it in 1:length(effSizes)) {
  tempEff <- effSizes[eff_it]
  # weird variable naming so we don't have to change a bunch of code 
  allResNew <- fullResNew %>% filter(minEff1 == tempEff)  %>%
    as.data.frame(.) %>%
    mutate(fdpNew = ifelse(nRejNew == 0, 0, fdpNew))
  allRes7df <- fullRes7df %>% filter(minEff1 == tempEff) %>%
    as.data.frame(.) %>%
    mutate(fdp7df = ifelse(nRej7df == 0, 0, fdp7df))
  allRes50df <- fullRes50df %>% filter(minEff1 == tempEff) %>%
    as.data.frame(.) %>%
    mutate(fdp50df = ifelse(nRej50df == 0, 0, fdp50df))
  allResKernel <- fullResKernel %>% filter(minEff1 == tempEff) %>%
    as.data.frame(.) %>%
    mutate(fdpKernel = ifelse(nRejKernel == 0, 0, fdpKernel))

  summaryOut <- data.frame(minEff1 = allResNew$minEff1[1], minEff2 = allResNew$minEff2[1],
                           Method = c("Kernel", "df7", "df50", "New"))
  summaryOut$nRej <- c(mean(allResKernel$nRejKernel, na.rm=T), mean(allRes7df$nRej7df, na.rm=T), mean(allRes50df$nRej50df, na.rm=T), mean(allResNew$nRejNew, na.rm=T))
  summaryOut$Power <- c(mean(allResKernel$powerKernel, na.rm=T),
                      mean(allRes7df$power7df, na.rm=T), mean(allRes50df$power50df, na.rm=T), mean(allResNew$powerNew, na.rm=T))
  summaryOut$FDP <- c(mean(allResKernel$fdpKernel, na.rm=T), 
                    mean(allRes7df$fdp7df, na.rm=T), mean(allRes50df$fdp50df, na.rm=T), mean(allResNew$fdpNew, na.rm=T))
  summaryOut$pi0a <- c(mean(allResKernel$pi0aKernel, na.rm=T),
                      mean(allRes7df$pi0aDf7, na.rm=T), mean(allRes50df$pi0aDf50, na.rm=T), mean(allResNew$pi0aNew, na.rm=T))
  summaryOut$pi0b <- c(mean(allResKernel$pi0bKernel, na.rm=T),
                                             mean(allRes7df$pi0bDf7, na.rm=T), mean(allRes50df$pi0bDf50, na.rm=T), mean(allResNew$pi0bNew, na.rm=T))                                  
  summaryOut$pi0aTrue <- rep(1 - mean(allResNew$pi0aTrue, na.rm=T) / (10^5), 4)
  summaryOut$pi0bTrue <- rep(1 - mean(allResNew$pi0bTrue, na.rm=T) / (10^5), 4)
  summaryOut$numNA <- c(length(which(is.na(allResKernel$powerKernel))),
                      length(which(is.na(allRes7df$power7df))), length(which(is.na(allRes50df$power50df))), length(which(is.na(allResNew$powerNew))))          

  outDF <- rbind(outDF, summaryOut)
} 

# save results 9a
setwd(outputDir)
write.table(outDF, "med2d_changeeff_correctpi_5d.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

# read others
res9b <- c()
for (file_it in 1:length(names9b)) {
  tempRes <- tryCatch(fread(names9b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res9b, tempRes)
  res9b <- rbindlist(tempList)
}

# Read 9c
res9c <- c()
for (file_it in 1:length(names9c)) {
  tempRes <- tryCatch(fread(names9c[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res9c, tempRes)
  res9c <- rbindlist(tempList)
}

# Read 9d
res9d <- c()
for (file_it in 1:length(names9d)) {
  tempRes <- tryCatch(fread(names9d[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res9d, tempRes)
  res9d <- rbindlist(tempList)
}

# summarize
summary9b <- summarize_raw(res9b)
summary9c <- summarize_raw(res9c)
summary9d <- summarize_raw(res9d)

# save
setwd(outputDir)
write.table(summary9b, "med2d_changepi0_5d.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary9c, "med2d_changeeff_ind6d_v2.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary9d, "ind6d_changepi0.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

#-------------------------------------------------------#
# start plotting

# s fig 9a
setwd(outputDir)
ind5d_changeeff <- fread("med2d_changeeff_correctpi_5d.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 <= 1)

# plot
ind5d_changeeff_fdp_plot <- ggplot(data=ind5d_changeeff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (K=5)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[c(1, 3:6)]) +
  scale_linetype_manual(values=c(1, 3:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind5d_changeeff_fdp_plot

# s fig 9b
setwd(outputDir)
ind5d_changepi0 <- fread("med2d_changepi0_5d.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 <= 1)

# plot
ind5d_changepi0_fdp_plot <- ggplot(data=ind5d_changepi0, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression("Proportion of Variants with One Association")) +
  ylab("FDP (K=5)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind5d_changepi0_fdp_plot


# s fig 9c
setwd(outputDir)
ind6d_changeeff <- fread("med2d_changeeff_ind6d_v2.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 <= 1) %>%
  filter(!(Method %in% c("DACT", "HDMT")))


# plot 
ind6d_changeeff_fdp_plot <- ggplot(data=ind6d_changeeff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (K=6)") +
  ylim(c(0, 0.4)) +
  #xlim(c(0.25, 0.35)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[c(1, 3:6)]) +
  scale_linetype_manual(values=c(1, 3:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind6d_changeeff_fdp_plot

# s fig 9d
setwd(outputDir)
ind6d_changepi0 <- fread("ind6d_changepi0.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 <= 1)

# plot
ind6d_changepi0_fdp_plot <- ggplot(data=ind6d_changepi0, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Proportion of Variants with One Association") +
  ylab("FDP (K=6)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind6d_changepi0_fdp_plot


# save 
ind56d_fdp_plot <- plot_grid(ind5d_changeeff_fdp_plot + theme(legend.position = "none"),
                             ind5d_changepi0_fdp_plot + theme(legend.position = "none"),
                             ind6d_changeeff_fdp_plot + theme(legend.position = "none"),
                             ind6d_changepi0_fdp_plot + theme(legend.position = "none"),
                                 labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
ind56d_fdp_legend <- get_legend(ind5d_changeeff_fdp_plot +  theme(legend.direction="horizontal",
                                                                 legend.justification="center",
                                                                 legend.box.just="bottom"))
plot_grid(ind56d_fdp_plot, ind56d_fdp_legend, ncol=1, rel_heights=c(1, 0.15))
setwd(outputDir)
ggsave('ind56d_fdp.pdf', width=20, height=12)









