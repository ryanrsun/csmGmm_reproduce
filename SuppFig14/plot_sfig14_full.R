# Collect results and plot Supp Figure 14

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig14/plot_sfig14_full.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig14/plot_sfig14_full.R")

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
outputDir <- here::here("SuppFig14", "output/")
names14a <- paste0(outputDir, "SFig14A_aID", 1:2500, ".txt")
names14b <- paste0(outputDir, "SFig14B_aID", 1:800, ".txt")
names14c <- paste0(outputDir, "SFig14C_aID", 1:2000, ".txt")
names14d <- paste0(outputDir, "SFig14D_aID", 1:1000, ".txt")

names14a_maf10 <- paste0(outputDir, "SFig14A1_aID", 1:2500, ".txt")
names14b_maf10 <- paste0(outputDir, "SFig14B1_aID", 1:800, ".txt")
names14c_maf10 <- paste0(outputDir, "SFig14C1_aID", 1:2000, ".txt")
names14d_maf10 <- paste0(outputDir, "SFig14D1_aID", 1:1000, ".txt")
#-----------------------------------------#


# read raw output files
res14a <- c()
for (file_it in 1:length(names14a)) {
  tempRes <- tryCatch(fread(names14a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res14a, tempRes)
  res14a <- rbindlist(tempList)
}

# Read 14b
res14b <- c()
for (file_it in 1:length(names14b)) {
  tempRes <- tryCatch(fread(names14b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res14b, tempRes)
  res14b <- rbindlist(tempList)
}

# Read 14c
res14c <- c()
for (file_it in 1:length(names14c)) {
  tempRes <- tryCatch(fread(names14c[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res14c, tempRes)
  res14c <- rbindlist(tempList)
}

# Read 14d
res14d <- c()
for (file_it in 1:length(names14d)) {
  tempRes <- tryCatch(fread(names14d[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res14d, tempRes)
  res14d <- rbindlist(tempList)
}

# summarize
summary14a <- summarize_raw(res14a)
summary14b <- summarize_raw(res14b)
summary14c <- summarize_raw(res14c)
summary14d <- summarize_raw(res14d)

# save summaries
write.table(summary14a, paste0(outputDir, "ind3d_changeeff_maf01.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary14b, paste0(outputDir, "ind3d_changepi0_maf01.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary14c, paste0(outputDir, "bincor_changeeff_maf01.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary14d, paste0(outputDir, "rep2d_changeeff_maf01.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

#------------------------------------------------#

# MAF 10%
# read raw output files
res14a <- c()
for (file_it in 1:length(names14a_maf10)) {
  tempRes <- tryCatch(fread(names14a_maf10[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res14a, tempRes)
  res14a <- rbindlist(tempList)
}

# Read 14b
res14b <- c()
for (file_it in 1:length(names14b_maf10)) {
  tempRes <- tryCatch(fread(names14b_maf10[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res14b, tempRes)
  res14b <- rbindlist(tempList)
}

# Read 14c
res14c <- c()
for (file_it in 1:length(names14c_maf10)) {
  tempRes <- tryCatch(fread(names14c_maf10[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res14c, tempRes)
  res14c <- rbindlist(tempList)
}

# Read 14d
res14d <- c()
for (file_it in 1:length(names14d_maf10)) {
  tempRes <- tryCatch(fread(names14d_maf10[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res14d, tempRes)
  res14d <- rbindlist(tempList)
}

# summarize
summary14a <- summarize_raw(res14a)
summary14b <- summarize_raw(res14b)
summary14c <- summarize_raw(res14c)
summary14d <- summarize_raw(res14d)

# save summaries
write.table(summary14a, paste0(outputDir, "ind3d_changeeff_maf10.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary14b, paste0(outputDir, "ind3d_changepi0_maf10.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary14c, paste0(outputDir, "bincor_changeeff_maf10.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary14d, paste0(outputDir, "rep2d_changeeff_maf10.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')


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


# S Fig 14A
ind3d_changeeff_maf01 <- fread(paste(outputDir, "ind3d_changeeff_maf01.txt"), data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(!(Method %in% c("DACT", "HDMT")))

# plot 14A
ind3d_changeeff_fdp_maf01_plot <- ggplot(data=ind3d_changeeff_maf01, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  ylab("FDP (3D Pleiotropy)") +
  xlab("Min Effect Magnitude") +
 	ylim(c(0, 0.3)) +  
	xlim(c(1.6, 1.8)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))
ind3d_changeeff_fdp_maf01_plot

# S Fig 14B
ind3d_changepi0_maf01 <- fread(paste0(outputDir, "ind3d_changepi0_maf01.txt"), data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08) %>%
  filter(!(Method %in% c("DACT", "HDMT")))

# plot S Fig 14B
ind3d_changepi0_fdp_maf01_plot <- ggplot(data=ind3d_changepi0_maf01, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
	ylim(c(0, 0.3)) +   
	ylab("FDP (3D Pleiotropy)") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))
ind3d_changepi0_fdp_maf01_plot

# S Fig 14C
bincor_changeeff_maf01 <- fread(paste0(outputDir, "bincor_changeeff_maf01.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot s fig 14c
bincor_changeeff_fdp_maf01_plot <- ggplot(data=bincor_changeeff_maf01, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Correlated Test Stat.)") +
	ylim(c(0, 0.3)) +   
	xlim(c(1.6, 1.8)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changeeff_fdp_maf01_plot

# s fig 14 d
rep2d_changeeff_maf01 <- fread(paste0(outputDir, "rep2d_changeeff_maf01.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(!(Method %in% c("DACT", "HDMT")))

# plot s fig 14 d
rep2d_changeeff_fdp_maf01_plot <- ggplot(data=rep2d_changeeff_maf01 %>% mutate(FDP = ifelse(Method %in% c("HDMT", "DACT"), 1, FDP)),
                                         aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Replication)") +
  xlim(c(1.6, 1.8)) +
	ylim(c(0, 0.3)) + 
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
rep2d_changeeff_fdp_maf01_plot


#-----------------------------------------------------------#
# Now MAF 10%
# S Fig 14A
ind3d_changeeff_maf10 <- fread(paste0(outputDir, "ind3d_changeeff_maf10.txt"), data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(!(Method %in% c("DACT", "HDMT")))

# plot 14A
ind3d_changeeff_fdp_maf10_plot <- ggplot(data=ind3d_changeeff_maf10, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  ylab("FDP (3D Pleiotropy)") +
  xlab("Min Effect Magnitude") +
  xlim(c(0.4, 0.6)) +
	ylim(c(0, 0.3)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

# S Fig 14B
ind3d_changepi0_maf10 <- fread(paste0(outputDir, "ind3d_changepi0_maf10.txt"), data.table=F) %>%
  filter(Method != "DACTb" & Method != "DACTorig" & Method != "HDMTorig") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08) %>%
  filter(!(Method %in% c("DACT", "HDMT")))

# plot S Fig 14B
ind3d_changepi0_fdp_maf10_plot <- ggplot(data=ind3d_changepi0_maf10, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  ylab("FDP (3D Pleiotropy)") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

# S Fig 14C
bincor_changeeff_maf10 <- fread(paste0(outputDir, "bincor_changeeff_maf10.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))

# plot s fig 14c
bincor_changeeff_fdp_maf10_plot <- ggplot(data=bincor_changeeff_maf10, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Correlated Test Stat.)") +
  xlim(c(0.4, 0.6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))

# s fig 14 d
rep2d_changeeff_maf10 <- fread(paste0(outputDir, "rep2d_changeeff_maf10.txt"), data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(!(Method %in% c("DACT", "HDMT")))

# plot s fig 14 d
rep2d_changeeff_fdp_maf10_plot <- ggplot(data=rep2d_changeeff_maf10 %>% mutate(FDP = ifelse(Method %in% c("HDMT", "DACT"), 1, FDP)),
                                         aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Replication)") +
  xlim(c(0.4, 0.6)) +
	ylim(c(0, 0.3)) +  
	theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))


# put together figure 
ind3d_maf01_plot <- plot_grid(ind3d_changeeff_fdp_maf01_plot + theme(legend.position = "none"),
                        ind3d_changepi0_fdp_maf01_plot + theme(legend.position = "none"),
                        bincor_changeeff_fdp_maf01_plot + theme(legend.position = "none"),
                        rep2d_changeeff_fdp_maf01_plot + theme(legend.position = "none"),
												ind3d_changeeff_fdp_maf10_plot + theme(legend.position = "none"),
                        ind3d_changepi0_fdp_maf10_plot + theme(legend.position = "none"),
                        bincor_changeeff_fdp_maf10_plot + theme(legend.position = "none"),
                        rep2d_changeeff_fdp_maf10_plot + theme(legend.position = "none"),
                        labels=c("A", "B", "C", "D", "E", "F", "G", "H"), ncol=2, label_size=22)
ind3d_maf01_legend <- get_legend(bincor_changeeff_fdp_maf01_plot +  theme(legend.direction="horizontal",
                                                              legend.justification="center",
                                                              legend.box.just="bottom"))
plot_grid(ind3d_maf01_plot, ind3d_maf01_legend, ncol=1, rel_heights=c(1, 0.1))

#ggsave('SFig14_full.pdf', width=18, height=20)





















