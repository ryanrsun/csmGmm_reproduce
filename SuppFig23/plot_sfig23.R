# Collect results and plot Supp Figure 23
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig23/origOutput"
names23a <- paste0("sim_n1k_j100k_ind3d_changeeff_noAss_aID", 1:400, ".txt")
names23b <- paste0("sim_n1k_j100k_ind3d_changepi0_noAss_aID", 1:160, ".txt")
#names23c <- paste0("sim_n1k_j100k_bincor_changeeff_noass_aID", 1:400, ".txt")
names23c <- paste0("SFig23C_aID", 1:400, ".txt")
names23d <- paste0("sim_n1k_j100k_rep2d_changeeff_noAss_aID", 1:400, ".txt")
#-----------------------------------------#

# read raw output files
setwd(outputDir)
res23a <- c()
for (file_it in 1:length(names23a)) {
  tempRes <- tryCatch(fread(names23a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res23a, tempRes)
  res23a <- rbindlist(tempList)
}

# Read 23b
res23b <- c()
for (file_it in 1:length(names23b)) {
  tempRes <- tryCatch(fread(names23b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res23b, tempRes)
  res23b <- rbindlist(tempList)
}

# Read 23c
res23c <- c()
for (file_it in 1:length(names23c)) {
  tempRes <- tryCatch(fread(names23c[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res23c, tempRes)
  res23c <- rbindlist(tempList)
}

# Read 23d
res23d <- c()
for (file_it in 1:length(names23d)) {
  tempRes <- tryCatch(fread(names23d[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res23d, tempRes)
  res23d <- rbindlist(tempList)
}

# summarize
summary23a <- summarize_raw(res23a)
summary23b <- summarize_raw(res23b)
summary23c <- summarize_raw(res23c)
summary23d <- summarize_raw(res23d)

# save summaries
setwd(outputDir)
write.table(summary23a, "ind3d_changeeff_noAssumption.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary23b, "ind3d_changepi0_noAssumption.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary23c, "cor2d_changeeff_noAssumption_final.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary23d, "rep2d_changeeff_noAssumption.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

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


# s fig 23 a
setwd(outputDir)
ind3d_changeeff_noAss <- fread("ind3d_changeeff_noAssumption.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  filter(Method != "DACT") %>%
  filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df"))) 


# plot s fig 23 a
ind3d_changeeff_fdp_noAss_plot <- ggplot(data=ind3d_changeeff_noAss, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  ylab("FDP (3D Pleiotropy)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.3)) + 
  #xlim(c(1, 1.2)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))
ind3d_changeeff_fdp_noAss_plot


# s fig 23 b
setwd(outputDir)
ind3d_changepi0_noAss <- fread("ind3d_changepi0_noAssumption.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  filter(Method != "DACT") %>%
  filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  

# plot
ind3d_changepi0_fdp_noAss_plot <- ggplot(data=ind3d_changepi0_noAss, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  ylab("FDP (3D Pleiotropy)") +
  ylim(c(0, 0.3)) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))
ind3d_changepi0_fdp_noAss_plot


# s fig 23 c
setwd(outputDir)
bincor_changeeff_noAssumption <- fread("cor2d_changeeff_noAssumption_final.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  #filter(Method != "New") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "c-csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("c-csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) 
# plot sfig23 c
bincor_changeeff_fdp_noAss_plot <- ggplot(data=bincor_changeeff_noAssumption, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Correlated Test Stat.)") +
  ylim(c(0, 0.45)) +
  #xlim(c(1.8, 2)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changeeff_fdp_noAss_plot

# s fig 23 d
setwd(outputDir)
rep2d_changeeff_noAss <- fread("rep2d_changeeff_noAssumption.txt", data.table=F) %>%
  #filter(Method != "DACTb" & Method != "DACT" & Method != "HDMT") %>%
  filter(Method != "DACTb") %>%
  filter(Method != "DACT") %>%
  filter(Method != "HDMT") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "r-csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("r-csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df"))) 

# plot rep2d changeeff fdp
rep2d_changeeff_fdp_noAss_plot <- ggplot(data=rep2d_changeeff_noAss %>% mutate(FDP = ifelse(Method %in% c("HDMT", "DACT"), 1, FDP)),
                                         aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Replication)") +
  ylim(c(0, 0.45)) + 
  #xlim(c(1.8, 2)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
rep2d_changeeff_fdp_noAss_plot

# put together figure
ind3d_noAss_plot <- plot_grid(ind3d_changeeff_fdp_noAss_plot + theme(legend.position = "none"),
                              ind3d_changepi0_fdp_noAss_plot + theme(legend.position = "none"),
                              bincor_changeeff_fdp_noAss_plot + theme(legend.position = "none"),
                              rep2d_changeeff_fdp_noAss_plot + theme(legend.position = "none"),
                              labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
ind3d_noAss_legend <- get_legend(ind3d_changeeff_fdp_noAss_plot +  theme(legend.direction="horizontal",
                                                                         legend.justification="center",
                                                                         legend.box.just="bottom"))
plot_grid(ind3d_noAss_plot, ind3d_noAss_legend, ncol=1, rel_heights=c(1, 0.1))
setwd(outputDir)
ggsave('ind3d_correp_noAssumption.pdf', width=18, height=12)























