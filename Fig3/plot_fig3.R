# Collect results and plot Figure 3
library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(csmGmm)

#-----------------------------------------#
# change to where the output files are stored

# source the .R scripts from the supportingCode/ folder in the csmGmm_reproduce repository
setwd('/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SupportingCode/')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/test/output"
names3a <- paste0("Fig3A_aID", 1:620, ".txt")
names3b <- paste0("Fig3B_aID", 1:160, ".txt")
names3c <- paste0("Fig3C_aID", 401:1240, ".txt")
names3d <- paste0("Fig3D_aID", 1:1000, ".txt")
#-----------------------------------------#

# read raw output files
setwd(outputDir)
res3a <- c()
for (file_it in 1:length(names3a)) {
  tempRes <- tryCatch(fread(names3a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res3a, tempRes)
  res3a <- rbindlist(tempList)
}

# Read 3b
res3b <- c()
for (file_it in 1:length(names3b)) {
  tempRes <- tryCatch(fread(names3b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res3b, tempRes)
  res3b <- rbindlist(tempList)
}

# Read 3c
res3c <- c()
for (file_it in 1:length(names3c)) {
  tempRes <- tryCatch(fread(names3c[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res3c, tempRes)
  res3c <- rbindlist(tempList)
}

# Read 3d
res3d <- c()
for (file_it in 1:length(names3d)) {
  tempRes <- tryCatch(fread(names3d[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res3d, tempRes)
  res3d <- rbindlist(tempList)
}

# summarize
summary3a <- summarize_raw(res3a)
summary3b <- summarize_raw(res3b)
#summary3c <- summarize_raw(res3c, cor=TRUE)
summary3c <- summarize_raw(res3c)
summary3d <- summarize_raw(res3d)

# save summaries
setwd(outputDir)
write.table(summary3a, "Fig3a_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary3b, "Fig3b_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary3c, "Fig3c_summary_final.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary3d, "Fig3d_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

#-----------------------------------------#
# start plotting

# define colors
gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"

# Figure 3A
setwd(outputDir)
Fig3A_data <- fread("Fig3a_summary.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 >= 0.25 & minEff1 <= 0.44) %>%
  filter(!is.na(Method))

Fig3A_plot <- ggplot(data=Fig3A_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (3D Pleiotropy)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.3)) + xlim(c(0.25, 0.45)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

# Figure 3B
setwd(outputDir)
Fig3B_data <- fread("Fig3b_summary.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 <= 0.08) %>%
  filter(!is.na(Method))

Fig3B_plot <- ggplot(data=Fig3B_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (3D Pleiotropy)") +
  ylim(c(0, 0.3)) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

# Figure 3C
setwd(outputDir)
Fig3C_data <- fread("Fig3c_summary_final.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "c-csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("c-csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 >= 0.2 & minEff1 <= 0.39) %>%
  filter(!is.na(Method))

Fig3C_plot <- ggplot(data=Fig3C_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Correlated Test Stat.)") +
  ylim(c(0, 0.45)) + xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))

# Figure 3D
setwd(outputDir)
Fig3D_data <- fread("Fig3d_summary.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "r-csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("r-csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df"))) %>%
  filter(minEff1 >= 0.30 & minEff1 <= 0.49) %>%
  filter(!is.na(Method))

Fig3D_plot <- ggplot(data=Fig3D_data %>% mutate(FDP = ifelse(Method %in% c("HDMT", "DACT"), 1, FDP)), aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Replication)") +
  ylim(c(0, 0.45)) + xlim(c(0.3, 0.5)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))


# put together figure 3
Fig3_plot <- plot_grid(Fig3A_plot + theme(legend.position = "none"),
                        Fig3B_plot + theme(legend.position = "none"),
                        Fig3C_plot + theme(legend.position = "none"),
                        Fig3D_plot + theme(legend.position = "none"),
                                           labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
Fig3_legend <- get_legend(Fig3C_plot +  theme(legend.direction="horizontal",
                                                                                        legend.justification="center",
                                                                                        legend.box.just="bottom"))
plot_grid(Fig3_plot, Fig3_legend, ncol=1, rel_heights=c(1, 0.1))
setwd(outputDir)
ggsave('Fig3.pdf', width=18, height=12)











