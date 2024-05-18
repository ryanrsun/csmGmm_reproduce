# Collect results and plot Figure 1
library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
devtools::install.packages("ryanrsun/csmGmm")
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#-----------------------------------------#
# change to where the output files are stored
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/Fig1/origOutput"
names1a <- paste0("Fig1A_aID", 1:800, ".txt")
names1b <- paste0("Fig1B_aID", 1:320, ".txt")
names1c <- paste0("Fig1C_aID", 1:800, ".txt")
names1d <- paste0("Fig1D_aID", 1:320, ".txt")
#-----------------------------------------#

# read raw output files
setwd(outputDir)
res1a <- c()
for (file_it in 1:length(names1a)) {
  tempRes <- tryCatch(fread(names1a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res1a, tempRes)
  res1a <- rbindlist(tempList)
}

# Read 1b
res1b <- c()
for (file_it in 1:length(names1b)) {
  tempRes <- tryCatch(fread(names1b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res1b, tempRes)
  res1b <- rbindlist(tempList)
}

# Read 1c
res1c <- c()
for (file_it in 1:length(names1c)) {
  tempRes <- tryCatch(fread(names1c[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res1c, tempRes)
  res1c <- rbindlist(tempList)
}

# Read 1d
res1d <- c()
for (file_it in 1:length(names1d)) {
  tempRes <- tryCatch(fread(names1d[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res1d, tempRes)
  res1d <- rbindlist(tempList)
}

# summarize
summary1a <- summarize_raw(res1a)
summary1b <- summarize_raw(res1b)
summary1c <- summarize_raw(res1c)
summary1d <- summarize_raw(res1d)

# save summaries
setwd(outputDir)
write.table(summary1a, "Fig1a_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary1b, "Fig1b_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary1c, "Fig1c_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary1d, "Fig1d_summary.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

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

# Fig 1A
setwd(outputDir)
Fig1a_data <- fread("Fig1a_summary.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.19)

Fig1a_plot <- ggplot(data=Fig1a_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (Mediation, No Alt)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.2)) + xlim(c(0, 0.2)) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# Fig 1B
setwd(outputDir)
Fig1b_data <- fread("Fig1b_summary.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

Fig1b_plot <-ggplot(data=Fig1b_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP (Mediation, No Alt)") +
  ylim(c(0, 0.2)) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))


# Figure 1C
setwd(outputDir)
Fig1c_data <- fread("Fig1c_summary.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 <= 0.19)

Fig1c_plot <- ggplot(data=Fig1c_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (Mediation, With Alt)") +
  ylim(c(0, 0.2)) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))

# Figure 1D
setwd(outputDir)
Fig1d_data <- fread("Fig1d_summary.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "DACT", "DACT", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  %>%
  filter(minEff1 <= 0.08)

Fig1d_plot <- ggplot(data=Fig1d_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylab("FDP (Mediation, With Alt)") +
  ylim(c(0, 0.45)) +
  geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))


