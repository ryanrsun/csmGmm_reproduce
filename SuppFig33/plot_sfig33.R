# Collect results and plot Supp Figure 33
library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
setwd('/rsrch3/home/biostatistics/rsun3/github/ancillaryFunctions')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#-----------------------------------------#
# change to where the output files are stored
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig33/output"
names33a <- paste0("SFig33A_aID", 1:400, ".txt")
names33b <- paste0("SFig33B_aID", 1:160, ".txt")
#-----------------------------------------#

# read raw output files
setwd(outputDir)
res33a <- c()
for (file_it in 1:length(names33a)) {
  tempRes <- tryCatch(fread(names33a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  if (ncol(tempRes) == 27) {tempRes = tempRes %>% mutate(nRej50df = NA)} 
  tempList <- list(res33a, tempRes)
  res33a <- rbindlist(tempList, use.names=TRUE)
}

res33b <- c()
for (file_it in 1:length(names33b)) {
  tempRes <- tryCatch(fread(names33b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  if (ncol(tempRes) == 27) {tempRes = tempRes %>% mutate(nRej50df = NA)}
  tempList <- list(res33b, tempRes)
  res33b <- rbindlist(tempList, use.names=TRUE)
}

# summarize
summary33a <- summarize_raw(res33a)
summary33b <- summarize_raw(res33b)

# save summaries
setwd(outputDir)
write.table(summary33a, "med2d_changeeff_fullLik.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary33b, "med2d_changepi0_fullLik.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')


#------------------------------------------------#
# plotting starts

# define colors
gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[2] <- "brown"
mycols[4] <- "black"
mycols[5] <- "blue"


# full likelihood change eff
setwd(outputDir)
ind2d_changeeff_fullLik <- fread("med2d_changeeff_fullLik.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "DACT", "Full", Method)) %>% 
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Full", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  

# plot changeeff FDP
ind2d_changeeff_fullLik_fdp_plot <- ggplot(data=ind2d_changeeff_fullLik, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (2D Pleiotropy)") +
  xlab("Mean Magnitude") +
  ylim(c(0, 0.25)) + 
  xlim(c(1.25, 2.5)) + 
  #xlim(c(0.4, 0.6)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_fullLik_fdp_plot


# plot changeeff power
ind2d_changeeff_fullLik_power_plot <- ggplot(data=ind2d_changeeff_fullLik, aes(x=minEff1, y=Power, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("Power (2D Pleiotropy)") +
  xlab("Mean Magnitude") +
  ylim(c(0, 1)) + 
  xlim(c(1.25, 2.5)) + 
  #xlim(c(0.4, 0.6)) +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_fullLik_power_plot


#-------------------------------------------------------------------------------------------#
# full likelihood change pi0
setwd(outputDir)
ind2d_changepi0_fullLik <- fread("med2d_changepi0_fullLik.txt", data.table=F) %>%
  mutate(Method = ifelse(Method == "DACT", "Full", Method)) %>% 
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Full", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  

# plot changepi0 FDP
ind2d_changepi0_fullLik_fdp_plot <- ggplot(data=ind2d_changepi0_fullLik, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (2D Pleiotropy)") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylim(c(0, 0.2)) + 
  #xlim(c(0.4, 0.6)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changepi0_fullLik_fdp_plot


# plot changepi0 power
ind2d_changepi0_fullLik_power_plot <- ggplot(data=ind2d_changepi0_fullLik, aes(x=minEff1, y=Power, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("Power (2D Pleiotropy)") +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0) = Proportion of (", alpha[j], "=0", ", ", beta[j] != 0, ")"))) +
  ylim(c(0, 1)) + 
  #xlim(c(0.4, 0.6)) +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changepi0_fullLik_power_plot


# save
ind2d_fulllik_plot <- plot_grid(ind2d_changeeff_fullLik_fdp_plot + theme(legend.position = "none"),
                                ind2d_changeeff_fullLik_power_plot + theme(legend.position = "none"),
                            ind2d_changepi0_fullLik_fdp_plot + theme(legend.position = "none"),
                            ind2d_changepi0_fullLik_power_plot + theme(legend.position = "none"),
                            labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
ind2d_fulllik_legend <- get_legend(ind2d_changeeff_fullLik_fdp_plot +  theme(legend.direction="horizontal",
                                                                   legend.justification="center",
                                                                   legend.box.just="bottom"))
plot_grid(ind2d_fulllik_plot, ind2d_fulllik_legend, ncol=1, rel_heights=c(1, 0.15))
setwd(outputDir)
ggsave('ind2d_fulllik_fixedeff.pdf', width=20, height=12)











