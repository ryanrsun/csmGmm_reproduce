# Collect results and plot Supp Figure 32
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig32/output"
names32a <- paste0("SFig32A_aID", 1:400, ".txt")
names32b <- paste0("SFig32B_aID", 1:160, ".txt")
#-----------------------------------------#

# read raw output files
setwd(outputDir)
res32a <- c()
for (file_it in 1:length(names32a)) {
  tempRes <- tryCatch(fread(names32a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res32a, tempRes)
  res32a <- rbindlist(tempList)
}

res32b <- c()
for (file_it in 1:length(names32b)) {
  tempRes <- tryCatch(fread(names32b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res32b, tempRes)
  res32b <- rbindlist(tempList)
}

# summarize
summary32a <- summarize_raw(res32a)
summary32b <- summarize_raw(res32b)

# save summaries
setwd(outputDir)
write.table(summary32a, "med2d_changeeff_fullLik_randomeff.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary32b, "med2d_changepi0_fullLik_randomeff.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')


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


# full likelihood random effects change eff
setwd(outputDir)
ind2d_changeeff_fullLik_randomeff <- fread("med2d_changeeff_fullLik_randomeff.txt", data.table=F) %>%
  filter(Method != "Known") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
	mutate(Method = ifelse(Method == "DACT", "Full", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Full", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  

# plot changeeff FDP
ind2d_changeeff_fullLik_randomeff_fdp_plot <- ggplot(data=ind2d_changeeff_fullLik_randomeff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (2D Pleiotropy)") +
  xlab("Mean Magnitude") +
  #ylim(c(0, 0.2)) + 
  #xlim(c(0.4, 0.6)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_fullLik_randomeff_fdp_plot


# plot changeeff power
ind2d_changeeff_fullLik_randomeff_power_plot <- ggplot(data=ind2d_changeeff_fullLik_randomeff, aes(x=minEff1, y=Power, group=Method)) +
  geom_smooth(aes(linetype = Method, color=Method),lwd=1.2, se=FALSE)+
  #geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("Power (2D Pleiotropy)") +
  xlab("Mean Magnitude") +
  ylim(c(0, 1)) + 
  #xlim(c(0.4, 0.6)) +
  scale_linetype_manual(values=1:6) +
  scale_color_manual(values=mycols) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind2d_changeeff_fullLik_randomeff_power_plot


# full likelihood random eff change pi0
setwd(outputDir)
ind2d_changepi0_fullLik_randomeff <- fread("med2d_changepi0_fullLik_randomeff.txt", data.table=F) %>%
  filter(Method != "Known") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
	mutate(Method = ifelse(Method == "DACT", "Full", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Full", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT")))  

# plot changepi0 FDP
ind2d_changepi0_fullLik_randomeff_fdp_plot <- ggplot(data=ind2d_changepi0_fullLik_randomeff, aes(x=minEff1, y=FDP, group=Method)) +
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
ind2d_changepi0_fullLik_randomeff_fdp_plot


# plot changepi0 power
ind2d_changepi0_fullLik_randomeff_power_plot <- ggplot(data=ind2d_changepi0_fullLik_randomeff, aes(x=minEff1, y=Power, group=Method)) +
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
ind2d_changepi0_fullLik_randomeff_power_plot


# save
fulllik_randomeff_plot <- plot_grid(ind2d_changeeff_fullLik_randomeff_fdp_plot + theme(legend.position = "none"),
                                    ind2d_changeeff_fullLik_randomeff_power_plot + theme(legend.position = "none"),
                                ind2d_changepi0_fullLik_randomeff_fdp_plot + theme(legend.position = "none"),
                                ind2d_changepi0_fullLik_randomeff_power_plot + theme(legend.position = "none"),
                                labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
fulllik_randomeff_legend <- get_legend(ind2d_changeeff_fullLik_randomeff_fdp_plot +  theme(legend.direction="horizontal",
                                                                             legend.justification="center",
                                                                             legend.box.just="bottom"))
plot_grid(fulllik_randomeff_plot, fulllik_randomeff_legend, ncol=1, rel_heights=c(1, 0.15))
setwd(outputDir)
ggsave('ind2d_fulllik_randomeff.pdf', width=20, height=12)












