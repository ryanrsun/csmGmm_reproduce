# Collect results and plot Supp Figure 18
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
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig18/output"
names18a <- paste0("SFig18A_aID", 201:700, ".txt")
names18b <- paste0("SFig18B_aID", 201:700, ".txt")
#-----------------------------------------#

# read raw output files
setwd(outputDir)
res18a <- c()
for (file_it in 1:length(names18a)) {
  tempRes <- tryCatch(fread(names18a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res18a, tempRes)
  res18a <- rbindlist(tempList)
}

# Read 18b
res18b <- c()
for (file_it in 1:length(names18b)) {
  tempRes <- tryCatch(fread(names18b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res18b, tempRes)
  res18b <- rbindlist(tempList)
}

# summarize
summary18a <- summarize_raw(res18a)
summary18b <- summarize_raw(res18b)

# save summaries
setwd(outputDir)
write.table(summary18a, "bincor_changeeff_cor5_large.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary18b, "bincor_changeeff_cor8_large.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')

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


#---------------------------------------------------------------------------#
# read sfig 18 A
setwd(outputDir)
bincor_changeeff5 <- fread("bincor_changeeff_cor5_large.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "c-csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("c-csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 >= 0.2 & minEff1 <= 0.39)

# plot sfig 18a
bincor_changeeff5_fdp_plot <- ggplot(data=bincor_changeeff5, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Correlation = 0.5)") +
  #ylim(c(0, 0.45)) + 
  xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changeeff5_fdp_plot

# plot sfig18c
bincor_changeeff5_power_plot <- ggplot(data=bincor_changeeff5 %>% filter(Method == "c-csmGmm"),
                                       aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  #geom_hline(yintercept = 0.1, linetype=2, color="grey") +
  scale_color_manual(values=mycols[1]) +
  scale_linetype_manual(values=1) +
  xlab("Min Effect Magnitude") + ylab("Power (2D Correlation = 0.5)") +
  xlim(c(0.2, 0.4)) +
  ylim(c(0, 0.8)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changeeff5_power_plot


# read sfig 18 b
setwd(outputDir)
bincor_changeeff8 <- fread("bincor_changeeff_cor8_large.txt", data.table=F) %>%
  filter(Method != "DACTb") %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "c-csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("c-csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) %>%
  filter(minEff1 >= 0.2 & minEff1 <= 0.39)

# plot sfig 18 b
bincor_changeeff8_fdp_plot <- ggplot(data=bincor_changeeff8,
                                     aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[1:6]) +
  scale_linetype_manual(values=1:6) +
  xlab("Min Effect Magnitude") + ylab("FDP (2D Correlation = 0.8)") +
  #ylim(c(0, 0.45)) + 
  xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changeeff8_fdp_plot

# plot sfig 18d
bincor_changeeff8_power_plot <- ggplot(data=bincor_changeeff8 %>% filter(Method %in% c("c-csmGmm")),
                                     aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method), lwd=1.2) +
  scale_color_manual(values=mycols[1:2]) +
  scale_linetype_manual(values=1:2) +
  xlab("Min Effect Magnitude") + ylab("Power (2D Correlation = 0.8)") +
  ylim(c(0, 0.8)) + 
  xlim(c(0.2, 0.4)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20))+
  theme(legend.key.size = unit(3,"line"))
bincor_changeeff8_power_plot


# save
morecore_changeeeff_fdp_plot <- plot_grid(bincor_changeeff5_fdp_plot + theme(legend.position = "none"),
                            bincor_changeeff8_fdp_plot + theme(legend.position = "none"),
                            bincor_changeeff5_power_plot + theme(legend.position = "none"),
                            bincor_changeeff8_power_plot + theme(legend.position = "none"),
                            labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
morecore_changeeeff_fdp_legend <- get_legend(bincor_changeeff5_fdp_plot +  theme(legend.direction="horizontal",
                                                                   legend.justification="center",
                                                                   legend.box.just="bottom"))
plot_grid(morecore_changeeeff_fdp_plot, morecore_changeeeff_fdp_legend, ncol=1, rel_heights=c(1, 0.15))
setwd(outputDir)
ggsave('morecor_changeeff.pdf', width=20, height=12)












