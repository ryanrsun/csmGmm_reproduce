# Plot supp fig 9
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggformula)
library(devtools)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

outputDir <- here::here("SuppFig9", "output/")
names9a <- paste0(outputDir, "sim_n1k_j100k_ind5d_changeeff2_combined_aID", 1:2000, ".txt")
names9b <- paste0(outputDir, "sim_n1k_j100k_ind5d_changepi0_combined_aID", 1:800, ".txt")
names9c <- paste0(outputDir, "sim_n1k_j100k_ind6d_changeeff2_combined_aID", 1:2000, ".txt")
names9d <- paste0(outputDir, "sim_n1k_j100k_ind6d_changepi0_combined_aID", 1:800, ".txt")

# read data
res9a <- c()
for (file_it in 1:length(names9a)) {
  tempRes <- tryCatch(fread(names9a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res9a, tempRes)
  res9a <- rbindlist(tempList)
}

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
summary9a <- summarize_raw(res9a)
summary9b <- summarize_raw(res9b)
summary9c <- summarize_raw(res9c)
summary9d <- summarize_raw(res9d)

# save
write.table(summary9b, paste0(outputDir, "med2d_changepi0_5d.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary9b, paste0(outputDir, "med2d_changepi0_5d.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary9c, paste0(outputDir, "med2d_changeeff_ind6d_v2.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary9d, paste0(outputDir, "ind6d_changepi0.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

#-------------------------------------------------------#
# start plotting

# s fig 9a
ind5d_changeeff <- fread(paste0(outputDir, "med2d_changeeff_correctpi_5d.txt"), data.table=F) %>%
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
ind5d_changepi0 <- fread(paste0(outputDir, "med2d_changepi0_5d.txt"), data.table=F) %>%
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
ind6d_changeeff <- fread(paste0(outputDir, "med2d_changeeff_ind6d_v2.txt"), data.table=F) %>%
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
ind6d_changepi0 <- fread(paste0(outputDir, "ind6d_changepi0.txt"), data.table=F) %>%
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
#ggsave('ind56d_fdp.pdf', width=20, height=12)









