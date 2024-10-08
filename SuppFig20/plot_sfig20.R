# Collect results and plot Supp Figure 20

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig20/plot_sfig20.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig20/plot_sfig20.R")

library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(devtools)
library(csmGmm)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("SuppFig20", "output/")
names20a <- paste0(outputDir, "SFig20A_aID", 201:700, ".txt")
names20b <- paste0(outputDir, "SFig20B_aID", 201:900, ".txt")
#-----------------------------------------#

# read raw output files
res20a <- c()
for (file_it in 1:length(names20a)) {
  tempRes <- tryCatch(fread(names20a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next} 
  if (ncol(tempRes) > 28) {
    tempRes$powerHDMT = tempRes$powerHMDT
    tempRes <- tempRes %>% select(-powerHMDT)
  }
  tempList <- list(res20a, tempRes)
  res20a <- rbindlist(tempList)
}

# Read 20b
res20b <- c()
for (file_it in 1:length(names20b)) {
  tempRes <- tryCatch(fread(names20b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  if (ncol(tempRes) > 28) {
    tempRes$powerHDMT = tempRes$powerHMDT
    tempRes <- tempRes %>% select(-powerHMDT)
  }  
  tempList <- list(res20b, tempRes)
  res20b <- rbindlist(tempList)
}

# summarize
summary20a <- summarize_raw(res20a)
summary20b <- summarize_raw(res20b)

# save summaries
write.table(summary20a, paste0(outputDir, "bincor_changeeff_emrho.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary20b, paste0(outputDir, "bincor_changeeff_emrho_ez.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

#------------------------------------------------#
# plotting starts

# define colors
gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(7)
mycols[4] <- "black"
mycols[5] <- "blue"

# read sfig 20a
emrho_changeeff <- fread(paste0(outputDir, "bincor_changeeff_emrho.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "New", "c-csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "HDMT", "c-csmGmm-rho", Method)) %>%
  filter(Method %in% c("c-csmGmm", "c-csmGmm-rho")) %>%
  mutate(FDP = ifelse(is.na(FDP), 0, FDP))

# plot sfig20a
emrho_changeeff_fdp_plot <- ggplot(data=emrho_changeeff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[c(1, 7)]) +
  scale_linetype_manual(values=1:2) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
emrho_changeeff_fdp_plot

# plot sfig20c
emrho_changeeff_power_plot <- ggplot(data=emrho_changeeff,
                                          aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("Power") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=mycols[c(1,7)]) +
  scale_linetype_manual(values=c(1,2)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
emrho_changeeff_power_plot


# read sfig20b
emrho_ez <- fread(paste0(outputDir, "bincor_changeeff_emrho_ez.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "New", "c-csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "HDMT", "c-csmGmm-rho", Method))  %>%
  filter(Method %in% c("c-csmGmm", "c-csmGmm-rho")) %>%
  mutate(FDP = ifelse(is.na(FDP), 0, FDP))

# plot sfig20b
emrho_ez_fdp_plot <- ggplot(data=emrho_ez, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[c(1, 7)]) +
  scale_linetype_manual(values=1:2) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
emrho_ez_fdp_plot


# plot sfig20d
emrho_ez_power_plot <- ggplot(data=emrho_ez %>% filter(!(Method %in% c("DACT"))),
                                          aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("Power") +
  #ylim(c(0, 0.4)) +
  scale_color_manual(values=mycols[c(1, 7)]) +
  scale_linetype_manual(values=1:2) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
emrho_ez_power_plot

# em rho
emrho_plot <- plot_grid(emrho_changeeff_fdp_plot + theme(legend.position = "none"),
                       emrho_ez_fdp_plot + theme(legend.position = "none"),
                       emrho_changeeff_power_plot + theme(legend.position = "none"),
                       emrho_ez_power_plot + theme(legend.position = "none"),
                       labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
emrho_legend <- get_legend(emrho_changeeff_fdp_plot +  theme(legend.direction="horizontal",
                                                                 legend.justification="center",
                                                                 legend.box.just="bottom"))
plot_grid(emrho_plot, emrho_legend, ncol=1, rel_heights=c(1, 0.15))
#ggsave('emrho.pdf', width=20, height=12)










