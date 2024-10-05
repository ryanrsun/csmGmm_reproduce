# Collect results and plot Supp Figure 11
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

# change to where the output files are stored
outputDir <- here::here("SuppFig11/output/")
newFiles <- paste0(outputDir, "sim_n1k_j100k_ind7d_hard_changeeff_new_aID", 1:2000, ".txt")
df7Files <- paste0(outputDir, "sim_n1k_j100k_ind7d_hard_changeeff_7df_aID", 1:2000, ".txt")
df50Files <- paste0(outputDir, "sim_n1k_j100k_ind7d_hard_changeeff_50df_aID", 1:2000, ".txt")
kernelFiles <- paste0(outputDir, "sim_n1k_j100k_ind7d_hard_changeeff_kernel_aID", 1:2000, ".txt")

fullResNew <- c()
fullRes7df <- c()
fullRes50df <- c()
fullResKernel <- c()

# read raw output files
for (file_it in 2:max(length(newFiles), length(df7Files), length(df50Files), length(kernelFiles))) {
  if (length(newFiles) >= file_it) {
    tempResNew <- tryCatch(fread(newFiles[file_it]), error=function(e) e)
    if (class(tempResNew)[1] %in% c("simpleError", "simpleWarning")) {next} 
    tempListNew <- list(fullResNew, tempResNew)
    fullResNew <- rbindlist(tempListNew)
  }

  if (length(df7Files) >= file_it) {
    tempRes7df <- fread(df7Files[file_it])
    tempList7df <- list(fullRes7df, tempRes7df)
    fullRes7df <- rbindlist(tempList7df)
  }

  if (length(df50Files) >= file_it) {
    tempRes50df <- fread(df50Files[file_it])
    tempList50df <- list(fullRes50df, tempRes50df)
    fullRes50df <- rbindlist(tempList50df)
  }

  if (length(kernelFiles) >= file_it) {
    tempResKernel <- fread(kernelFiles[file_it])
    tempListKernel <- list(fullResKernel, tempResKernel)
    fullResKernel <- rbindlist(tempListKernel)
  }
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

# plot
ind7d_changeeff <- outDF  %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "DACT", "Kernel", "locfdr7df",
                                          "locfdr50df", "HDMT"))) 

ind7d_changeeff_fdp_plot <- ggplot(data=ind7d_changeeff, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  xlab("Min Effect Magnitude") +
  ylab("FDP (K=7)") +
  ylim(c(0, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_color_manual(values=mycols[c(1, 3:6)]) +
  scale_linetype_manual(values=c(1, 3:6)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
ind7d_changeeff_fdp_plot

ggsave("ind7d_changeeff_fdp.pdf")

