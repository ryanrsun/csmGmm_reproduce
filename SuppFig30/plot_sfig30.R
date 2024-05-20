# Collect results and plot Supp Figure 30
library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
devtools::install_github("ryanrsun/csmGmm")
setwd('../supportingCode')
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#-----------------------------------------#
# change to where the output files are stored
outputDir <- "/rsrch3/home/biostatistics/rsun3/empBayes/reproduce/SuppFig30/output"
names30f1 <- paste0("Fig30A_aID", 1:400, "_high.txt")
names30f2 <- paste0("Fig30A_aID", 1:400, "_med.txt")
names30f3 <- paste0("Fig30A_aID", 1:400, "_low.txt")

#-----------------------------------------#

# read raw output files
setwd(outputDir)
res30f1 <- c()
for (file_it in 1:length(names30f1)) {
  tempRes <- tryCatch(fread(names30f1[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res30f1, tempRes)
  res30f1 <- rbindlist(tempList)
}

# Read 30f2
res30f2 <- c()
for (file_it in 1:length(names30f2)) {
  tempRes <- tryCatch(fread(names30f2[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res30f2, tempRes)
  res30f2 <- rbindlist(tempList)
}

# Read 30f3
res30f3 <- c()
for (file_it in 1:length(names30f3)) {
  tempRes <- tryCatch(fread(names30f3[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res30f3, tempRes)
  res30f3 <- rbindlist(tempList)
}

# summarize
summary30f1 <- summarize_raw(res30f1)
summary30f2 <- summarize_raw(res30f2)
summary30f3 <- summarize_raw(res30f3)

# save summaries
setwd(outputDir)
write.table(summary30f1, "bincor_changeeff_onemix_high.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary30f2, "bincor_changeeff_onemix_med.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary30f3, "bincor_changeeff_onemix_low.txt", append=F,quote=F, row.names=F, col.names=T, sep='\t')


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

# read data
setwd(outputDir)
bincor_changeeff_use1_high <- fread("bincor_changeeff_onemix_high.txt", data.table=F) %>%
  filter(Method == "New") %>%
  mutate(Method = "c-csmGmm-1-high") %>%
  select(minEff1, Method, Power, FDP)
bincor_changeeff_use1_low <- fread("bincor_changeeff_onemix_low.txt", data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "c-csmGmm-1-low") %>%
  select(minEff1, Method, Power, FDP)
bincor_changeeff_use1_orig <- fread("bincor_changeeff_onemix_med.txt", data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "c-csmGmm-1-med") %>%
  select(minEff1, Method, Power, FDP)
bincor_changeeff_standard <- fread("Fig3c_summary_final.txt", data.table=F) %>%
  filter(Method == "New")  %>%
  mutate(Method = "c-csmGmm-standard") %>%
  select(minEff1, Method, Power, FDP)
bincor_changeeff_OM <- rbind(bincor_changeeff_use1_high,
                            bincor_changeeff_use1_low,
                            bincor_changeeff_use1_orig,
                            bincor_changeeff_standard)

# plot bincor
bincor_changeeff_OM_plot <- ggplot(data=bincor_changeeff_OM, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  ylab("FDP (2D Correlated Test Statistics)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.2)) + 
  xlim(c(0.2, 0.4)) +
  geom_hline(yintercept = 0.1, linetype=2, color="darkgrey", lwd=1.4) +
  scale_linetype_manual(values=c(4,3,2,1)) +
  scale_color_manual(values=c("purple", "orange", "darkgreen", mycols[1])) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=20)) +
  theme(legend.key.size = unit(3,"line"))
bincor_changeeff_OM_plot


# save
setwd(outputDir)
ggsave("bincor_changeeff_onemix.pdf", width=12, height=6)











