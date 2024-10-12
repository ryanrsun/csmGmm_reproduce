# Supp Fig 34

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig34/plot_sfig34.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("SuppFig34/plot_sfig34.R")

library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(data.table)
library(xtable)

# 1 is CAD and BMI
# 2 is ILCCO overall and Cardiogram CAD
# 3 is ILCCO overall and UKB BMI
# 4 is replication ILCCO overall and UKB lc
# 5 is replication CAD
# 6 is three way ILCCO overall, Cardiogram CAD, UKB BMI

outputDir <- here::here("Fig4", "output/")

# for colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# plot manhattan function
plotManhattan <- function(plotRes, chrCounts, colValues, shapeValues, ylimits, legName) {
  # arrange data by chromosome
  plotRes <- plotRes %>% arrange(Chr)
  uniqueChrs <- sort(unique(plotRes$Chr))

  # add true positions
  truePos <- rep(NA, nrow(plotRes))
  counter <- 1
  for (chr_it in 1:length(uniqueChrs)) {
    tempChr <- uniqueChrs[chr_it]
    tempDat <- plotRes %>% filter(Chr == tempChr)
    truePos[counter:(counter + nrow(tempDat) - 1)] <- rep(sum(chrCounts[1:tempChr]), nrow(tempDat)) + tempDat$BP
    counter <- counter + nrow(tempDat)
  }

  # plot
  #xBreaks <- cumsum(chrCounts[-1])
  xBreaks <- c(249182887,  492035665,  689855565,  880699603, 1061374625, 1232263370, 1391342188, 1537613317,
               1678651335, 1814062324,1948990547, 2082797100, 2197858639, 2305109387, 2407535086, 2497679391,
               2578726035, 2656695761, 2715744515, 2778629720,2826695719, 2877914096)
  xBreaksLabs <- 1:22
  xBreaksLabs[c(9, 11, 13, 15, 16, 17, 19, 20, 21)] <- ""

  plotDat <- plotRes %>% mutate(truePos = truePos)

  returnPlot <- ggplot(plotDat, aes(x=truePos, y=-log10(newLfdr), color=as.factor(cat), shape=as.factor(cat))) +
    geom_point() +
    xlab("Chromosome") + ylab("-log10(lfdr)") +
    
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(name="Chr", breaks=xBreaks, labels=xBreaksLabs, limits=c(0,2877914096 )) +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=18),
          legend.title=element_text(size=18), legend.text=element_text(size=18)) +
    guides(colour = guide_legend(override.aes = list(size=4)))


  return(returnPlot)
}


# add position information to data
s2 <- fread(paste0(outputDir, "reject_bmi_with_overall_neg5_reject_aID2.txt"))
s4 <- fread(paste0(outputDir, "reject_bmi_with_overall_neg5_reject_aID4.txt"))
s5 <- fread(paste0(outputDir, "reject_bmi_with_overall_neg5_reject_aID5.txt"))
s2new <- s2 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))
s4new <- s4 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))
s5new <- s5 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))

# for plotting axes
allZ <- fread(paste0(outputDir, "bmi_with_overall.txt")) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP))
chrCounts <- rep(0, 23)
for (chr_it in 1:22) {
  tempDat <- allZ %>% filter(Chr == chr_it)
  maxPos <- max(tempDat$BP)
  chrCounts[chr_it + 1] <- maxPos
}

# data for manhattan plot - replication
manDataRep <- rbind(s5new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(pheno = "CAD"),
                    s4new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(pheno = "LC"),
                    s2new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(pheno = "Both")) %>%
  as.data.frame(.) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP)) %>%
  arrange(newLfdr) %>%
  # distinct() keeps the first one, so we arrange first
  distinct(., chrpos, .keep_all = TRUE) %>%
  # overlap with LC, CAD pleiotropy
  mutate(Pleio = ifelse(chrpos %in% s2new$chrpos, 1, 0)) %>%
  mutate(RepCAD = ifelse(chrpos %in% s5new$chrpos, 1, 0)) %>%
  mutate(RepLC = ifelse(chrpos %in% s4new$chrpos, 1, 0)) %>%
  mutate(cat = ifelse(Pleio == 1 & RepCAD == 1 & RepLC == 1, "Pleio + Rep", "Rep Only")) %>%
  mutate(cat = ifelse(Pleio == 1 & (RepCAD == 0 | RepLC == 0), "Pleio Only", cat)) 

# separate plot
manPlotRep1 <- plotManhattan(plotRes = manDataRep %>% filter(cat == "Pleio + Rep"), chrCounts,
                            colValues=c(gg_color_hue(3)[3], "black", "darkorange"), shapeValues=c(17, 18, 15),
                            ylimits=c(0, 12.5), legName="Lung Cancer")
manPlotRep1

manPlotRep2 <- plotManhattan(plotRes = manDataRep %>% filter(cat == "Pleio Only"), chrCounts,
                             colValues=c("black", "black", "darkorange"), shapeValues=c(18, 18, 15),
                             ylimits=c(0, 12.5), legName="Lung Cancer")
manPlotRep2

manPlotRep3 <- plotManhattan(plotRes = manDataRep %>% filter(cat == "Rep Only"), chrCounts,
                             colValues=c("darkorange", "black", "darkorange"), shapeValues=c(15, 18, 15),
                             ylimits=c(0, 12.5), legName="Lung Cancer")
manPlotRep3


# cowplot together
application_plot <- plot_grid(manPlotRep1, manPlotRep2, manPlotRep3,
                              labels=c("A", "B", "C"), nrow=3, label_size=24)
application_plot
#ggsave("Fig4b_decon.pdf", width=9, height=12)




