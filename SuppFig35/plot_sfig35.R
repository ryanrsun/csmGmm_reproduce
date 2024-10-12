# Supp Fig 35
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
s7 <- fread(paste0(outputDir, "reject_bmi_with_overall_neg5_reject_S_1_aID7.txt"))
s8 <- fread(paste0(outputDir, "reject_bmi_with_overall_neg5_reject_S_1_aID8.txt"))
s9 <- fread(paste0(outputDir, "reject_bmi_with_overall_neg5_reject_S_1_aID9.txt"))
s7new <- s7 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))
s8new <- s8 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))
s9new <- s9 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))

# for plotting axes
allZ <- fread(paste0(here::here("Data", "bmi_with_overall.txt"))) %>%
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

# data for manhattan plot - correlated pleiotropy
manDataCor <- rbind(s7new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(cat = "LC,CAD"),
                    s8new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(cat = "LC,BMI"),
                    s9new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(cat = "CAD,BMI")) %>%
  as.data.frame(.) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP)) %>%
  arrange(newLfdr) %>%
  # distinct() keeps the first one, so we arrange first
  distinct(., chrpos, .keep_all = TRUE)

manPlotCor <- plotManhattan(plotRes = manDataCor, chrCounts, colValues=gg_color_hue(3),
                            shapeValues=c(8,16,17), ylimits=c(0, 8), legName="Pleiotropy \n(Within UKB)")
manPlotCor

#ggsave("s8_withinukb.pdf", width=18, height=12)








