# make boxplot of human and mouse microglianess and rotate axis
## SETUP ----
library(xlsx)
library(dplyr)
library(graphics)
library(tidyr)
library(stats)

## OPEN HUMAN MICROGLIANESS + MOUSE MICROGLIANESS AND MERGE ----
human <- read.xlsx('h_microglia_rank.xlsx', sheetName = 'h_microglia_rank')
human <- as.data.frame(human$Micro_C1QC)
mouse <- read.xlsx('m_microglia_ranked.xlsx', sheetName = 'm_microglia_ranked')
mouse <- as.data.frame(mouse$Micro.PVM)

# merge the microglianess values
tgt<- merge(data.frame(human$`human$Micro_C1QC`, row.names=NULL), data.frame(mouse$`mouse$Micro.PVM`, row.names=NULL), by = 0, all = TRUE)[-1]
tgt <- tgt %>% rename('human' = 'human..human.Micro_C1QC.') %>% rename('mouse' = 'mouse..mouse.Micro.PVM.')

# gather the data frame into a long data frame (https://www.youtube.com/watch?v=DiflCZDncOE)
tgt2<- gather(tgt, organism, microglianess, 1:2)

# remove NAs (https://statisticsglobe.com/r-remove-data-frame-rows-with-some-or-all-na)
tgt3 <- na.omit(tgt2)

## MAKE BOXPLOT ----
# (https://www.statology.org/multiple-boxplots-r/)
p1 <- boxplot(microglianess ~ organism, data = tgt3, main = 'microglia in mouse and human patch-seq datasets', xlab = 'microglianess', ylab = 'organism', col = 'coral', horizontal = T)


