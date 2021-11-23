# mgp microglianess vs qcmetrics macrophage (mouse)
library(xlsx)
library(dplyr)
library(ggplot2)

# open qcmetrics output
new <- read.csv("https://raw.githubusercontent.com/keon-arbabi/patch-seq-microglia/main/output/qcMetrics.csv")
new$sample_id = sub("-",".", new$sample_id)
new$sample_id = sub("-",".", new$sample_id)

# open microglianess
old <- read.xlsx("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/patch_seq_microglia/m_microglia_ranked.xlsx", sheetName = "m_microglia_ranked") # did not use here because didn't want to change directory
old <- dplyr::rename(old, sample_id = transcriptomics_sample_id)

# inner_join
tgt <- inner_join(new, old, by = "sample_id")
tgt <- dplyr::select(tgt, sample_id, Macrophage, Micro.PVM)

# make plot
plot1 <- ggplot(
  tgt,
  aes(x = Micro.PVM, y = Macrophage)) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) + 
  geom_point() +
  ggtitle('microglianess vs. macrophage') +
  ylab('qcmetrics macrophage') + xlab ('microglianess') +
  theme_bw()

plot2 <- ggplot(
  tgt,
  aes(x = Micro.PVM, y = Macrophage)) +
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1) + 
  geom_point() +
  ggtitle('microglianess vs. macrophage') +
  ylab('qcmetrics macrophage') + xlab ('microglianess') +
  theme_bw()

plot1
plot2

