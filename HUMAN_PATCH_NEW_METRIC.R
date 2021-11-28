# 3 TOP 25 HUMAN MICROGLIA ----

library(patchwork)
library(tibble)
library(dplyr)
library(janitor)
library(xlsx)
library(ggplot2)
library(cowplot)
library(ggrepel)
theme_set(theme_cowplot())
library(stringr)

## SPLIT THE MOUSE DATA INTO HIGH AND LOW MICROGLIANESS ----
# load mouse dataset and reshape to fit H_microglia_ttype_unranked.xlsx
human_data <- read.csv("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data/20200512_Human_PatchSeq_Release_count.csv")

mouse_data = mouse_data %>% distinct(X, .keep_all = T)
transcript = mouse_data$X
rownames(mouse_data) = rownames(mouse_data) %>% make.names(unique=T)
mouse_data_trans = mouse_data[-1] %>% t() %>% as.data.frame()
colnames(mouse_data_trans) = transcript
mouse_data_trans # now we have the count matrix rotated 90 degrees

human_data = human_data %>% distinct(X, .keep_all = T)
transcript = human_data$X
rownames(human_data) = rownames(human_data) %>% make.names(unique=T)
human_data_trans = human_data[-1] %>% t() %>% as.data.frame()
colnames(human_data_trans) = transcript
human_data_trans = human_data_trans %>% row_to_names(row_number = 1)
human_data_trans # now we have the count matrix rotated 90 degrees

# open qcmetrics output
qc <- read.csv("human_qcMetrics.csv")
qc$sample_id = sub("-",".", qc$sample_id)
qc$sample_id = sub("-",".", qc$sample_id)
qc$sample_id = sub("-",".", qc$sample_id)
microglia_ranked <- qc %>% dplyr::select(c(sample_id, Microglia)) %>% arrange(desc(Microglia))

top25 <- slice_max(microglia_ranked, order_by = Microglia, prop = 0.25) # microglianess cutoff is 704.15163
bottom25 <- slice_min(microglia_ranked, order_by = Microglia, prop = 0.25) # microglianess cutoff is 3059.118

microglia_ranked$group <- "low"
microglia_ranked$group[microglia_ranked$Microglia > 704.15163 & microglia_ranked$Microglia < 3059.118] <- "med"
microglia_ranked$group[microglia_ranked$Microglia >= 3059.118] <- "high"

write.xlsx(microglia_ranked, file = "h_QC_microglia_ranked.xlsx",
           sheetName = "h_QC_microglia_ranked", append = FALSE)

top_25 <- filter(microglia_ranked, group == 'high')
bottom_25 <- filter(microglia_ranked, group == 'low')

top_25 <- dplyr::rename(top_25, transcriptomics_sample_id = sample_id)
bottom_25 <- dplyr::rename(bottom_25, transcriptomics_sample_id = sample_id)

write.xlsx(top_25, file = "h_QC_top_25_microglia.xlsx",
           sheetName = "h_QC_top_25_microglia", append = FALSE)
write.xlsx(bottom_25, file = "h_QC_bottom_25_microglia.xlsx",
           sheetName = "h_QC_bottom_25_microglia", append = FALSE)

# inner_join the transcriptomics_sample_id of top_25 and bottom_25 with human_data_trans, respectively (save as csv)
top_25_id <- dplyr::select(top_25, transcriptomics_sample_id)
bottom_25_id <- dplyr::select(bottom_25, transcriptomics_sample_id)

human_data_trans <- human_data_trans %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'transcriptomics_sample_id')

micro_top25 <- inner_join(human_data_trans, top_25_id, by = 'transcriptomics_sample_id', copy = T)
micro_bottom25 <- inner_join(human_data_trans, bottom_25_id, by = 'transcriptomics_sample_id', copy = T)

write.csv(micro_top25, file = "h_QC_top_25_micro_count_trans.csv")
write.csv(micro_bottom25, file = "h_QC_bottom_25_micro_count_trans.csv")

# reshape micro_top20 and micro_bottom20 and save them as xlsx files (these are the files to work with for the integration)
micro_top25 <- micro_top25 %>% as.data.frame()
micro_top25 = micro_top25 %>% distinct(transcriptomics_sample_id, .keep_all = T)
transcriptomics_sample_id = micro_top25$transcriptomics_sample_id
rownames(micro_top25) = rownames(micro_top25) %>% make.names(unique=T)
human_t = micro_top25[-1] %>% t() %>% as.data.frame()
colnames(human_t) = transcriptomics_sample_id

micro_bottom25 <- micro_bottom25 %>% as.data.frame()
micro_bottom25 = micro_bottom25 %>% distinct(transcriptomics_sample_id, .keep_all = T)
transcriptomics_sample_id = micro_bottom25$transcriptomics_sample_id
rownames(micro_bottom25) = rownames(micro_bottom25) %>% make.names(unique=T)
human_b = micro_bottom25[-1] %>% t() %>% as.data.frame()
colnames(human_b) = transcriptomics_sample_id

write.csv(human_t, file = "h_QC_micro_top25_count.csv")
write.csv(human_b, file = "h_QC_micro_bottom25_count.csv")


## JOIN METADATA WITH TOP25 AND BOTTOM25 ----
human_top <- dplyr::select(micro_top25, transcriptomics_sample_id)
human_bottom <- dplyr::select(micro_bottom25, transcriptomics_sample_id)

# open metadata
meta <- read.csv("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data/20200625_patchseq_metadata_human.csv")

meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)

# inner_join top20 and bottom20 with meta; save as csv
top25sub <- inner_join(huma_top, meta, by = 'transcriptomics_sample_id', copy = T)
bottom25sub <- inner_join(human_bottom, meta, by = 'transcriptomics_sample_id', copy = T)

write.csv(top25sub, file = "h_QC_metadata_top25.csv")
write.csv(bottom25sub, file = "h_QC_metadata_bottom25.csv")


# 5 REGRESSION ----
library(dplyr)
library(ggplot2)
library(xlsx)
library(stats)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(data.table)

meta$s1 <- gsub( "^\\S+\\s+", "", meta$corresponding_AIT2.3.1_alias)
meta$s2 <- gsub( "^\\S+\\s+", "", meta$s1)
meta$subclass <- as.factor(gsub( "^\\S+\\s+", "", meta$s2))
meta <- meta %>% dplyr::select(-c('s1', 's2'))

write.csv(meta, file = "h_QC_metadata_ttype.csv")

micro <- microglia_ranked
micro <- dplyr::rename(micro, transcriptomics_sample_id = sample_id)
data <- inner_join(micro, meta, by = "transcriptomics_sample_id")

data$subclass <- as.factor(data$subclass)
data$donor_id <- as.factor(data$donor_id)

## MAKE REGRESSION ----
# subclass
p1 <- ggplot(
  data,
  aes(x = subclass, y = Microglia, fill = subclass)) + geom_point() +
  ggtitle('Neuronal subclass vs. microglia score\nin the human patch-seq dataset') +
  ylab('Microglia score') + xlab ('Neuronal subclass') +
  theme_bw() + geom_boxplot() + scale_fill_brewer(palette="Spectral")

microglianess <- data$Microglia
subclass <- data$subclass
m1 <- lm(microglianess ~ subclass)
summary(m1) # Multiple R-squared:  0.02053,	Adjusted R-squared:  0.006184 

# cell soma normalized depth
p2 <- ggplot(
  data,
  aes(x = cell_soma_normalized_depth, y = Microglia)) +
  geom_smooth(method = "lm", se = F, color = 'red4') + 
  geom_point(color = "skyblue4") + 
  ggtitle('Cell soma normalized depth vs. microglia score\nin the human patch-seq dataset') + 
  ylab('Microglia score') + xlab ('Cell soma normalized depth') +
  theme_bw()

cell_soma_norm_depth <- data$cell_soma_normalized_depth
m2 <- lm(microglianess ~ cell_soma_norm_depth)
summary(m2) # Multiple R-squared:  0.03725,	Adjusted R-squared:  0.02711

p1 + p2 + plot_layout(ncol = 2, widths=c(4, 3))

# donor ID
cols <- 48
mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(cols)

p3 <- ggplot(
  data,
  aes(x = donor_id, y = Microglia, fill = donor_id)) + geom_point() +
  ggtitle('Donor ID vs. microglia score\nin the human patch-seq dataset') +
  ylab('Microglia score') + xlab ('Donor ID') +
  theme_bw() + geom_boxplot() + scale_fill_manual(values = mycolors) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

donor_id <- data$donor_id
m3 <- lm(microglianess ~ donor_id)
summary(m3) # Multiple R-squared:  0.3715,	Adjusted R-squared:  0.2431 

## IBA1 ----
# https://github.com/AllenInstitute/patchseq_human_L23/blob/master/data/pathology_scoring.csv
pathology <- fread('https://raw.githubusercontent.com/AllenInstitute/patchseq_human_L23/master/data/pathology_scoring.csv') %>% as.data.frame()
pathology <- dplyr::select(pathology, donor_id, 'Iba1 Cortex')

write.xlsx(pathology, file = "h_IBA1.xlsx", sheetName = "h_IBA1", append = FALSE)
# then adjust donor ID in excel :(

#after editing on excel, continue
pathology <- read.xlsx('h_IBA1.xlsx', sheetName = 'h_IBA1')

join <- inner_join(data, pathology, by = "donor_id", copy = FALSE)

cpm <- read.csv('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data/20200512_Human_PatchSeq_Release_cpm.csv')
cpm <- cpm %>% as.data.frame()
cpm = cpm %>% distinct(X, .keep_all = T)
X = cpm$X
rownames(cpm) = rownames(cpm) %>% make.names(unique=T)
cpm_trans = cpm[-1] %>% t() %>% as.data.frame()
colnames(cpm_trans) = X

cpm_aif1 <- cpm_trans %>% tibble::rownames_to_column(var = 'transcriptomics_sample_id') %>% dplyr::select(c('transcriptomics_sample_id', 'AIF1'))

iba.aif.micro <- inner_join(join, cpm_aif1, by = "transcriptomics_sample_id")
write.xlsx(iba.aif.micro, file = "h_IBA1_AIF1_microglianess.xlsx",
           sheetName = "h_IBA1_AIF1_microglianess", append = FALSE)

# make plot
p4 <- ggplot(
  iba.aif.micro,
  aes(x = Iba1.Cortex, y = Microglia)) +
  geom_smooth(method = "lm", se = F, color = 'red4') + 
  geom_point(color = "skyblue4") + 
  ggtitle('IBA1 protein level vs. microglia score\nin the human patch-seq dataset') + 
  ylab('Microglia score') + xlab ('IBA1 level') +
  theme_bw()

iba1 <- iba.aif.micro$Iba1.Cortex
micro <- iba.aif.micro$Microglia
m4 <- lm(micro ~ iba1)
summary(m4) # Multiple R-squared:  0.0001673,	Adjusted R-squared:  -0.003577

p5 <- ggplot(
  iba.aif.micro,
  aes(x = Iba1.Cortex, y = AIF1)) +
  geom_smooth(method = "lm", se = F, color = 'red4') + 
  geom_point(color = "skyblue4") + 
  ggtitle('IBA1 protein level vs. AIF1 expression\nin the human patch-seq dataset') + 
  ylab('AIF1 expression') + xlab ('IBA1 level') +
  theme_bw()

aif1 <- iba.aif.micro$AIF1
m5 <- lm(iba1 ~ aif1)
summary(m5) # Multiple R-squared:  0.0008681,	Adjusted R-squared:  -0.002874

p3 + p4 + plot_layout(ncol = 2, widths=c(4, 3))
p3 + p4

