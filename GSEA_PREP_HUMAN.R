# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# make normalized count for human data for GSEA analysis
# make expression dataset and phenotype labels files
# create own gene sets from microglia subtype markers

library(tibble)
library(dplyr)
library(DESeq2)

#### MATCH THE METADATA AND COUNTS DATA ####
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data')
data <- read.csv('20200512_Human_PatchSeq_Release_count.csv')

top25_meta <- read.csv('patchseq_metadata_human_ttype_top25.csv')
top25_meta$microglianess <- 'high'
bottom25_meta <- read.csv('patchseq_metadata_human_ttype_bottom25.csv')
bottom25_meta$microglianess <- 'low'
meta <- bind_rows(top25_meta, bottom25_meta)
meta <- as.data.frame(meta[-1])
meta <- dplyr::rename(meta, 'transcriptomics.sample.id'= 'transcriptomics_sample_id')

# no need but just in case
 # meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
 # meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
 # meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)

# the transcriptomics_sample_id don't match, so rotate the count matrix and inner_join with meta
data <- data %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'transcript')

data <- data %>% distinct(transcript, .keep_all = T)
transcript = data$transcript
rownames(data) = rownames(data) %>% make.names(unique=T)
data_trans = data[-1] %>% t() %>% as.data.frame() # integers became string :(
colnames(data_trans) = transcript
data_trans = data_trans %>% janitor::row_to_names(row_number = 1)
data_trans[] <- lapply(data_trans, function(x) as.numeric(as.character(x)))
View(data_trans) # now we have the count matrix rotated 90 degrees

data_trans <- tibble::rownames_to_column(data_trans, var = 'transcriptomics.sample.id')
data_id <- as.data.frame(data_trans$transcriptomics.sample.id)
data_id <- dplyr::rename(data_id, 'transcriptomics.sample.id'= 'data_trans$transcriptomics.sample.id')

meta_joined <- inner_join(data_id, meta, by = 'transcriptomics.sample.id')
meta_joined$rowname <- meta_joined$transcriptomics.sample.id
meta_joined <- column_to_rownames(meta_joined, var = 'rowname')

data_trans <- inner_join(data_trans, meta_joined, by = 'transcriptomics.sample.id')
data_trans <- select(data_trans, -c(subclass, microglianess))

# reshape data_trans
data_trans <- data_trans %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'transcript')

data_trans <- data_trans %>% distinct(transcript, .keep_all = T)
transcript = data_trans$transcript
rownames(data_trans) = rownames(data_trans) %>% make.names(unique=T)
data2 = data_trans[-1] %>% t() %>% as.data.frame()
colnames(data2) = transcript
data2 = data2 %>% janitor::row_to_names(row_number = 1)
data2[] <- lapply(data2, function(x) as.numeric(as.character(x)))
View(data2)

all(colnames(data2) %in% rownames(meta_joined))
all(colnames(data2) == rownames(meta_joined))


#### CREATE DESEq2 OBJECT ####
## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data2, colData = meta_joined, design = ~ microglianess)
View(counts(dds))


#### GENERATE THE MOV10 NORMALIZED COUNTS ####
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="berg_normalized_counts.txt", sep="\t", quote=F, col.names=NA)


#### MAKE PHENOTYPE LABELS FILES ####
# skipped expression dataset because did it on Excel
# also skipped gene set because did it on Excel
# for the phenotype labels files, take microglianess row in meta_joined and reshape (already checked to make sure that the order is the same as expression dataset)
group <- as.data.frame(meta_joined$microglianess)
group_trans <- group %>% t() %>% as.data.frame()

write.xlsx(group_trans, file = "phenotype_labels.xlsx",
           sheetName = "phenotype_labels", append = FALSE)