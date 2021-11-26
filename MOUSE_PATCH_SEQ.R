## MOUSE PATCH-SEQ ##
# 1 SEURAT TUTORIAL + DIMENSIONALITY REDUCTION ----
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(dplyr)
library(Seurat)
library(patchwork)

## SET UP THE SEURAT OBJECT ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/mouse_patch_seq/data')
mouse_data <- read.csv("20200513_Mouse_PatchSeq_Release_count.v2.csv")
rownames(mouse_data) <- mouse_data$X
mouse_seurat <- CreateSeuratObject(counts = mouse_data, project = "mouse", min.cells = 3, min.features = 200) # min.cells = include features detected in at least 3 cells; min.features = include cells where at least 200 features are detected
mouse_seurat # 39182 features (transcripts) across 4435 samples (cells) within 1 assay

## STANDARD PRE-PROCESSING WORKFLOW ----
# visualize QC metrics as a violin plot
VlnPlot(mouse_seurat, features = c("nFeature_RNA", "nCount_RNA"), cols = "lightsteelblue3", ncol = 2)

# feature scatter
FeatureScatter(mouse_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

mouse_seurat <- subset(mouse_seurat, subset = nFeature_RNA > 3000 & nFeature_RNA < 15000)

## NORMALIZING THE DATA ----
mouse_seurat <- NormalizeData(mouse_seurat, normalization.method = "LogNormalize", scale.factor = 1000000)

## IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION) ----
# calculate a subset of features that exhibit high cell-to-cell variation in the dataset
mouse_seurat <- FindVariableFeatures(mouse_seurat, selection.method = "vst", nfeatures = 5000)

# identify the 15 most highly variable genes
top15 <- head(VariableFeatures(mouse_seurat), 15)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mouse_seurat)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
plot1 + plot2

## SCALING THE DATA ----
# mean is 0, variance is 1
all.genes <- rownames(mouse_seurat)
mouse_seurat <- ScaleData(mouse_seurat, features = all.genes)


## PERFORM LINEAR DIMENSIONAL REDUCTION ----
mouse_seurat <- RunPCA(mouse_seurat, features = VariableFeatures(object = mouse_seurat))
DimHeatmap(mouse_seurat, dims = 1:10, cells = 500, balanced = TRUE)

## DETERMINE THE ‘DIMENSIONALITY’ OF THE DATASET ----
mouse_seurat <- JackStraw(mouse_seurat, num.replicate = 100)
mouse_seurat <- ScoreJackStraw(mouse_seurat, dims = 1:20)
JackStrawPlot(mouse_seurat, dims = 1:20)
ElbowPlot(mouse_seurat)

## CLUSTER THE CELLS ----
mouse_seurat <- FindNeighbors(mouse_seurat, dims = 1:20)
mouse_seurat <- FindClusters(mouse_seurat, resolution = 0.5)

## RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/tSNE) ----
reticulate::py_install(packages = 'umap-learn')
mouse_seurat <- RunUMAP(mouse_seurat, dims = 1:20)

m_seurat_plot = DimPlot(mouse_seurat, reduction = "umap", label = FALSE)

# to save
saveRDS(mouse_seurat, file = "/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/patch_seq_microglia/m_seurat.rds")

# to read if needed
mouse_seurat <- readRDS('m_seurat.rds')

## FINDING DIFFERENTIALLY EXPRESSED FEATURES (CLUSTER BIOMARKERS) ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')

# find markers for every cluster compared to all remaining clusters, report positive markers
mouse.markers <- FindAllMarkers(mouse_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mouse.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.xlsx(mouse.markers, file = "m_all_markers.xlsx",
           sheetName = "m_all_markers", append = FALSE)

# find markers distinguishing cluster x from all other clusters
cluster0.markers <- FindMarkers(mouse_seurat, ident.1 = 0, min.pct = 0.25)
cluster0.markers.genes <- rownames_to_column(cluster0.markers)
head(cluster0.markers.genes[1], n = 100)
write.xlsx(cluster0.markers.genes, file = "m_cluster0_markers.xlsx",
           sheetName = "m_cluster0_markers", append = FALSE)

cluster1.markers <- FindMarkers(mouse_seurat, ident.1 = 1, min.pct = 0.25)
cluster1.markers.genes <- rownames_to_column(cluster1.markers)
write.xlsx(cluster1.markers.genes, file = "m_cluster1_markers.xlsx",
           sheetName = "m_cluster1_markers", append = FALSE)

cluster2.markers <- FindMarkers(mouse_seurat, ident.1 = 2, min.pct = 0.25)
cluster2.markers.genes <- rownames_to_column(cluster2.markers)
write.xlsx(cluster2.markers.genes, file = "m_cluster2_markers.xlsx",
           sheetName = "m_cluster2_markers", append = FALSE)

cluster3.markers <- FindMarkers(mouse_seurat, ident.1 = 3, min.pct = 0.25)
cluster3.markers.genes <- rownames_to_column(cluster3.markers)
write.xlsx(cluster3.markers.genes, file = "m_cluster3_markers.xlsx",
           sheetName = "m_cluster3_markers", append = FALSE)

cluster4.markers <- FindMarkers(mouse_seurat, ident.1 = 4, min.pct = 0.25)
cluster4.markers.genes <- rownames_to_column(cluster4.markers)
write.xlsx(cluster4.markers.genes, file = "m_cluster4_markers.xlsx",
           sheetName = "m_cluster4_markers", append = FALSE)

cluster5.markers <- FindMarkers(mouse_seurat, ident.1 = 5, min.pct = 0.25)
cluster5.markers.genes <- rownames_to_column(cluster5.markers)
write.xlsx(cluster5.markers.genes, file = "m_cluster5_markers.xlsx",
           sheetName = "m_cluster5_markers", append = FALSE)

cluster6.markers <- FindMarkers(mouse_seurat, ident.1 = 6, min.pct = 0.25)
cluster6.markers.genes <- rownames_to_column(cluster6.markers)
write.xlsx(cluster6.markers.genes, file = "m_cluster6_markers.xlsx",
           sheetName = "m_cluster6_markers", append = FALSE)

cluster7.markers <- FindMarkers(mouse_seurat, ident.1 = 7, min.pct = 0.25)
cluster7.markers.genes <- rownames_to_column(cluster7.markers)
write.xlsx(cluster7.markers.genes, file = "m_cluster7_markers.xlsx",
           sheetName = "m_cluster7_markers", append = FALSE)

cluster8.markers <- FindMarkers(mouse_seurat, ident.1 = 8, min.pct = 0.25)
cluster8.markers.genes <- rownames_to_column(cluster8.markers)
write.xlsx(cluster8.markers.genes, file = "m_cluster8_markers.xlsx",
           sheetName = "m_cluster8_markers", append = FALSE)

cluster9.markers <- FindMarkers(mouse_seurat, ident.1 = 9, min.pct = 0.25)
cluster9.markers.genes <- rownames_to_column(cluster9.markers)
write.xlsx(cluster9.markers.genes, file = "m_cluster9_markers.xlsx",
           sheetName = "m_cluster9_markers", append = FALSE)

cluster10.markers <- FindMarkers(mouse_seurat, ident.1 = 10, min.pct = 0.25)
cluster10.markers.genes <- rownames_to_column(cluster10.markers)
write.xlsx(cluster10.markers.genes, file = "m_cluster10_markers.xlsx",
           sheetName = "m_cluster10_markers", append = FALSE)

cluster11.markers <- FindMarkers(mouse_seurat, ident.1 = 11, min.pct = 0.25)
cluster11.markers.genes <- rownames_to_column(cluster11.markers)
write.xlsx(cluster11.markers.genes, file = "m_cluster11_markers.xlsx",
           sheetName = "m_cluster11_markers", append = FALSE)

cluster12.markers <- FindMarkers(mouse_seurat, ident.1 = 12, min.pct = 0.25)
cluster12.markers.genes <- rownames_to_column(cluster12.markers)
write.xlsx(cluster12.markers.genes, file = "m_cluster12_markers.xlsx",
           sheetName = "m_cluster12_markers", append = FALSE)

cluster13.markers <- FindMarkers(mouse_seurat, ident.1 = 13, min.pct = 0.25)
cluster13.markers.genes <- rownames_to_column(cluster13.markers)
write.xlsx(cluster13.markers.genes, file = "m_cluster13_markers.xlsx",
           sheetName = "m_cluster13_markers", append = FALSE)

cluster14.markers <- FindMarkers(mouse_seurat, ident.1 = 14, min.pct = 0.25)
cluster14.markers.genes <- rownames_to_column(cluster14.markers)
write.xlsx(cluster14.markers.genes, file = "m_cluster14_markers.xlsx",
           sheetName = "m_cluster14_markers", append = FALSE)

# visualization
mouse_violin <- VlnPlot(mouse_seurat, features = c("Cd74", "Tspo", "C1qc", "Aif1"), ncol = 2)
mouse_violin2 <- VlnPlot(mouse_seurat, features = c("C1qc", "Aif1","Cd74", "Tmem119", "Tspo"), ncol = 2)
FeaturePlot(mouse_seurat, features = c("C1qc", "C1qa", "Aif1","Cd74", "Spi1", "Tspo"))

m_seurat_plot + mouse_violin

# can plot raw counts as well
VlnPlot(mouse_seurat, features = c("C1qc", "C1qa"), slot = "counts", log = TRUE)


# DoHeatmap() generates an expression heatmap for given cells and features
# plot top 10 markers (or all markers if less than 10) for each cluster
top10 <- mouse.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(mouse_seurat, features = top10$gene, group.bar = TRUE, angle = 0) + scale_fill_gradientn(colors = c("lightsteelblue3", "black", "salmon2")) + theme(text = element_text(size = 3))

# did not label the clusters here


# 2 MGP + LABONTE ----
library(here)
library(tidyverse)
library(matrixStats)
library(cowplot)
library(broom)
ggplot2::theme_set(theme_cowplot())
library(edgeR)
library(markerGeneProfile) 
library(RCurl)
library(xlsx)


## READ MOUSE CELL TYPE SPECIFIC MARKERS FROM EXTERNAL CSV FILE ----
# load mouse markers from https://github.com/keon-arbabi/patch-seq-microglia/blob/main/shared/AIBS-mouse-VIS-smartseq-markers.csv
mouse_markers = read.csv("https://raw.githubusercontent.com/keon-arbabi/patch-seq-microglia/main/shared/AIBS-mouse-VIS-smartseq-markers.csv?token=AUB7ZH65DEWQKSJT5Q5NPA3BCHRHM")
cell_types = mouse_markers$subclass %>% unique()
knitr::kable(head(mouse_markers))


## PROCESS MARKERS DATA FRAME FOR MGP FUNCTION ----
marker_list <- lapply(cell_types, function(cell_types){
  return(mouse_markers %>% filter(subclass_label == cell_types) %>% pull(gene) %>% unlist() %>% str_to_title())
})
names(marker_list) <- cell_types
print(cell_types)


## LIST NAMES OF MARKER GENES OF THE 'MICROGLIA' CELL TYPE ----
print(marker_list$'Micro-PVM')
micro_markers_save <- as.data.frame(marker_list$'Micro-PVM')
write.xlsx(micro_markers_save, file = "m_microglia_markers.xlsx",
           sheetName = "m_microglia_markers", append = FALSE)


## READ AND PROCESS TRANSCRIPTOMIC DATA OF MOUSE PATCH SEQ ----
# Read in the metadata and gene expression count matrix
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/mouse_patch_seq/data")
mouse_meta <- read.csv('20200711_patchseq_metadata_mouse.csv') %>%
  dplyr::select(transcriptomics_sample_id, 
         t_type_label)
mouse_expr <- read.csv("20200513_Mouse_PatchSeq_Release_count.v2.csv") %>%
  tibble::column_to_rownames("X") %>% as.data.frame()
knitr::kable(mouse_expr[1:5, 1:5])


## NORMALIZE THE GENE EXPRESSION DATASET INTO CPM ----
mouse_cpm <- edgeR::cpm(mouse_expr, log = F, prior.count = 0)
knitr::kable(mouse_cpm[1:5, 1:5])

# or directly open the CPM file
# setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/mouse_patch_seq/data")
# mouse_cpm <- read.csv("20200513_Mouse_PatchSeq_Release_cpm.v2.csv") %>% tibble::column_to_rownames("X") %>% as.data.frame()


## PROCESS EXPRESSION MATRIX PRIOR TO MGP ANALYSIS ----
use_gene_list <- mouse_cpm %>% rownames() %>% unlist %>% unique()
gene_mat <- mouse_cpm[use_gene_list, ]


## CALCULATE GENES WITH VERY LOW SD AND REMOVE THEM FROM THE EXPRESSION MATRIX ----
gene_sds <- rowSds(gene_mat %>% as.matrix(), na.rm = T) 
gene_mat <- gene_mat[gene_sds > .1, ]


## RESHAPE GENE EXPRESSION MATRIX ----
gene_mat <- gene_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'X')

gene_mat = gene_mat %>% distinct(X, .keep_all = T)
X = gene_mat$X
rownames(gene_mat) = rownames(gene_mat) %>% make.names(unique=T)
gene_mat_trans = gene_mat[-1] %>% t() %>% as.data.frame()
colnames(gene_mat_trans) = X


## MERGE SAMPLE METADATA MATRIX WITH GENE EXPRESSION MATRIX ----
gene_mat_trans <- gene_mat_trans %>% 
  tibble::rownames_to_column(var = 'transcriptomics_sample_id')

mouse_meta$transcriptomics_sample_id = sub("-",".", mouse_meta$transcriptomics_sample_id)
mouse_meta$transcriptomics_sample_id = sub("-",".", mouse_meta$transcriptomics_sample_id)

# merge gene expression and meta data frames
gene_mat_comb <- inner_join(mouse_meta, gene_mat_trans, by = 'transcriptomics_sample_id')


## RUN CELL TYPE PROPORTION ESTIMATION ----
# run MGP analysis
estimations <-  mgpEstimate(
  exprData=gene_mat,
  genes=marker_list,
  geneColName='X',
  outlierSampleRemove=F, # should outlier samples removed. this is done using boxplot stats.
  geneTransform = NULL, # this is the default option for geneTransform
  groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus=FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority=TRUE)


## MERGE RELATIVE CELL TYPE PROPORTIONS PER CELL TYPE WITH SAMPLE METADATA ----
mgp_estimates <- as.data.frame(estimations$estimates) %>%
  tibble::rownames_to_column(var = 'transcriptomics_sample_id')

mgp_df <- inner_join(mouse_meta, mgp_estimates, by = "transcriptomics_sample_id")
knitr::kable(mgp_df[1:5, 1:12])


## MAKE FILES CONTAINING CELLS AND THEIR 'MICROGLIANESS' ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
# ranked
microglia_ranked <- mgp_df %>% select(c(transcriptomics_sample_id, Micro.PVM)) %>% arrange(desc(Micro.PVM))
ttype <- mouse_meta %>% select(c(transcriptomics_sample_id, t_type_label))
microglia_ttype <- inner_join(microglia_ranked, ttype, by = "transcriptomics_sample_id")
top20 <- slice_max(microglia_ttype, order_by = Micro.PVM, prop = 0.20) # microglianess cutoff is 0.3353369
bottom20 <- slice_min(microglia_ttype, order_by = Micro.PVM, prop = 0.20) # microglianess cutoff is -3.676738

microglia_ttype$group <- "low"
microglia_ttype$group[microglia_ttype$Micro.PVM > -3.676738 & microglia_ttype$Micro.PVM < 0.3353369] <- "med"
microglia_ttype$group[microglia_ttype$Micro.PVM >= 0.3353369] <- "high"

write.xlsx(microglia_ttype, file = "m_microglia_ttype_ranked.xlsx",
           sheetName = "m_microglia_ttype_ranked", append = FALSE)

# unranked
microglia_unranked <- mgp_df %>% select(c(transcriptomics_sample_id, Micro.PVM))
order <- as.data.frame(microglia_unranked$transcriptomics_sample_id)
order <- rename(order, transcriptomics_sample_id = 'microglia_unranked$transcriptomics_sample_id')
microglia_ttype_unranked <- inner_join(order, microglia_ttype, by = "transcriptomics_sample_id")
write.xlsx(microglia_ttype_unranked, file = "m_microglia_ttype_unranked.xlsx",
           sheetName = "m_microglia_ttype_unranked", append = FALSE)

# if needed, take column of mouse_cpm (store as new data) and inner join to the first column of microglia_ranked to get a ranked data frame
microglia_ranked_full <- inner_join(microglia_ranked, gene_mat_trans, by = "transcriptomics_sample_id")
write.xlsx(microglia_ranked_full[, 1:2], file = "m_microglia_ranked.xlsx",
           sheetName = "m_microglia_ranked", append = FALSE)

# rm(list = ls()) to remove everything


# 3 SEURAT INTEGRATION ----
# https://satijalab.org/seurat/articles/integration_introduction.html
# https://rpubs.com/mathetal/integratedanalysis

library(Seurat)
library(SeuratData)
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
# load mouse dataset and reshape to fit m_microglia_ttype_unranked.xlsx
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/mouse_patch_seq/data')
mouse_data <- read.csv("20200513_Mouse_PatchSeq_Release_count.v2.csv")

mouse_data <- mouse_data %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'transcript')

mouse_data = mouse_data %>% distinct(transcript, .keep_all = T)
transcript = mouse_data$transcript
rownames(mouse_data) = rownames(mouse_data) %>% make.names(unique=T)
mouse_data_trans = mouse_data[-1] %>% t() %>% as.data.frame()
colnames(mouse_data_trans) = transcript
mouse_data_trans = mouse_data_trans %>% row_to_names(row_number = 1)
mouse_data_trans # now we have the count matrix rotated 90 degrees

# load m_microglia_ttype_unranked.xlsx
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
micro <- read.xlsx('m_microglia_ttype_unranked.xlsx', sheetName = 'm_microglia_ttype_unranked')
micro <- select(micro, -c(NA.)) # removed first column of numbers

# create new dataframes: top_20 and bottom_20 (and save as xlsx just in case)
top_20 <- filter(micro, group == 'high')
bottom_20 <- filter(micro, group == 'low')

write.xlsx(top_20, file = "m_top_20_microglia.xlsx",
           sheetName = "m_top_20_microglia", append = FALSE)
write.xlsx(bottom_20, file = "m_bottom_20_microglia.xlsx",
           sheetName = "m_bottom_20_microglia", append = FALSE)

# inner_join the transcriptomics_sample_id of top_20 and bottom_20 with mouse_data_trans, respectively (save as csv)
top_20 <- as.data.frame(top_20[, 1])
top_20 <- rename(top_20, 'transcriptomics_sample_id' = 'top_20[, 1]')

bottom_20 <- as.data.frame(bottom_20[, 1])
bottom_20 <- rename(bottom_20, 'transcriptomics_sample_id' = 'bottom_20[, 1]')

mouse_data_trans <- mouse_data_trans %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'transcriptomics_sample_id')

micro_top20 <- inner_join(mouse_data_trans, top_20, by = 'transcriptomics_sample_id', copy = T)
micro_bottom20 <- inner_join(mouse_data_trans, bottom_20, by = 'transcriptomics_sample_id', copy = T)

write.csv(micro_top20, file = "./ignored/m_top_20_micro_count_trans.csv")
write.csv(micro_bottom20, file = "./ignored/m_bottom_20_micro_count_trans.csv")

# reshape micro_top20 and micro_bottom20 and save them as xlsx files (these are the files to work with for the integration)
micro_top20 <- micro_top20 %>% as.data.frame()
micro_top20 = micro_top20 %>% distinct(transcriptomics_sample_id, .keep_all = T)
transcriptomics_sample_id = micro_top20$transcriptomics_sample_id
rownames(micro_top20) = rownames(micro_top20) %>% make.names(unique=T)
mouse_top = micro_top20[-1] %>% t() %>% as.data.frame()
colnames(mouse_top) = transcriptomics_sample_id

micro_bottom20 <- micro_bottom20 %>% as.data.frame()
micro_bottom20 = micro_bottom20 %>% distinct(transcriptomics_sample_id, .keep_all = T)
transcriptomics_sample_id = micro_bottom20$transcriptomics_sample_id
rownames(micro_bottom20) = rownames(micro_bottom20) %>% make.names(unique=T)
mouse_bottom = micro_bottom20[-1] %>% t() %>% as.data.frame()
colnames(mouse_bottom) = transcriptomics_sample_id

write.csv(mouse_top, file = "./ignored/m_micro_top20_count.csv")
write.csv(mouse_bottom, file = "./ignored/m_micro_bottom20_count.csv")


## JOIN METADATA + TTYPE WITH TOP20 AND BOTTOM20 ----
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final")
mouse_top <- read.csv("m_top_20_micro_count_trans.csv")
mouse_bottom <- read.csv("m_bottom_20_micro_count_trans.csv")
mouse_top <- as.data.frame(mouse_top[, 2])
mouse_bottom <- as.data.frame(mouse_bottom[, 2])
mouse_top <- rename(mouse_top, "transcriptomics_sample_id" = "mouse_top[, 2]")
mouse_bottom <- rename(mouse_bottom, "transcriptomics_sample_id" = "mouse_bottom[, 2]")

# open metadata and add ttype
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/mouse_patch_seq/data')
meta <- read.csv("20200625_patchseq_metadata_mouse.csv")
meta$subclass <- gsub( "\\s.*", "", meta$corresponding_AIT2.3.1_alias) # take everything before first space

meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)

setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
write.csv(meta, file = "m_metadata_ttype.csv")

meta <- meta %>% dplyr::select(transcriptomics_sample_id, subclass)

# inner_join top20 and bottom20 with meta; save as csv
top20sub <- inner_join(mouse_top, meta, by = 'transcriptomics_sample_id', copy = T)
bottom20sub <- inner_join(mouse_bottom, meta, by = 'transcriptomics_sample_id', copy = T)

write.csv(top20sub, file = "m_metadata_ttype_top20.csv")
write.csv(bottom20sub, file = "m_metadata_ttype_bottom20.csv")

# remove everything to clear up space: rm(list = ls())


## SET UP HIGH MICROGLIANESS SEURAT OBJECT AND JOIN WITH METADATA ----
# open m_micro_top20_count.csv 
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
top <- read.csv('m_micro_top20_count.csv')
rownames(top) <- top$X

top_seurat <- CreateSeuratObject(counts = top, project = "high", min.cells = 3, min.features = 200)

# add ttype/subclass
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final")
top20_meta <- read.csv("m_metadata_ttype_top20.csv")
subclass_f <- as.factor(top20_meta$subclass)
top_seurat <- AddMetaData(top_seurat, col.name = "ttype", metadata = subclass_f)
# don't do Idents(top_seurat) <- "ttype", because it'll change "high" to "ttypes"

# normalization
top_seurat <- NormalizeData(top_seurat, normalization.method = "LogNormalize", scale.factor = 1000000)

# highly variable features
top_seurat <- FindVariableFeatures(top_seurat, selection.method = "vst", nfeatures = 5000)

# scaling the data
all.genes <- rownames(top_seurat)
top_seurat <- ScaleData(top_seurat, features = all.genes)

# add metadata
top_seurat@meta.data[, "microglianess"] <- "high"

# linear dimensional reduction
top_seurat <- RunPCA(top_seurat, features = VariableFeatures(object = top_seurat))

# determine dimensionality of dataset
top_seurat <- JackStraw(top_seurat, num.replicate = 100)
top_seurat <- ScoreJackStraw(top_seurat, dims = 1:20)
JackStrawPlot(top_seurat, dims = 1:20)
ElbowPlot(top_seurat)

# cluster the cells
top_seurat <- FindNeighbors(top_seurat, dims = 1:20)
top_seurat <- FindClusters(top_seurat, resolution = 0.5)
head(Idents(top_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
top_seurat <- RunUMAP(top_seurat, dims = 1:20)

DimPlot (top_seurat, reduction = "umap", label = TRUE)
DimPlot(top_seurat, reduction = "umap", pt.size = 1, group.by = "microglianess", label = TRUE)


## SET UP LOW MICROGLIANESS SEURAT OBJECT ----
# open m_micro_bottom20_count.csv 
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
bottom <- read.csv('m_micro_bottom20_count.csv')
rownames(bottom) <- bottom$X

bottom_seurat <- CreateSeuratObject(counts = bottom, project = "low", min.cells = 3, min.features = 200)

# add ttype/subclass
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final")
bottom20_meta <- read.csv("m_metadata_ttype_bottom20.csv")
subclass_f2 <- as.factor(bottom20_meta$subclass)
bottom_seurat <- AddMetaData(bottom_seurat, col.name = "ttype", metadata = subclass_f2)

# normalization
bottom_seurat <- NormalizeData(bottom_seurat, normalization.method = "LogNormalize", scale.factor = 1000000)

# highly variable features
bottom_seurat <- FindVariableFeatures(bottom_seurat, selection.method = "vst", nfeatures = 5000)

# scaling the data
all.genes <- rownames(bottom_seurat)
bottom_seurat <- ScaleData(bottom_seurat, features = all.genes)

# add metadata
bottom_seurat@meta.data[, "microglianess"] <- "low"

# linear dimensional reduction
bottom_seurat <- RunPCA(bottom_seurat, features = VariableFeatures(object = bottom_seurat))

# determine dimensionality of dataset
bottom_seurat <- JackStraw(bottom_seurat, num.replicate = 100)
bottom_seurat <- ScoreJackStraw(bottom_seurat, dims = 1:20)
JackStrawPlot(bottom_seurat, dims = 1:20)
ElbowPlot(bottom_seurat)

# cluster the cells
bottom_seurat <- FindNeighbors(bottom_seurat, dims = 1:20)
bottom_seurat <- FindClusters(bottom_seurat, resolution = 0.5)
head(Idents(bottom_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
bottom_seurat <- RunUMAP(bottom_seurat, dims = 1:20)

DimPlot (bottom_seurat, reduction = "umap", label = TRUE)
DimPlot(bottom_seurat, reduction = "umap", pt.size = 1, group.by = "microglianess", label = TRUE)


## COMBINE AND SPLIT SEURAT OBJECTS ----
top_bottom_combined = merge(top_seurat, y = bottom_seurat, add.cell.ids = c("high", "low"), project = "microglianess")

# split the dataset into a list of two seurat objects (high and low)
microglia.list <- SplitObject(top_bottom_combined, split.by = "microglianess")
reference.list <- microglia.list[c("high", "low")]


## PERFORM INTEGRATION ----
microglia.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)

# create 'integrated' data assay
microglia.combined <- IntegrateData(anchorset = microglia.anchors, dims = 1:20)


## PERFORM AN INTEGRATED ANALYSIS ----
DefaultAssay(microglia.combined) <- "integrated"

# run the standard workflow for visualization and clustering
microglia.combined <- ScaleData(microglia.combined, verbose = FALSE)
microglia.combined <- RunPCA(microglia.combined, npcs = 20, verbose = FALSE)
microglia.combined <- RunUMAP(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindNeighbors(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindClusters(microglia.combined, resolution = 0.5)

saveRDS(microglia.combined, file = "/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final/./ignored/m_microglia_combined.rds")

# if want to read it
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
microglia.combined <- readRDS('m_microglia_combined.rds')

# visualization
p1 <- DimPlot(microglia.combined, reduction = "umap", group.by = "microglianess")
p2 <- DimPlot(microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(microglia.combined, reduction = "umap", split.by = "microglianess")
p1 + p2
p3

# save figure
ggsave(
  "m_microglia_combined_umap.jpeg",
  plot = p1,
  device = jpeg,
  path = "/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final",
  width = 15,
  height = 10,
  units = "cm",
)

# more visualization
DefaultAssay(microglia.combined) <- "RNA"
FeaturePlot(microglia.combined, features = c("C1qc", "C1qa", "Aif1", "Spi1"), min.cutoff = "q9")
markers.to.plot <- c("C1qc", "C1qa", "Aif1", "Spi1", "Tspo", "Cd74", "Adap2", "C3ar1", "Ccr5", "Cd14")
DotPlot(microglia.combined, features = markers.to.plot, cols = c("mediumturquoise", "coral2"), dot.scale = 8, split.by = "microglianess") + RotatedAxis()


## IDENTIFY DIFERENTIAL EXPRESSED GENES ACROSS CONDITIONS ----
# make high vs. low expression (for graphs)
DefaultAssay(microglia.combined) <- "integrated"
Idents(microglia.combined) <- "microglianess" # make sure that the assay is 'integrated'
avg.microglia.combined <- as.data.frame(log1p(AverageExpression(microglia.combined, verbose = FALSE)$RNA))
avg.microglia.combined$gene <- rownames(avg.microglia.combined)

# open external mouse microglia markers (for graphs)
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/kcni_summer")
micro_markers <- read.xlsx("microglia_markers_mouse.xlsx", sheetName = "microglia_markers")
micro_markers <- as.data.frame(micro_markers[, 2])
micro_markers <- dplyr::rename(micro_markers, 'marker.gene' = "micro_markers[, 2]")

# open Olah's markers (for graphs)
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/olah_suppl_data')
m1 <- read.xlsx('m1_de.xlsx', sheetName = 'm1_de')
m1$gene <- str_to_title(m1$gene)
m2 <- read.xlsx('m2_de.xlsx', sheetName = 'm2_de')
m2$gene <- str_to_title(m2$gene)
m3 <- read.xlsx('m3_de.xlsx', sheetName = 'm3_de')
m3$gene <- str_to_title(m3$gene)
m4 <- read.xlsx('m4_de.xlsx', sheetName = 'm4_de')
m4$gene <- str_to_title(m4$gene)
m5 <- read.xlsx('m5_de.xlsx', sheetName = 'm5_de')
m5$gene <- str_to_title(m5$gene)
m6 <- read.xlsx('m6_de.xlsx', sheetName = 'm6_de')
m6$gene <- str_to_title(m6$gene)
m7 <- read.xlsx('m7_de.xlsx', sheetName = 'm7_de')
m7$gene <- str_to_title(m7$gene)
m8 <- read.xlsx('m8_de.xlsx', sheetName = 'm8_de')
m8$gene <- str_to_title(m8$gene)
m9 <- read.xlsx('m9_de.xlsx', sheetName = 'm9_de')
m9$gene <- str_to_title(m9$gene)

options(ggrepel.max.overlaps = Inf)

# plots
p_basic <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("all top and bottom 20%") + geom_text(aes(label=ifelse(high/low > 2 & high > 3, as.character(gene),'')),hjust=0,vjust=0)
p_basic # only labelled high > 3 to prevent laptop dying

p_colored <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("all top and bottom 20%") + geom_label_repel(aes(fill = ifelse(high/low > 2 & high > 3 & gene %in% micro_markers$marker.gene,'red','green'), label = ifelse(high/low > 2 & high > 3,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
p_colored

p_zoomed <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("all top and bottom 20%") + ylim(3, 7.3) + xlim(0, 3.2) + geom_label_repel(aes(fill = ifelse(high/low > 2 & high > 3 & gene %in% micro_markers$marker.gene,'red','green'), label = ifelse(high/low > 2 & high > 3,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
p_zoomed

# plots + Olah's markers
pm1 <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("microglia_subtype1") + geom_label_repel(aes( label = ifelse(high/low > 1.5 & high > 2 & gene %in% m1$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
pm1
pm1_color <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point(aes(x = low, y = high, colour = high/low > 2 & gene %in% m1$gene)) + ggtitle("microglia subtype 1") + scale_colour_manual(name = 'microglia type 1 markers', values = setNames(c('grey','red'),c(F, T)))
pm1_color

pm2 <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("microglia_subtype2") + geom_label_repel(aes( label = ifelse(high/low > 1.5 & high > 2 & gene %in% m2$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
pm2_color <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point(aes(x = low, y = high, colour = high/low > 2 & gene %in% m2$gene)) + ggtitle("microglia subtype 2") + scale_colour_manual(name = 'microglia type 2 markers', values = setNames(c('grey','red'),c(F, T)))

pm3 <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("microglia_subtype3") + geom_label_repel(aes( label = ifelse(high/low > 1.5 & high > 2 & gene %in% m3$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
pm3_color <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point(aes(x = low, y = high, colour = high/low > 2 & gene %in% m3$gene)) + ggtitle("microglia subtype 3") + scale_colour_manual(name = 'microglia type 3 markers', values = setNames(c('grey','red'),c(F, T)))

pm4 <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("microglia_subtype4") + geom_label_repel(aes( label = ifelse(high/low > 1.5 & high > 2 & gene %in% m4$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
pm4_color <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point(aes(x = low, y = high, colour = high/low > 2 & gene %in% m4$gene)) + ggtitle("microglia subtype 4") + scale_colour_manual(name = 'microglia type 4 markers', values = setNames(c('grey','red'),c(F, T)))

pm5 <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("microglia_subtype5") + geom_label_repel(aes( label = ifelse(high/low > 1.5 & high > 2 & gene %in% m5$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
pm5_color <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point(aes(x = low, y = high, colour = high/low > 2 & gene %in% m5$gene)) + ggtitle("microglia subtype 5") + scale_colour_manual(name = 'microglia type 5 markers', values = setNames(c('grey','red'),c(F, T)))

pm6 <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("microglia_subtype6") + geom_label_repel(aes( label = ifelse(high/low > 1.5 & high > 2 & gene %in% m6$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
pm6_color <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point(aes(x = low, y = high, colour = high/low > 2 & gene %in% m6$gene)) + ggtitle("microglia subtype 6") + scale_colour_manual(name = 'microglia type 6 markers', values = setNames(c('grey','red'),c(F, T)))

pm7 <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("microglia_subtype7") + geom_label_repel(aes( label = ifelse(high/low > 1.5 & high > 2 & gene %in% m7$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
pm7_color <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point(aes(x = low, y = high, colour = high/low > 2 & gene %in% m7$gene)) + ggtitle("microglia subtype 7") + scale_colour_manual(name = 'microglia type 7 markers', values = setNames(c('grey','red'),c(F, T)))

pm8 <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("microglia_subtype8") + geom_label_repel(aes( label = ifelse(high/low > 1.5 & high > 2 & gene %in% m8$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
pm8_color <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point(aes(x = low, y = high, colour = high/low > 2 & gene %in% m8$gene)) + ggtitle("microglia subtype 8") + scale_colour_manual(name = 'microglia type 8 markers', values = setNames(c('grey','red'),c(F, T)))

pm9 <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("microglia_subtype9") + geom_label_repel(aes( label = ifelse(high/low > 1.5 & high > 2 & gene %in% m9$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
pm9_color <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point(aes(x = low, y = high, colour = high/low > 2 & gene %in% m9$gene)) + ggtitle("microglia subtype 9") + scale_colour_manual(name = 'microglia type 9 markers', values = setNames(c('grey','red'),c(F, T)))

# use FindMarkers to find DE genes
DefaultAssay(microglia.combined) <- "RNA"
Idents(microglia.combined) <- "microglianess" # assay should be RNA
hi_vs_lo <- FindMarkers(microglia.combined, ident.1 = 'high', ident.2 = 'low', verbose = FALSE)

setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
write.xlsx(hi_vs_lo, file= "m_de.xlsx", sheetName = "m_de", append = FALSE)

# visualize DE using feature plot
FeaturePlot(microglia.combined, features = c("C1qc", "C1qa", "Aif1"), split.by = "microglianess", max.cutoff = 3, 
            cols = c("grey", "red"))

# visualize DE using violin plot
DefaultAssay(microglia.combined) <- "RNA"
microglia.combined$celltype.microglianess <- paste(Idents(microglia.combined), microglia.combined$microglianess, sep = "_")
microglia.combined$celltype <- Idents(microglia.combined)
Idents(microglia.combined) <- "celltype.microglianess"
plots <- VlnPlot(microglia.combined, features = c("C1qc", "C1qa", "Aif1"), split.by = "microglianess", group.by = "celltype", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# UMAP projection + ttype
Idents(microglia.combined) <- "ttype"
microglia.combined$subtype.microglianess <- paste(Idents(microglia.combined), microglia.combined$microglianess, sep = "_")
microglia.combined$subtype <- Idents(microglia.combined)
Idents(microglia.combined) <- "subtype.microglianess"
DimPlot(microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)

# top vs bottom heatmap
DefaultAssay(microglia.combined) <- "integrated"
Idents(microglia.combined) <- "microglianess" # assay should be integrated here
mouse.markers <- FindAllMarkers(microglia.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- mouse.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(microglia.combined, features = top10$gene, group.bar = TRUE, angle = 0) + scale_fill_gradientn(colors = c("lightsteelblue3", "white", "salmon2")) + theme(text = element_text(size = 12))


# 4 HYPERGEOMETRIC TEST ----
library(dplyr)
library(stats)
library(xlsx)

# find each term of the hypergeometric test: q, n, m, k
# perform hypergeometric test on each of the 9 up_types

## OPEN UP_TYPES ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/olah_suppl_data')
m1 <- read.xlsx('m1_de.xlsx', sheetName = 'm1_de')
m2 <- read.xlsx('m2_de.xlsx', sheetName = 'm2_de')
m3 <- read.xlsx('m3_de.xlsx', sheetName = 'm3_de')
m4 <- read.xlsx('m4_de.xlsx', sheetName = 'm4_de')
m5 <- read.xlsx('m5_de.xlsx', sheetName = 'm5_de')
m6 <- read.xlsx('m6_de.xlsx', sheetName = 'm6_de')
m7 <- read.xlsx('m7_de.xlsx', sheetName = 'm7_de')
m8 <- read.xlsx('m8_de.xlsx', sheetName = 'm8_de')
m9 <- read.xlsx('m9_de.xlsx', sheetName = 'm9_de')

m1$gene <- str_to_title(m1$gene)
m2$gene <- str_to_title(m2$gene)
m3$gene <- str_to_title(m3$gene)
m4$gene <- str_to_title(m4$gene)
m5$gene <- str_to_title(m5$gene)
m6$gene <- str_to_title(m6$gene)
m7$gene <- str_to_title(m7$gene)
m8$gene <- str_to_title(m8$gene)
m9$gene <- str_to_title(m9$gene)

setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/olah_suppl_data')
write.xlsx(m1, file = "m_m1.xlsx", sheetName = "m_m1", append = FALSE)
write.xlsx(m2, file = "m_m2.xlsx", sheetName = "m_m2", append = FALSE)
write.xlsx(m3, file = "m_m3.xlsx", sheetName = "m_m3", append = FALSE)
write.xlsx(m4, file = "m_m4.xlsx", sheetName = "m_m4", append = FALSE)
write.xlsx(m5, file = "m_m5.xlsx", sheetName = "m_m5", append = FALSE)
write.xlsx(m6, file = "m_m6.xlsx", sheetName = "m_m6", append = FALSE)
write.xlsx(m7, file = "m_m7.xlsx", sheetName = "m_m7", append = FALSE)
write.xlsx(m8, file = "m_m8.xlsx", sheetName = "m_m8", append = FALSE)
write.xlsx(m9, file = "m_m9.xlsx", sheetName = "m_m9", append = FALSE)


## OPEN MOUSE DE AND ADJUST ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
mouse_DE <- read.xlsx('m_de.xlsx', sheetName = 'm_de')
mouse_DE$fdr <- p.adjust(mouse_DE$p_val, "fdr")
mouse_DE <- filter(mouse_DE, fdr < 0.05)
mouse_DE <- filter(mouse_DE, avg_log2FC > 0) # to get positive markers only
mouse_DE <- dplyr::rename(mouse_DE, 'gene' = 'NA.')
write.xlsx(mouse_DE, file = "m_de_FDR.xlsx", sheetName = "m_de_FDR", append = FALSE)


## FIND q: INNER_JOIN EACH UP_TYPE AND DE BASED ON GENE NAME ----
q1 <- dplyr::inner_join(m1, mouse_DE, by = 'gene', copy = FALSE) # 24
q2 <- dplyr::inner_join(m2, mouse_DE, by = 'gene', copy = FALSE) # 152
q3 <- dplyr::inner_join(m3, mouse_DE, by = 'gene', copy = FALSE) # 78
q4 <- dplyr::inner_join(m4, mouse_DE, by = 'gene', copy = FALSE) # 113
q5 <- dplyr::inner_join(m5, mouse_DE, by = 'gene', copy = FALSE) # 254
q6 <- dplyr::inner_join(m6, mouse_DE, by = 'gene', copy = FALSE) # 208
q7 <- dplyr::inner_join(m7, mouse_DE, by = 'gene', copy = FALSE) # 122
q8 <- dplyr::inner_join(m8, mouse_DE, by = 'gene', copy = FALSE) # 319
q9 <- dplyr::inner_join(m9, mouse_DE, by = 'gene', copy = FALSE) # 385


## FIND n: (OVERLAP BETWEEN OLAH GENES AND OUR HUMAN PATCH SEQ GENES) - m ----
# from Olah et al. 2020, Supplementary Data 3
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/olah_suppl_data")
external_genes = read.xlsx('41467_2020_19737_MOESM5_ESM.xls', sheetName = 'cluster_mean_tpm_20190425') # 25914 genes
external_genes$NA. <- str_to_title(external_genes$NA.)
# from mouse patch seq:
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/mouse_patch_seq/data')
our_genes = read.csv('20200513_Mouse_PatchSeq_Release_cpm.v2.csv') %>% as.data.frame()
our_genes = unlist(our_genes[1])
# find overlapping genes
overlapping_genes = intersect(our_genes, external_genes$NA.)
# get n for hypergeometric testing 
n = length(overlapping_genes) # n = 15942


## PERFORM HYPERGEOMETRIC TEST ----
# 1
phyper(24, 57, 15942-57, 4276, lower.tail = FALSE, log.p = FALSE) # 0.004050799
# 2
phyper(152, 514, 15942-514, 4276, lower.tail = FALSE, log.p = FALSE) # 0.07039337
# 3
phyper(78, 208, 15942-208, 4276, lower.tail = FALSE, log.p = FALSE) # 0.0002654318
# 4
phyper(113, 454, 15942-454, 4276, lower.tail = FALSE, log.p = FALSE) # 0.8125971
# 5
phyper(254, 575, 15942-575, 4276, lower.tail = FALSE, log.p = FALSE) # 2.595097e-20
# 6
phyper(208, 618, 15942-618, 4276, lower.tail = FALSE, log.p = FALSE) # 5.359315e-05
# 7
phyper(122, 318, 15942-318, 4276, lower.tail = FALSE, log.p = FALSE) # 2.181347e-06
# 8
phyper(319, 444, 15942-444, 4276, lower.tail = FALSE, log.p = FALSE) # 1.972736e-90
# 9
phyper(385, 1312, 15942-1312, 4276, lower.tail = FALSE, log.p = FALSE) # 0.01499625


# 5 LINEAR REGRESSION ----
library(dplyr)
library(ggplot2)
library(xlsx)
library(stats)
library(tidyverse)
library(fastDummies)


## LOAD DATA ----
# open metadata with ttype
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
meta <- read.csv("m_metadata_ttype.csv")

# open microglianess
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
micro <- read.xlsx('m_microglia_ranked.xlsx', sheetName = 'm_microglia_ranked')

data <- inner_join(micro, meta, by = "transcriptomics_sample_id")

# remove subtypes with very few cells
data <- data %>% filter(subclass != 'L2/3') %>% filter(subclass != 'Meis2')

data$subclass <- as.factor(data$subclass)
data$donor_id <- as.factor(data$donor_id)


## MAKE PLOTS ----
p1 <- ggplot(
  data,
  aes(x = subclass, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('subclass vs. microglianess') +
  ylab('microglianess') + xlab ('subclass') +
  theme_bw() + geom_boxplot()

p2 <- ggplot(
  data,
  aes(x = cell_soma_normalized_depth, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('cell soma normalized depth vs. microglianess') + ylim(0, 50) + 
  ylab('microglianess') + xlab ('cell soma normalized depth') +
  theme_bw()

p3 <- ggplot(
  data,
  aes(x = donor_id, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('donor ID vs. microglianess') + ylim(0, 50) + 
  ylab('microglianess') + xlab ('donor ID') +
  theme_bw() + geom_boxplot()
'plot_genes[i], " MGP"'

## REGRESSION ----
microglianess <- data$Micro.PVM
subclass <- data$subclass
cell_soma_norm_depth <- data$cell_soma_normalized_depth
donor_id <- data$donor_id

# subclasses only 0.02847
m1 <- lm(microglianess ~ subclass)
summary(m1)

# cell soma normalized depth
m2 <- lm(microglianess ~ cell_soma_norm_depth)
summary(m2)

m3 <- lm(microglianess ~ donor_id)
summary(m3)

ggsave(
  "cell_soma_mouse.jpeg",
  plot = last_plot(),
  device = jpeg,
  path = "/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final",
  width = 15,
  height = 10,
  units = "cm",
)
