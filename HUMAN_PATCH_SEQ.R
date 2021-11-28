## HUMAN PATCH-SEQ ##
# 1 SEURAT TUTORIAL + DIMENSIONALITY REDUCTION ----
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(dplyr)
library(Seurat)
library(patchwork)


## SET UP THE SEURAT OBJECT ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data')
human_data <- read.csv("20200512_Human_PatchSeq_Release_count.csv")
rownames(human_data) <- human_data$X
human_seurat <- CreateSeuratObject(counts = human_data, project = "human", min.cells = 3, min.features = 200)
human_seurat <- SetIdent(human_seurat, value = "human")


## STANDARD PRE-PROCESSING WORKFLOW ----
VlnPlot(human_seurat, features = c("nFeature_RNA", "nCount_RNA"), cols = "lightsteelblue3", ncol = 3)
human_seurat <- subset(human_seurat, subset = nFeature_RNA > 3000 & nFeature_RNA < 17000)


## NORMALIZING THE DATA ----
human_seurat <- NormalizeData(human_seurat, normalization.method = "LogNormalize", scale.factor = 1000000)


## IDENTIFICATION OF HIGHLY VARIABLE FEATURES (FEATURE SELECTION) ----
human_seurat <- FindVariableFeatures(human_seurat, selection.method = "vst", nfeatures = 5000)
top15 <- head(VariableFeatures(human_seurat), 15)
plot1 <- VariableFeaturePlot(human_seurat) + theme(legend.position="none")
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE) + theme(legend.position="none")
plot1 + plot2


## SCALING THE DATA ----
all.genes <- rownames(human_seurat)
human_seurat <- ScaleData(human_seurat, features = all.genes)


## PERFORM LINEAR DIMENSIONAL REDUCTION ----
human_seurat <- RunPCA(human_seurat, features = VariableFeatures(object = human_seurat))
DimHeatmap(human_seurat, dims = 1:10, cells = 278, balanced = TRUE)


## DETERMINE THE ‘DIMENSIONALITY’ OF THE DATASET ----
human_seurat <- JackStraw(human_seurat, num.replicate = 100)
human_seurat <- ScoreJackStraw(human_seurat, dims = 1:20)
JackStrawPlot(human_seurat, dims = 1:20)
ElbowPlot(human_seurat)


## CLUSTER THE CELLS ----
human_seurat <- FindNeighbors(human_seurat, dims = 1:20)
human_seurat <- FindClusters(human_seurat, resolution = 0.5)


## RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/tSNE) ----
reticulate::py_install(packages = 'umap-learn')
human_seurat <- RunUMAP(human_seurat, dims = 1:20)

h_seurat_plot <- DimPlot(human_seurat, reduction = "umap", label = FALSE)

# to save
saveRDS(human_seurat, file = "/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/patch_seq_microglia/h_seurat.rds")

# to read if needed
human_seurat <- readRDS('h_seurat.rds')


## FINDING DIFFERENTIALLY EXPRESSED FEATURES (CLUSTER BIOMARKERS) ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')

human.markers <- FindAllMarkers(human_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
human.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.xlsx(human.markers, file = "h_all_markers.xlsx",
           sheetName = "h_all_markers", append = FALSE)

cluster0.markers <- FindMarkers(human_seurat, ident.1 = 0, min.pct = 0.25)
cluster0.markers.genes <- rownames_to_column(cluster0.markers)
head(cluster0.markers.genes[1], n = 100)
write.xlsx(cluster0.markers, file = "h_cluster0_markers.xlsx",
           sheetName = "h_cluster0_markers", append = FALSE)

cluster1.markers <- FindMarkers(human_seurat, ident.1 = 1, min.pct = 0.25)
cluster1.markers.genes <- rownames_to_column(cluster1.markers)
head(cluster1.markers.genes[1], n = 100)
write.xlsx(cluster1.markers, file = "h_cluster1_markers.xlsx",
           sheetName = "h_cluster1_markers", append = FALSE)

cluster2.markers <- FindMarkers(human_seurat, ident.1 = 2, min.pct = 0.25)
cluster2.markers.genes <- rownames_to_column(cluster2.markers)
head(cluster2.markers.genes[1], n = 100)
write.xlsx(cluster2.markers, file = "h_cluster2_markers.xlsx",
           sheetName = "h_cluster2_markers", append = FALSE)

cluster3.markers <- FindMarkers(human_seurat, ident.1 = 3, min.pct = 0.25)
cluster3.markers.genes <- rownames_to_column(cluster3.markers)
head(cluster3.markers.genes[1], n = 100)
write.xlsx(cluster3.markers, file = "h_cluster3_markers.xlsx",
           sheetName = "h_cluster3_markers", append = FALSE)

# visualization
human_violin <- VlnPlot(human_seurat, features = c("CD74", "TSPO", "C1QC", "AIF1"), ncol = 2)
FeaturePlot(human_seurat, features = c("C1QC", "C1QA", "AIF1","CD74", "SPI1", "TSPO"))

h_seurat_plot + human_violin

top10 <- human.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(human_seurat, features = top10$gene, group.bar = TRUE, angle = 0) + scale_fill_gradientn(colors = c("cadetblue", "black", "salmon2")) + theme(text = element_text(size = 3))


## SAVE MICROGLIA CLUSTER FOR LATER COMPARISON ----
micro_cluster <- subset(human_seurat, idents = "2")
head(micro_cluster)
micro_cluster <- as.data.frame(as.matrix(GetAssayData(object = micro_cluster, slot = "scale.data")))

# reshape row <-> column
micro_cluster_switched <- micro_cluster %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'X')

micro_cluster_switched = micro_cluster_switched %>% distinct(X, .keep_all = T)
X = micro_cluster_switched$X
rownames(micro_cluster_switched) = rownames(micro_cluster_switched) %>% make.names(unique=T)
micro_switched = micro_cluster_switched[-1] %>% t() %>% as.data.frame()
colnames(micro_switched) = X
micro_switched <- tibble::rownames_to_column(micro_switched, var = "transcriptomics_sample_id")

write.csv(micro_switched, file = "h_seurat_micro_switched.csv")


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
human_markers = read.csv('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/kcni_summer/marker_genes/ACC_results.csv')
cell_types = human_markers$cluster %>% unique()
knitr::kable(head(human_markers))


## PROCESS MARKERS DATA FRAME FOR MGP FUNCTION ----
marker_list <- lapply(cell_types, function(cell_type){
  return(human_markers %>% filter(cluster == cell_type) %>% pull(gene) %>% unlist())
})
names(marker_list) <- cell_types
print(cell_types)


## LIST NAMES OF MARKER GENES OF THE 'MICROGLIA' CELL TYPE ----
print(marker_list$'Micro_C1QC')
micro_markers_save <- as.data.frame(marker_list$'Micro_C1QC')
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
write.xlsx(micro_markers_save, file = "h_microglia_markers.xlsx",
           sheetName = "h_microglia_markers", append = FALSE)

## READ AND PROCESS SCRNA SEQ PORTION OF HUMAN PATCH SEQ ----
# Read in the metadata and gene expression count matrix
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data")
human_meta <- read.csv('20200625_patchseq_metadata_human.csv') %>%
  dplyr::select(transcriptomics_sample_id, 
                ttype = "corresponding_AIT2.3.1_alias")

human_expr <- read.csv("20200512_Human_PatchSeq_Release_count.csv") %>%
  tibble::column_to_rownames("X") %>% as.data.frame()
knitr::kable(human_expr[1:5, 1:5])


## NORMALIZE THE GENE EXPRESSION DATASET INTO CPM ----
human_cpm <- edgeR::cpm(human_expr, log = F, prior.count = 0)
knitr::kable(human_cpm[1:5, 1:5])

# or
# setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data")
# human_cpm <- read.csv("20200512_Human_PatchSeq_Release_cpm.csv") %>% tibble::column_to_rownames("X") %>% as.data.frame()


## PROCESS EXPRESSION MATRIX PRIOR TO MGP ANALYSIS ----
use_gene_list <- human_cpm %>% rownames() %>% unlist %>% unique()
gene_mat <- human_cpm[use_gene_list, ]


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

human_meta$transcriptomics_sample_id = sub("-",".", human_meta$transcriptomics_sample_id)
human_meta$transcriptomics_sample_id = sub("-",".", human_meta$transcriptomics_sample_id)
human_meta$transcriptomics_sample_id = sub("-",".", human_meta$transcriptomics_sample_id)

# merge gene expression and metadata frames
gene_mat_comb <- inner_join(human_meta, gene_mat_trans, by = 'transcriptomics_sample_id')


## RUN CELL TYPE PROPORTION ESTIMATION ----
# run MGP analysis
estimations <-  mgpEstimate(
  exprData=gene_mat,
  genes=marker_list,
  geneColName='X',
  outlierSampleRemove=F, # outlier samples removed. This is done using boxplot stats.
  geneTransform = NULL, # this is the default option for geneTransform
  groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus=FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority=TRUE)


## MERGE RELATIVE CELL TYPE PROPORTIONS PER CELL TYPE WITH SAMPLE METADATA ----
mgp_estimates <- as.data.frame(estimations$estimates) %>%
  tibble::rownames_to_column(var = 'transcriptomics_sample_id')

mgp_df <- inner_join(human_meta, mgp_estimates, by = "transcriptomics_sample_id")
knitr::kable(mgp_df[1:5, 1:12])


## MAKE FILES CONTAINING CELLS AND THEIR 'MICROGLIANESS' ----
# ranked
microglia_ranked <- mgp_df %>% dplyr::select(c(transcriptomics_sample_id, Micro_C1QC)) %>% arrange(desc(Micro_C1QC))
ttype <- human_meta %>% dplyr::select(c(transcriptomics_sample_id, ttype))
microglia_ttype <- inner_join(microglia_ranked, ttype, by = "transcriptomics_sample_id")
top25 <- slice_max(microglia_ttype, order_by = Micro_C1QC, prop = 0.25) # microglianess cutoff is 2.191213 inclusive
bottom25 <- slice_min(microglia_ttype, order_by = Micro_C1QC, prop = 0.25) # microglianess cutoff is -4.765124 inclusive

microglia_ttype$group <- "low"
microglia_ttype$group[microglia_ttype$Micro_C1QC > -4.765124 & microglia_ttype$Micro_C1QC < 2.191213] <- "med"
microglia_ttype$group[microglia_ttype$Micro_C1QC >= 2.191213] <- "high"

setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
write.xlsx(microglia_ttype, file = "h_microglia_ttype.xlsx",
           sheetName = "h_microglia_ttype", append = FALSE)

# unranked
microglia_unranked <- mgp_df %>% dplyr::select(c(transcriptomics_sample_id, Micro_C1QC))
order <- as.data.frame(microglia_unranked$transcriptomics_sample_id)
order <- rename(order, 'transcriptomics_sample_id' = 'microglia_unranked$transcriptomics_sample_id')
microglia_ttype_unranked <- inner_join(order, microglia_ttype, by = "transcriptomics_sample_id")
write.xlsx(microglia_ttype_unranked, file = "h_microglia_ttype_unranked.xlsx",
           sheetName = "h_microglia_ttype_unranked", append = FALSE)

# now, take column of mouse_cpm (store as new data) and inner join it to the first column of microglia_ranked to get a ranked data frame
microglia_ranked_full <- inner_join(microglia_ranked, gene_mat_trans, by = "transcriptomics_sample_id")
write.xlsx(microglia_ranked_full[, 1:2], file = "h_microglia_rank.xlsx",
           sheetName = "h_microglia_rank", append = FALSE)


## LABONTE CLUSTER vs SEURAT CLUSTER ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
micro_switched <- read.csv('h_seurat_micro_switched.csv') # has 73 cells
micro_switched <- as.data.frame(micro_switched[2])
micro_overlap <- inner_join(micro_switched, microglia_ttype_unranked, by = "transcriptomics_sample_id", copy = FALSE) # 20 in med, 53 in high (there are 69 high cells in total)

# rm(list = ls())


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
# load mouse dataset and reshape to fit h_microglia_ttype_unranked.xlsx
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data')
human_data <- read.csv("20200512_Human_PatchSeq_Release_count.csv")

human_data <- human_data %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'transcript')

human_data = human_data %>% distinct(transcript, .keep_all = T)
transcript = human_data$transcript
rownames(human_data) = rownames(human_data) %>% make.names(unique=T)
human_data_trans = human_data[-1] %>% t() %>% as.data.frame()
colnames(human_data_trans) = transcript
human_data_trans = human_data_trans %>% row_to_names(row_number = 1)
human_data_trans # now we have the count matrix rotated 90 degrees

# load microglia_ttype_unranked.xlsx
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
micro <- read.xlsx('h_microglia_ttype_unranked.xlsx', sheetName = 'h_microglia_ttype_unranked')
micro <- dplyr::select(micro, -c(NA.)) # removed first column of numbers

# create new dataframes: top_25 and bottom_25 (and save as xlsx)
top_25 <- filter(micro, group == 'high')
bottom_25 <- filter(micro, group == 'low')

write.xlsx(top_25, file = "h_top_25_microglia.xlsx",
           sheetName = "h_top_25_microglia", append = FALSE)
write.xlsx(bottom_25, file = "h_bottom_20_microglia.xlsx",
           sheetName = "h_bottom_25_microglia", append = FALSE)

# inner_join the transcriptomics_sample_id of top_20 and bottom_20 with mouse_data_trans, respectively (save as csv)
top_25 <- as.data.frame(top_25[, 1])
top_25 <- rename(top_25, 'transcriptomics_sample_id' = 'top_25[, 1]')

bottom_25 <- as.data.frame(bottom_25[, 1])
bottom_25 <- rename(bottom_25, 'transcriptomics_sample_id' = 'bottom_25[, 1]')

human_data_trans <- human_data_trans %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'transcriptomics_sample_id')

micro_top25 <- inner_join(human_data_trans, top_25, by = 'transcriptomics_sample_id', copy = T)
micro_bottom25 <- inner_join(human_data_trans, bottom_25, by = 'transcriptomics_sample_id', copy = T)

write.csv(micro_top25, file = "h_top_25_micro_count_trans.csv")
write.csv(micro_bottom25, file = "h_bottom_25_micro_count_trans.csv")

# reshape micro_top25 and micro_bottom25 and save them as xlsx files (these are the files to work with for the integration)
micro_top25 <- micro_top25 %>% as.data.frame()
micro_top25 = micro_top25 %>% distinct(transcriptomics_sample_id, .keep_all = T)
transcriptomics_sample_id = micro_top25$transcriptomics_sample_id
rownames(micro_top25) = rownames(micro_top25) %>% make.names(unique=T)
human_top = micro_top25[-1] %>% t() %>% as.data.frame()
colnames(human_top) = transcriptomics_sample_id

micro_bottom25 <- micro_bottom25 %>% as.data.frame()
micro_bottom25 = micro_bottom25 %>% distinct(transcriptomics_sample_id, .keep_all = T)
transcriptomics_sample_id = micro_bottom25$transcriptomics_sample_id
rownames(micro_bottom25) = rownames(micro_bottom25) %>% make.names(unique=T)
human_bottom = micro_bottom25[-1] %>% t() %>% as.data.frame()
colnames(human_bottom) = transcriptomics_sample_id

write.csv(human_top, file = "h_micro_top25_count.csv")
write.csv(human_bottom, file = "h_micro_bottom25_count.csv")


## JOIN METADATA + TTYPE WITH TOP20 AND BOTTOM20 ----
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/kcni_summer")
human_top <- read.csv("h_top_25_micro_count_trans.csv")
human_bottom <- read.csv("h_bottom_25_micro_count_trans.csv")
human_top <- as.data.frame(human_top[, 2])
human_bottom <- as.data.frame(human_bottom[, 2])
human_top <- rename(human_top, "transcriptomics_sample_id" = "human_top[, 2]")
human_bottom <- rename(human_bottom, "transcriptomics_sample_id" = "human_bottom[, 2]")

# open metadata and add ttype
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data")
human_meta <- read.csv("20200625_patchseq_metadata_human.csv")
human_meta$s1 <- gsub( "^\\S+\\s+", "", human_meta$corresponding_AIT2.3.1_alias)
human_meta$s2 <- gsub( "^\\S+\\s+", "", human_meta$s1)
human_meta$subclass <- as.factor(gsub( "^\\S+\\s+", "", human_meta$s2))
human_meta <- human_meta %>% dplyr::select(-c('s1', 's2'))

setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
write.csv(human_meta, file = "h_metadata_ttype.csv")

# open metadata with ttype
meta <- human_meta %>% dplyr::select(transcriptomics_sample_id, subclass)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)

# inner_join top25 and bottom25 with meta; save as csv
top25sub <- inner_join(human_top, meta, by = 'transcriptomics_sample_id', copy = T)
bottom25sub <- inner_join(human_bottom, meta, by = 'transcriptomics_sample_id', copy = T)

write.csv(top25sub, file = "h_metadata_ttype_top25.csv")
write.csv(bottom25sub, file = "h_metadata_ttype_bottom25.csv")

# remove everything to clear up space: rm(list = ls())


## SET UP HIGH MICROGLIANESS SEURAT OBJECT AND JOIN WITH METADATA ----
# open h_micro_top25_count.csv 
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
top <- read.csv('h_micro_top25_count.csv')
rownames(top) <- top$X

top_seurat <- CreateSeuratObject(counts = top, project = "high", min.cells = 3, min.features = 200)

# add ttype/subclass
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final")
top25_meta <- read.csv("h_metadata_ttype_top25.csv")
subclass_f <- as.factor(top25_meta$subclass)
top_seurat <- AddMetaData(top_seurat, col.name = "ttype", metadata = subclass_f)
# don't do Idents(top_seurat) <- "ttype", because it'll change "high" to ttypes

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
JackStrawPlot(top_seurat, dims = 1:20) # pick 17 dims... because dim17 is on the line
ElbowPlot(top_seurat)

# cluster the cells
top_seurat <- FindNeighbors(top_seurat, dims = 1:17)
top_seurat <- FindClusters(top_seurat, resolution = 0.5)
head(Idents(top_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
top_seurat <- RunUMAP(top_seurat, dims = 1:17)

DimPlot (top_seurat, reduction = "umap", label = TRUE) # there's only 1 cluster
DimPlot(top_seurat, reduction = "umap", pt.size = 1, group.by = "microglianess", label = TRUE)


## SET UP LOW MICROGLIANESS SEURAT OBJECT ----
# open h_micro_bottom25_count.csv 
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
bottom <- read.csv('h_micro_bottom25_count.csv')
rownames(bottom) <- bottom$X

bottom_seurat <- CreateSeuratObject(counts = bottom, project = "low", min.cells = 3, min.features = 200)

# add ttype/subclass
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final")
bottom25_meta <- read.csv("h_metadata_ttype_bottom25.csv")
subclass_f <- as.factor(bottom25_meta$subclass)
bottom_seurat <- AddMetaData(bottom_seurat, col.name = "ttype", metadata = subclass_f)
# don't do Idents(top_seurat) <- "ttype", because it'll change "low" to ttypes

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
JackStrawPlot(bottom_seurat, dims = 1:20) # take 16 dims
ElbowPlot(bottom_seurat)

# cluster the cells
bottom_seurat <- FindNeighbors(bottom_seurat, dims = 1:16)
bottom_seurat <- FindClusters(bottom_seurat, resolution = 0.5)
head(Idents(bottom_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
bottom_seurat <- RunUMAP(bottom_seurat, dims = 1:16)

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
microglia.combined <- IntegrateData(anchorset = microglia.anchors, dims = 1:20, k.weight = 65)


## PERFORM AN INTEGRATED ANALYSIS ----
# specify that we will perform downstream analysis on the corrected data note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(microglia.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
microglia.combined <- ScaleData(microglia.combined, verbose = FALSE)
microglia.combined <- RunPCA(microglia.combined, npcs = 20, verbose = FALSE)
microglia.combined <- RunUMAP(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindNeighbors(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindClusters(microglia.combined, resolution = 0.5)

saveRDS(microglia.combined, file = "/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final/h_microglia_combined.rds")

# if want to read it
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
microglia.combined <- readRDS('h_microglia_combined.rds')

# Visualization
p1 <- DimPlot(microglia.combined, reduction = "umap", group.by = "microglianess")
p2 <- DimPlot(microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(microglia.combined, reduction = "umap", split.by = "microglianess")
p1 + p2
p3

# save figure
ggsave(
  "h_microglia_combined_umap.jpeg",
  plot = p1,
  device = jpeg,
  path = "/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final",
  width = 15,
  height = 10,
  units = "cm",
)

# more visualization
DefaultAssay(microglia.combined) <- "RNA"
FeaturePlot(microglia.combined, features = c("C1QC", "C1QA", "AIF1", "SPI1"), min.cutoff = "q9")
markers.to.plot <- c("C1QC", "C1QA", "AIF1", "SPI1", "TSPO", "CD74", "ADAP2", "C3AR1", "CCR5", "CD14")
DotPlot(microglia.combined, features = markers.to.plot, cols = c("mediumturquoise", "coral2"), dot.scale = 8, split.by = "microglianess") + RotatedAxis()


## IDENTIFY DIFERENTIAL EXPRESSED GENES ACROSS CONDITIONS ----
# make high vs. low expression (for graphs)
DefaultAssay(microglia.combined) <- "integrated"
Idents(microglia.combined2) <- "microglianess"
avg.microglia.combined <- as.data.frame(log1p(AverageExpression(microglia.combined2, verbose = FALSE)$RNA))
avg.microglia.combined$gene <- rownames(avg.microglia.combined)

# open external mouse microglia markers (for graphs)
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/kcni_summer")
micro_markers <- read.xlsx("microglia_markers_human.xlsx", sheetName = "microglia_markers_human")
micro_markers <- as.data.frame(micro_markers[, 2])
micro_markers <- dplyr::rename(micro_markers, 'marker.gene' = 'micro_markers[, 2]')

# open Olah's markers (for graphs)
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

options(ggrepel.max.overlaps = Inf)

# plots
p_basic <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("all top and bottom 25%") + geom_text(aes(label=ifelse(high/low > 2 & high > 3, as.character(gene),'')),hjust=0,vjust=0)
p_basic

p_colored <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("all top and bottom 25%") + geom_label_repel(aes(fill = ifelse(high/low > 2 & high > 3 & gene %in% micro_markers$marker.gene,'red','green'), label = ifelse(high/low > 2 & high > 3,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
p_colored

p_zoomed <- ggplot(avg.microglia.combined, aes(low, high)) + geom_point() + ggtitle("all top and bottom 25%") + ylim(3, 7.3) + xlim(0, 3.2) + geom_label_repel(aes(fill = ifelse(high/low > 2 & high > 3 & gene %in% micro_markers$marker.gene,'red','green'), label = ifelse(high/low > 2 & high > 3,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
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
write.xlsx(hi_vs_lo, file= "h_de.xlsx", sheetName = "h_de", append = FALSE)

# visualize DE using feature plot
FeaturePlot(microglia.combined, features = c("HLA-DRA", "C1QA", "AIF1"), split.by = "microglianess", max.cutoff = 3, 
            cols = c("grey", "red"))

# visualize DE using violin plot
DefaultAssay(microglia.combined) <- "RNA"
microglia.combined$celltype.microglianess <- paste(Idents(microglia.combined), microglia.combined$microglianess, sep = "_")
microglia.combined$celltype <- Idents(microglia.combined)
Idents(microglia.combined) <- "celltype.microglianess"
plots <- VlnPlot(microglia.combined, features = c("C1QC", "C1QA", "AIF1"), split.by = "microglianess", group.by = "celltype", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# UMAP projection + ttype
Idents(microglia.combined) <- "ttype"
microglia.combined$subtype.microglianess <- paste(Idents(microglia.combined), microglia.combined$microglianess, sep = "_")
microglia.combined$subtype <- Idents(microglia.combined)
Idents(microglia.combined) <- "subtype.microglianess"
#Idents(microglia.combined) <- microglia.combined$subtype.microglianess
DimPlot(microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)

# top vs bottom heatmap logFC
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


## OPEN MOUSE DE AND ADJUST ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
human_DE <- read.xlsx('h_de.xlsx', sheetName = 'h_de')
human_DE$fdr <- p.adjust(human_DE$p_val, "fdr")
human_DE <- filter(human_DE, fdr < 0.05)
human_DE <- filter(human_DE, avg_log2FC > 0)
human_DE <- dplyr::rename(human_DE, 'gene' = 'NA.')
write.xlsx(human_DE, file = "h_de_FDR.xlsx", sheetName = "h_de_FDR", append = FALSE)


## FIND q: INNER_JOIN EACH UP_TYPE AND DE BASED ON GENE NAME ----
q1 <- dplyr::inner_join(m1, human_DE, by = 'gene', copy = FALSE) # 39
q2 <- dplyr::inner_join(m2, human_DE, by = 'gene', copy = FALSE) # 171
q3 <- dplyr::inner_join(m3, human_DE, by = 'gene', copy = FALSE) # 78
q4 <- dplyr::inner_join(m4, human_DE, by = 'gene', copy = FALSE) # 132
q5 <- dplyr::inner_join(m5, human_DE, by = 'gene', copy = FALSE) # 294
q6 <- dplyr::inner_join(m6, human_DE, by = 'gene', copy = FALSE) # 225
q7 <- dplyr::inner_join(m7, human_DE, by = 'gene', copy = FALSE) # 166
q8 <- dplyr::inner_join(m8, human_DE, by = 'gene', copy = FALSE) # 360
q9 <- dplyr::inner_join(m9, human_DE, by = 'gene', copy = FALSE) # 327


## FIND n: (OVERLAP BETWEEN OLAH GENES AND OUR HUMAN PATCH SEQ GENES) - m ----
# from Olah et al. 2020, Supplementary Data 3
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/olah_suppl_data")
external_genes = read.xlsx('41467_2020_19737_MOESM5_ESM.xls', sheetName = 'cluster_mean_tpm_20190425')
# from Berg et al. 2020: https://github.com/keon-arbabi/berg_patchseq/blob/main/data/20200512_Human_PatchSeq_Release_count.csv
our_genes = fread("https://raw.githubusercontent.com/keon-arbabi/berg_patchseq/main/data/20200512_Human_PatchSeq_Release_cpm.csv?token=AUB7ZHZIO5RYZNN74WJSVILA2M3UO") %>% as.data.frame()
our_genes = unlist(our_genes[1])
# find overlapping genes
overlapping_genes = intersect(our_genes, external_genes$NA.)
# get n for hypergeometric testing 
n = length(overlapping_genes)


## PERFORM HYPERGEOMETRIC TEST ----
# 1
phyper(39, 57, 21698-57, 3495, lower.tail = FALSE, log.p = FALSE) # 1.248064e-19
# 2
phyper(171, 500, 21698-514, 3495, lower.tail = FALSE, log.p = FALSE) # 3.325731e-24
# 3
phyper(78, 208, 21698-208, 3495, lower.tail = FALSE, log.p = FALSE) # 1.797232e-14
# 4
phyper(132, 454, 21698-454, 3495, lower.tail = FALSE, log.p = FALSE) # 7.923177e-13
# 5
phyper(295, 500, 21698-575, 3495, lower.tail = FALSE, log.p = FALSE) # 9.605125e-109
# 6
phyper(225, 500, 21698-618, 3495, lower.tail = FALSE, log.p = FALSE) # 7.448115e-54
# 7
phyper(166, 318, 21698-318, 3495, lower.tail = FALSE, log.p = FALSE) # 2.139353e-51
# 8
phyper(360, 444, 21698-444, 3495, lower.tail = FALSE, log.p = FALSE) # 3.444141e-208
# 9
phyper(327, 500, 21698-1312, 3495, lower.tail = FALSE, log.p = FALSE) # 2.562925e-135


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
meta <- read.csv("h_metadata_ttype.csv")
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)

# open microglianess
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
micro <- read.xlsx('h_microglia_rank.xlsx', sheetName = 'h_microglia_rank')

data <- inner_join(micro, meta, by = "transcriptomics_sample_id", copy = F)

## FIX VARIABLES ----
# make age numeric
data$age_years <- as.numeric(gsub( " .*$", "", data$age))

# make sex, subclass, donor id and condition factor (but first make new subtypes according to paper) #
# sex
data$biological_sex <- as.factor(data$biological_sex)

# subclass
data$s1 <- gsub( "^\\S+\\s+", "", data$corresponding_AIT2.3.1_alias)
data$s2 <- gsub( "^\\S+\\s+", "", data$s1)
data$subclass <- as.factor(gsub( "^\\S+\\s+", "", data$s2))

# donor id
data$donor_id <- as.factor(data$donor_id)

# condition
data$medical_conditions <- as.factor(data$medical_conditions)

## MAKE PLOTS ----
p1 <- ggplot(
  data,
  aes(x = age_years, y = Micro_C1QC)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('age vs. microglianess') +
  ylab('microglianess') + xlab ('age (years)') +
  theme_bw()

p2 <- ggplot(
  data,
  aes(x = donor_id, y = Micro_C1QC)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('donor ID vs. microglianess') +
  ylab('microglianess') + xlab ('donor ID') +
  theme_bw() + geom_boxplot()

p3 <- ggplot(
  data,
  aes(x = medical_conditions, y = Micro_C1QC)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('medical condition vs. microglianess') +
  ylab('microglianess') + xlab ('medical condition') +
  theme_bw() + geom_boxplot()

p4 <- ggplot(
  data,
  aes(x = biological_sex, y = Micro_C1QC)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('biological sex vs. microglianess') +
  ylab('microglianess') + xlab ('biological sex') +
  theme_bw() + geom_boxplot()

p5 <- ggplot(
  data,
  aes(x = subclass, y = Micro_C1QC)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('subclass vs. microglianess') +
  ylab('microglianess') + xlab ('subclass') +
  theme_bw() + geom_boxplot()

p6 <- ggplot(
  data,
  aes(x = cell_soma_normalized_depth, y = Micro_C1QC)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('cell soma normalized depth vs. microglianess') +
  ylab('microglianess') + xlab ('cell soma normalized depth') +
  theme_bw()

p1 + p2 + p3 + p4 + p5 + p6


## REGRESSION ----
# make dummy variables for sex, conditions, subclass and donor_id
data <- dummy_cols(data, select_columns = c('biological_sex', 'medical_conditions', 'subclass', 'donor_id'))

microglianess <- data$Micro_C1QC
conditions <- data$medical_conditions_ # epilepsy is 0, tumor is 1
sex <- data$biological_sex_F # male is 0, female is 1
#subclass_CARM1P1 <- data$subclass_CARM1P1
#subclass_COL22A1 <- data$subclass_COL22A1
#subclass_FREM3 <- data$subclass_FREM3
#subclass_GLP2R <- data$subclass_GLP2R
#subclass_LTK <- data$subclass_LTK
subclass <- data$subclass
age <- data$age_years
donor <- data$donor_id
cell_soma_norm_depth <- data$cell_soma_normalized_depth

## MAKE MODELS ----
# all (none are significant, R squared = 0.04786)
m1 <- lm(microglianess ~ depth + conditions + sex + age + subclass_CARM1P1 + subclass_COL22A1 + subclass_FREM3 + subclass_GLP2R + subclass_LTK + donor)
summary(m1)
plot(m1)

# donor id only 0.3498
m2 <- lm(microglianess ~ donor)
summary(m2)

# age only 0.002278
m4 <- lm(microglianess ~ age)
summary(m4)

# sex only 0.001427
m5 <- lm(microglianess ~ sex)
summary(m5)

# conditions only 0.0008469
m6 <- lm(microglianess ~ conditions)
summary(m6)

# subclasses only 0.02847
m7 <- lm(microglianess ~ subclass)
summary(m7)

# cell soma normalized depth only 0.03267
m8 <- lm(microglianess ~ cell_soma_norm_depth)
summary(m8)


# 5 IBA1 PATHOLOGY VS MICROGLIANESS AND AIF1 MRNA ----
library(dplyr)
library(xlsx)
library(data.table)
library(ggplot2)


## JOIN RANKED HUMAN MICROGLIA CELLS (MGP) AND METADATA DONOR ID ----
setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
ranked_microglia <- read.xlsx('h_microglia_rank.xlsx', sheetName = 'h_microglia_rank')

setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data')
meta <- read.csv('20200625_patchseq_metadata_human.csv')
meta2 <- select(meta, donor_id, transcriptomics_sample_id)
meta2$transcriptomics_sample_id = sub("-",".", meta2$transcriptomics_sample_id)
meta2$transcriptomics_sample_id = sub("-",".", meta2$transcriptomics_sample_id)
meta2$transcriptomics_sample_id = sub("-",".", meta2$transcriptomics_sample_id)

tojoin <- inner_join(ranked_microglia, meta2, by = "transcriptomics_sample_id")
tojoin <- tojoin[-1]


## JOIN PATHOLOGY DONOR ID AND IBA1 CORTEX LEVELS ----
# https://github.com/AllenInstitute/patchseq_human_L23/blob/master/data/pathology_scoring.csv
pathology <- fread('https://raw.githubusercontent.com/AllenInstitute/patchseq_human_L23/master/data/pathology_scoring.csv') %>% as.data.frame()
pathology <- select(pathology, donor_id, 'Iba1 Cortex')

setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
write.xlsx(pathology, file = "h_patho.xlsx",
           sheetName = "h_patho", append = FALSE) # then adjust in excel :(

# don't know how to remove things after period '.'

#after editing on excel, continue
pathology <- read.xlsx('h_patho.xlsx', sheetName = 'h_patho')

tojoin2 <- inner_join(tojoin, pathology, by = "donor_id", copy = FALSE)
tojoin2 <- tojoin2[-4]


## JOIN WITH RESHAPED CPM WITH TRANSCRIPTOMICS SAMPLE ID AND IBA1 MRNA EXPRESSION ----
cpm <- read.csv('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/berg_patchseq-main/data20200512_Human_PatchSeq_Release_cpm.csv')
cpm <- cpm %>% as.data.frame()
cpm = cpm %>% distinct(X, .keep_all = T)
X = cpm$X
rownames(cpm) = rownames(cpm) %>% make.names(unique=T)
cpm_trans = cpm[-1] %>% t() %>% as.data.frame()
colnames(cpm_trans) = X

cpm_aif1 <- cpm_trans %>% tibble::rownames_to_column(var = 'transcriptomics_sample_id') %>% select(c('transcriptomics_sample_id', 'AIF1'))

joined <- inner_join(tojoin2, cpm_aif1, by = "transcriptomics_sample_id")

setwd('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/final')
write.xlsx(joined, file = "h_patho_microglianess_aif1.xlsx",
           sheetName = "h_patho_microglianess_aif1", append = FALSE)


## MAKE GRAPHS OF IBA1 VS. MICROGLIANESS AND IBA1 VS. AIF1 ----
iba1_vs_microglianess <- joined %>% 
  ggplot(aes(x = Iba1.Cortex, y = Micro_C1QC, group = 1, color = donor_id)) +  
  geom_smooth(method = "lm", se = F) + geom_point() + 
  ylab('microglianess') + xlab('IBA1 protein level') + ggtitle('IBA1 CORTEX VS. MICROGLIANESS')
iba1_vs_microglianess

iba1_vs_aif1 <- joined %>% 
  ggplot(aes(x = Iba1.Cortex, y = AIF1, group = 1, color = donor_id)) +  
  geom_smooth(method = "lm", se = F) + geom_point() + 
  ylab('AIF1 level') + xlab('IBA1 protein level') + ggtitle('IBA1 CORTEX VS. AIF1 LEVEL')
iba1_vs_aif1


