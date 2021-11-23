#### SEURAT TUTORIAL ####
library(dplyr)
library(Seurat)
library(patchwork)

#### SET UP THE SEURAT OBJECT ####

# open microglia count
setwd('/external/rprshnas01/kcni/yjiang')
scRNA <- read.csv('scRNA_microglia_count.csv')

# make seurat object of microglia cells
scRNA_seurat <- CreateSeuratObject(counts = scRNA, project = "scRNA", min.cells = 3, min.features = 200)

# normalization
scRNA_seurat <- NormalizeData(scRNA_seurat, normalization.method = "LogNormalize", scale.factor = 1000000)

# highly variable features
scRNA_seurat <- FindVariableFeatures(scRNA_seurat, selection.method = "vst", nfeatures = 5000)

# scaling the data
all.genes <- rownames(scRNA_seurat)
scRNA_seurat <- ScaleData(scRNA_seurat, features = all.genes)

# linear dimensional reduction
scRNA_seurat <- RunPCA(scRNA_seurat, features = VariableFeatures(object = scRNA_seurat))
DimHeatmap(scRNA_seurat, dims = 1:6, cells = 300, balanced = TRUE)

# determine dimensionality of dataset
scRNA_seurat <- JackStraw(scRNA_seurat, num.replicate = 100)
scRNA_seurat <- ScoreJackStraw(scRNA_seurat, dims = 1:20)
JackStrawPlot(scRNA_seurat, dims = 1:20) # pick 20 dims
ElbowPlot(scRNA_seurat)

# cluster the cells
scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:20)
scRNA_seurat <- FindClusters(scRNA_seurat, resolution = 0.5)
head(Idents(scRNA_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
scRNA_seurat <- RunUMAP(scRNA_seurat, dims = 1:20)

DimPlot (scRNA_seurat, reduction = "umap", label = TRUE) 

human_hi_vs_lo <- FindMarkers(microglia.combined, ident.1 = "0_high", ident.2 = "0_low", min.pct = 0.25, verbose = FALSE)


