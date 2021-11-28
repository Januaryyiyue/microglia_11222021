# 1 INTEGRATING PATCH AND SCRNA ----
## SETUP ----

library(Seurat)
library(dplyr)
library(xlsx)
library(janitor)
library(ggplot2)
library(cowplot)
library(ggrepel)
theme_set(theme_cowplot())
library(patchwork)

scRNA_seurat <- readRDS('yao_microglia_seurat.RDS')

scRNA_meta <- as.data.frame(scRNA_seurat@meta.data)
scRNA_temp <- as.data.frame(scRNA_seurat@assays$RNA@counts)

# rotate scRNA_temp
scRNA_count <- scRNA_temp %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'transcript')

scRNA_count = scRNA_count %>% distinct(transcript, .keep_all = T)
transcript = scRNA_count$transcript
rownames(scRNA_count) = rownames(scRNA_temp) %>% make.names(unique=T)
sc_count = scRNA_count[-1] %>% t() %>% as.data.frame()


## PATCH SEQ DATA ----
## set up patch-seq high microglia cells seurat object
# open h_micro_top25_count copy.csv (from my laptop)
patch <- read.csv("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/figures/files_too_large_for_git/m_QC_micro_top20_count.csv")
patch <- tibble::column_to_rownames(patch, var = 'X')

patch_seurat <- CreateSeuratObject(counts = patch, project = "patchseq", min.cells = 3, min.features = 200)

# normalization
patch_seurat <- NormalizeData(patch_seurat, normalization.method = "LogNormalize", scale.factor = 1000000)

# highly variable features
patch_seurat <- FindVariableFeatures(patch_seurat, selection.method = "vst", nfeatures = 5000)

# scaling the data
all.genes <- rownames(patch_seurat)
patch_seurat <- ScaleData(patch_seurat, features = all.genes)

# add metadata
patch_seurat@meta.data[, "source"] <- "patchseq"

# linear dimensional reduction
patch_seurat <- RunPCA(patch_seurat, features = VariableFeatures(object = patch_seurat))

# determine dimensionality of dataset
patch_seurat <- JackStraw(patch_seurat, num.replicate = 100)
patch_seurat <- ScoreJackStraw(patch_seurat, dims = 1:20)
JackStrawPlot(patch_seurat, dims = 1:20)
ElbowPlot(patch_seurat)

# cluster the cells
patch_seurat <- FindNeighbors(patch_seurat, dims = 1:20)
patch_seurat <- FindClusters(patch_seurat, resolution = 0.5)
head(Idents(patch_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
patch_seurat <- RunUMAP(patch_seurat, dims = 1:20)

DimPlot (patch_seurat, reduction = "umap", label = TRUE) 
DimPlot(patch_seurat, reduction = "umap", pt.size = 1, group.by = "source", label = TRUE)


## SCRNA SEQ DATA ----
## make seurat object of microglia cells
scRNA_seurat <- CreateSeuratObject(counts = scRNA_count, project = "scRNA", min.cells = 3, min.features = 200)

# normalization
scRNA_seurat <- NormalizeData(scRNA_seurat, normalization.method = "LogNormalize", scale.factor = 1000000)

# highly variable features
scRNA_seurat <- FindVariableFeatures(scRNA_seurat, selection.method = "vst", nfeatures = 5000)

# scaling the data
all.genes <- rownames(scRNA_seurat)
scRNA_seurat <- ScaleData(scRNA_seurat, features = all.genes)

# add metadata
scRNA_seurat@meta.data[, "source"] <- "scRNA"

# linear dimensional reduction
scRNA_seurat <- RunPCA(scRNA_seurat, features = VariableFeatures(object = scRNA_seurat))

# determine dimensionality of dataset
scRNA_seurat <- JackStraw(scRNA_seurat, num.replicate = 100)
scRNA_seurat <- ScoreJackStraw(scRNA_seurat, dims = 1:20)
JackStrawPlot(scRNA_seurat, dims = 1:20) # pick 20 dims... because dim17 is on the line
ElbowPlot(scRNA_seurat)

# cluster the cells
scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:20)
scRNA_seurat <- FindClusters(scRNA_seurat, resolution = 0.5)
head(Idents(scRNA_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
scRNA_seurat <- RunUMAP(scRNA_seurat, dims = 1:20)

DimPlot (scRNA_seurat, reduction = "umap", label = TRUE) 
DimPlot(scRNA_seurat, reduction = "umap", pt.size = 1, group.by = "source", label = TRUE)


## COMBINE AND SPLIT SEURAT OBJECTS ----
patch_scRNA_combined = merge(patch_seurat, y = scRNA_seurat, add.cell.ids = c("patchseq", "scRNA"), project = "source")

# split the dataset into a list of two seurat objects (patchseq and scRNA in terms of source)
microglia.list <- SplitObject(patch_scRNA_combined, split.by = "source")
reference.list <- microglia.list[c("patchseq", "scRNA")]


## PERFORM INTEGRATION ----
microglia.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)
microglia.combined <- IntegrateData(anchorset = microglia.anchors, dims = 1:20)


## PERFORM AN INTEGRATED ANALYSIS ----
DefaultAssay(microglia.combined) <- "integrated"

microglia.combined <- ScaleData(microglia.combined, verbose = FALSE)
microglia.combined <- RunPCA(microglia.combined, npcs = 20, verbose = FALSE)
microglia.combined <- RunUMAP(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindNeighbors(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindClusters(microglia.combined, resolution = 0.5)

saveRDS(microglia.combined, file = "m_patchseq_scRNA.rds")

# if want to read it
microglia.combined <- readRDS('h_patchseq_scRNA.rds')

# visualization
# Visualization
umap <- DimPlot(microglia.combined, reduction = "umap", group.by = "source")

## FIND DIFFERENTIALLY EXPRESSED GENES ----
DefaultAssay(microglia.combined) <- "integrated"
Idents(microglia.combined) <- "source"

patch_vs_scrna <- FindMarkers(microglia.combined, ident.1 = "patchseq", ident.2 = "scRNA", verbose = FALSE)
write.xlsx(patch_vs_scrna, file = "m_patch_scRNA_de.xlsx", sheetName = "m_patch_scRNA_de")


# 2 VOLCANO ----
# https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

de <- patch_vs_scrna
de <- tibble::rownames_to_column(de, var = "gene")

v0 <- ggplot(data = de, aes(x = avg_log2FC, y = p_val)) + geom_point()
v1 <- ggplot(data = de, aes(x = avg_log2FC, y = -log(p_val))) + geom_point() + theme_minimal()

v2 <- v1 + geom_vline(xintercept = c(-2, 2), col = "coral2") + geom_hline(yintercept = -log10(1e-50), col = "coral2")

# add a column of NAs
de$diffexpressed <- "no"
# if avg_log2FC > 2.5 and p_val < 0.05, set as "patchseq" 
de$diffexpressed[de$avg_log2FC > 2 & de$p_val < 1e-50] <- "patchseq"
# if avg_log2FC < -2.5 and p_val < 0.05, set as "scRNA"
de$diffexpressed[de$avg_log2FC < -2 & de$p_val < 1e-50] <- "scRNA"

v1 <- ggplot(data = de, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed)) + geom_point() + theme_minimal()

v2 <- v1 + geom_vline(xintercept = c(-2, 2), col = "coral2") +
  geom_hline(yintercept = -log10(1e-50), col = "coral2")

mycolors <- c("darkcyan", "coral2", "black")
names(mycolors) <- c("scRNA", "patchseq", "no")
v3 <- v2 + scale_colour_manual(values = mycolors)

# label
de$delabel <- NA

de$delabel[de$diffexpressed != "no"] <- de$gene[de$diffexpressed != "no"]

v4 <- ggplot(data = de, aes(x = avg_log2FC, y=-log10(p_val), col = diffexpressed, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

# plot adding up all layers we have seen so far
volcano <- ggplot(data = de, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed, label = delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = mycolors) +
  geom_vline(xintercept = c(-2, 2), col = "black") +
  geom_hline(yintercept = -log10(1e-50), col="black")


# 3 GO CLUSTER PROFILER ----
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html

## SETUP ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")

library("clusterProfiler")
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggnewscale)
library(topGO)
library(Rgraphviz)


## GO FOR MOUSE DE PATCH VS SCRNA (BP) ----
m_de <- patch_vs_scrna
m_de <- tibble::rownames_to_column(m_de, var = 'gene')
gene.list <- dplyr::filter(m_de, avg_log2FC >= 2) # get positive values for hi vs lo

ego <- enrichGO(gene          = gene.list$gene,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

ego_df <- as.data.frame(ego)
ego_sum <- summary(ego)


## DOTPLOT, ENRICHMAP AND CNETPLOT FOR MOUSE DE PATCH VS SCRNA ----
GO <- dotplot(ego, showCategory = 30)

# 4 PLOT ALL TOGETHER ----
umap + volcano + GO

(umap | GO)/volcano + plot_layout(heights = c(3, 5))


# 5 HYPERGEOMETRIC TEST ----
## SPLIT OLAH MARKERS ----
olah.markers <- read.xlsx('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/olah_suppl_data/microglia_subtype_markers.xls', sheetName = 'degenes_pairwise_microglia_only')

# make compatible with mouse data
olah.markers$gene <- stringr::str_to_title(olah.markers$gene)

# split based on up_type
olah.split <- split(olah.markers, olah.markers$up_type)
new_names <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9")
for (i in 1:length(olah.split)) {
  assign(new_names[i], olah.split[[i]])
}

## ADJUST DE ----
# m_de <- patch_vs_scrna
# m_de <- tibble::rownames_to_column(m_de, var = 'gene')

m_de$fdr <- p.adjust(m_de$p_val, "fdr")

m_de <- filter(m_de, fdr < 0.05) # 1373 genes
m_de <- filter(m_de, avg_log2FC > 0) # to get positive markers only, 190 genes

## FIND q: INNER_JOIN EACH UP_TYPE AND DE BASED ON GENE NAME ----

olah.list <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)

q_list = lapply(olah.list, function(m){
  return(m %>% as.data.frame() %>% inner_join(m_de, by = 'gene', copy = F))
})

new_names <- c("q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8", "q9")
for (i in 1:length(q_list)) {
  assign(new_names[i], q_list[[i]])
}


## FIND n: (OVERLAP BETWEEN OLAH GENES AND OUR MOUSE PATCH SEQ GENES) - m ----
# from Olah et al. 2020, Supplementary Data 3
external_genes = read.xlsx('41467_2020_19737_MOESM5_ESM.xls', sheetName = 'cluster_mean_tpm_20190425') # 25914 genes
external_genes$NA. <- str_to_title(external_genes$NA.)
# from mouse patch seq:
our_genes = read.csv('20200513_Mouse_PatchSeq_Release_cpm.v2.csv') %>% as.data.frame()
our_genes = unlist(our_genes[1])
# find overlapping genes
overlapping_genes = intersect(our_genes, external_genes$NA.)
# get n for hypergeometric testing 
l = length(overlapping_genes) # l = 15942

## PERFORM HYPERGEOMETRIC TEST ----
q = c(0, 0, 0, 3, 1, 3, 1, 1, 4)
m = c(57, 514, 208, 454, 575, 618, 318, 444, 1312)
#n = l - m[i]
k = 190

new_names <- c("h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", "h9")
for (i in 1:length(new_names)) {
  temp <- phyper(q[i], m[i], l - m[i], k, lower.tail = FALSE, log.p = FALSE)
  assign(new_names[i], temp)
}

