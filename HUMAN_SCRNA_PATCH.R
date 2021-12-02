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

scRNA_seurat <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds')


scRNA_meta <- as.data.frame(scRNA_seurat@meta.data)
scRNA_count <- as.data.frame(scRNA_seurat@assays$RNA@counts) # does not work

# used this instead
scRNA_count <- read_csv('/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/matrix.csv')

## PATCH SEQ DATA ----
## set up patch-seq high microglia cells seurat object
# open h_QC_micro_top25_count.csv
patch <- read.csv("/external/rprshnas01/kcni/yjiang/h_QC_micro_top25_count.csv")
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
JackStrawPlot(patch_seurat, dims = 1:20) # pick 15 dims... because dim15 is on the line
ElbowPlot(patch_seurat)

# cluster the cells
patch_seurat <- FindNeighbors(patch_seurat, dims = 1:15)
patch_seurat <- FindClusters(patch_seurat, resolution = 0.5)
head(Idents(patch_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
patch_seurat <- RunUMAP(patch_seurat, dims = 1:15)

DimPlot (patch_seurat, reduction = "umap", label = TRUE) 
DimPlot(patch_seurat, reduction = "umap", pt.size = 1, group.by = "source", label = TRUE)


## SCRNA SEQ DATA ----
## extract microglia cells from metadata and join with count matrix (matrix is already rotated)
scRNA_meta_micro <- filter(scRNA_meta, subclass_label_expanded == 'Microglia')

write.csv(scRNA_meta_micro, file = 'h_QC_scRNA_metadata_microglia.csv')

# extract sample name of metadata
scRNA_micro_samples <- dplyr::select(scRNA_meta_micro, sample_name)

scRNA_trans <- inner_join(scRNA_micro_samples, scRNA_count, by = 'sample_name')

# reshape
sample_id = scRNA_trans$sample_name
rownames(scRNA_trans) = rownames(scRNA_trans) %>% make.names(unique=T)
scRNA = scRNA_trans[-1] %>% t() %>% as.data.frame()
colnames(scRNA) = sample_id

write.csv(scRNA, file = "h_scRNA_microglia_count.csv")

# rm() some files

## make seurat object of microglia cells
scRNA_seurat <- CreateSeuratObject(counts = scRNA, project = "scRNA", min.cells = 3, min.features = 200)

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
JackStrawPlot(scRNA_seurat, dims = 1:20)
ElbowPlot(scRNA_seurat)

# cluster the cells
scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:20)
scRNA_seurat <- FindClusters(scRNA_seurat, resolution = 0.5)
head(Idents(scRNA_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
scRNA_seurat <- RunUMAP(scRNA_seurat, dims = 1:20)

# there's only one cluster, which is okay considering the fact that there are only 69 cells
DimPlot (scRNA_seurat, reduction = "umap", label = TRUE) 
DimPlot(scRNA_seurat, reduction = "umap", pt.size = 1, group.by = "source", label = TRUE)

## COMBINE AND SPLIT SEURAT OBJECTS ----
patch_scRNA_combined = merge(patch_seurat, y = scRNA_seurat, add.cell.ids = c("patchseq", "scRNA"), project = "source")

# split the dataset into a list of two seurat objects (patchseq and scRNA in terms of source)
microglia.list <- SplitObject(patch_scRNA_combined, split.by = "source")
reference.list <- microglia.list[c("patchseq", "scRNA")]


## PERFORM INTEGRATION ----
microglia.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)
microglia.combined <- IntegrateData(anchorset = microglia.anchors, dims = 1:20, k.weight = 65)


## PERFORM AN INTEGRATED ANALYSIS ----
DefaultAssay(microglia.combined) <- "integrated"

microglia.combined <- ScaleData(microglia.combined, verbose = FALSE)
microglia.combined <- RunPCA(microglia.combined, npcs = 20, verbose = FALSE)
microglia.combined <- RunUMAP(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindNeighbors(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindClusters(microglia.combined, resolution = 0.5)

saveRDS(microglia.combined, file = "/external/rprshnas01/kcni/yjiang/h_patchseq_scRNA.rds")

# if want to read it
microglia.combined <- readRDS('h_patchseq_scRNA.rds')

# visualization
umap <- DimPlot(microglia.combined, reduction = "umap", group.by = "source")


## FIND DIFFERENTIALLY EXPRESSED GENES ----
DefaultAssay(microglia.combined) <- "integrated"
Idents(microglia.combined) <- "source"

patch_vs_scrna <- FindMarkers(microglia.combined, ident.1 = "patchseq", ident.2 = "scRNA", verbose = FALSE)
write.xlsx(patch_vs_scrna, file = "h_patch_scRNA_de.xlsx", sheetName = "h_patch_scRNA_de")


# 2 VOLCANO ----
# https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

de <- patch_vs_scrna
de <- tibble::rownames_to_column(de, var = "gene")

v0 <- ggplot(data = de, aes(x = avg_log2FC, y = p_val)) + geom_point()

# add a column of NAs
de$diffexpressed <- "no"
# if avg_log2FC > 2 and p_val < 0.05, set as "patchseq" 
de$diffexpressed[de$avg_log2FC > 2 & de$p_val < 0.05] <- "patchseq"
# if avg_log2FC < -2 and p_val < 0.05, set as "scRNA"
de$diffexpressed[de$avg_log2FC < -2 & de$p_val < 0.05] <- "scRNA"

v1 <- ggplot(data = de, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed)) + geom_point() + theme_minimal()

v2 <- v1 + geom_vline(xintercept = c(-2, 2), col = "coral2") +
  geom_hline(yintercept = -log10(0.05), col = "coral2")

mycolors <- c("darkcyan", "coral2", "black")
names(mycolors) <- c("scRNA", "patchseq", "no")
v3 <- v2 + scale_colour_manual(values = mycolors)

# label
de$delabel <- NA

de$delabel[de$diffexpressed != "no"] <- de$gene[de$diffexpressed != "no"]

de$delabel[de$avg_log2FC > -2.5 & de$avg_log2FC < -2] <- NA

# plot adding up all layers we have seen so far
volcano <- ggplot(data = de, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed, label = delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size = 3, colour = 'black') +
  scale_color_manual(values = mycolors) +
  geom_vline(xintercept = c(-2, 2), col = "black") +
  geom_hline(yintercept = -log10(0.05), col="black")


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


## GO FOR HUMAN DE PATCH VS SCRNA (BP) ----
h_de <- patch_vs_scrna
h_de <- tibble::rownames_to_column(m_de, var = 'gene')
gene.list <- dplyr::filter(h_de, avg_log2FC >= 2) # get positive values for hi vs lo

ego <- enrichGO(gene          = gene.list$gene,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

ego_df <- as.data.frame(ego)
ego_sum <- summary(ego)


## DOTPLOT, ENRICHMAP AND CNETPLOT FOR MOUSE DE PATCH VS SCRNA ----
GO <- dotplot(ego, showCategory = 20, font.size = 7)

# 4 PLOT ALL TOGETHER ----
umap + volcano + GO

(umap + GO)/volcano + plot_layout(heights = c(1, 1))
# (umap + volcano + plot_layout(ncol = 2, widths = c(1, 3)))/GO + plot_layout(heights = c(1, 2))


# 5 HYPERGEOMETRIC TEST ----
## SPLIT OLAH MARKERS ----
olah.markers <- read.xlsx('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/olah_suppl_data/microglia_subtype_markers.xls', sheetName = 'degenes_pairwise_microglia_only')

# split based on up_type
olah.split <- split(olah.markers, olah.markers$up_type)
new_names <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9")
for (i in 1:length(olah.split)) {
  assign(new_names[i], olah.split[[i]])
}

## ADJUST DE ----
# h_de <- patch_vs_scrna
# h_de <- tibble::rownames_to_column(m_de, var = 'gene')

h_de$fdr <- p.adjust(h_de$p_val, "fdr")

h_de <- filter(h_de, fdr < 0.05) # 1234 genes
h_de <- filter(h_de, avg_log2FC > 0) # to get positive markers only, 287 genes

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
setwd("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/olah_suppl_data")
external_genes = read.xlsx('41467_2020_19737_MOESM5_ESM.xls', sheetName = 'cluster_mean_tpm_20190425')
# from Berg et al. 2020: https://github.com/keon-arbabi/berg_patchseq/blob/main/data/20200512_Human_PatchSeq_Release_count.csv
our_genes = fread("https://raw.githubusercontent.com/keon-arbabi/berg_patchseq/main/data/20200512_Human_PatchSeq_Release_cpm.csv?token=AUB7ZHZIO5RYZNN74WJSVILA2M3UO") %>% as.data.frame()
our_genes = unlist(our_genes[1])
# find overlapping genes
overlapping_genes = intersect(our_genes, external_genes$NA.)
# get n for hypergeometric testing 
l = length(overlapping_genes) # 21698

## PERFORM HYPERGEOMETRIC TEST ----
q = c(3, 19, 23, 52, 38, 21, 9, 11, 80)
m = c(57, 514, 208, 454, 575, 618, 318, 444, 1312)
#n = l - m[i]
k = 287

new_names <- c("h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", "h9")
for (i in 1:length(new_names)) {
  temp <- phyper(q[i], m[i], l - m[i], k, lower.tail = FALSE, log.p = FALSE)
  assign(new_names[i], temp)
}