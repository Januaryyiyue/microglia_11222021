## SETUP ----

library(Seurat)
library(dplyr)
library(xlsx)

# current working directory is: /external/rprshnas01/kcni/yjiang
setwd('/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects')

scRNA_seurat <- readRDS('Seu_AIBS_obj_update_07JUN21.rds')

scRNA_meta <- as.data.frame(scRNA_seurat@meta.data)
scRNA_count <- as.data.frame(scRNA_seurat@assays$RNA@counts) # DOES NOT WORK

# use this instead
scRNA_count <- read_csv('/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/matrix.csv')

## PATCH SEQ DATA ----
## set up patch-seq high microglia cells seurat object
# open h_micro_top25_count copy.csv (from my laptop)
patch <- read.csv("/external/rprshnas01/kcni/yjiang/h_micro_top25_count copy.csv")
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
JackStrawPlot(patch_seurat, dims = 1:20) # pick 17 dims... because dim17 is on the line
ElbowPlot(patch_seurat)

# cluster the cells
patch_seurat <- FindNeighbors(patch_seurat, dims = 1:17)
patch_seurat <- FindClusters(patch_seurat, resolution = 0.5)
head(Idents(patch_seurat), 5)

reticulate::py_install(packages = 'umap-learn')
patch_seurat <- RunUMAP(patch_seurat, dims = 1:17)

# there's only one cluster, which is okay considering the fact that there are only 69 cells
DimPlot (patch_seurat, reduction = "umap", label = TRUE) 
DimPlot(patch_seurat, reduction = "umap", pt.size = 1, group.by = "source", label = TRUE)


## SCRNA SEQ DATA ----
## extract microglia cells from metadata and join with count matrix (matrix is already rotated)
scRNA_meta_micro <- filter(scRNA_meta, subclass_label_expanded == 'Microglia')

setwd('/external/rprshnas01/kcni/yjiang')
write.csv(scRNA_meta_micro, file = 'scRNA_metadata_microglia.csv')

# extract sample name of metadata
scRNA_micro_samples <- dplyr::select(scRNA_meta_micro, sample_name)

scRNA_trans <- inner_join(scRNA_micro_samples, scRNA_count, by = 'sample_name')

# reshape
sample_id = scRNA_trans$sample_name
rownames(scRNA_trans) = rownames(scRNA_trans) %>% make.names(unique=T)
scRNA = scRNA_trans[-1] %>% t() %>% as.data.frame()
colnames(scRNA) = sample_id

setwd('/external/rprshnas01/kcni/yjiang')
write.csv(scRNA, file = "scRNA_microglia_count.csv")

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
JackStrawPlot(scRNA_seurat, dims = 1:20) # pick 20 dims... because dim17 is on the line
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
microglia.combined <- IntegrateData(anchorset = microglia.anchors, dims = 1:20)


## PERFORM AN INTEGRATED ANALYSIS ----
DefaultAssay(microglia.combined) <- "integrated"

microglia.combined <- ScaleData(microglia.combined, verbose = FALSE)
microglia.combined <- RunPCA(microglia.combined, npcs = 20, verbose = FALSE)
microglia.combined <- RunUMAP(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindNeighbors(microglia.combined, reduction = "pca", dims = 1:20)
microglia.combined <- FindClusters(microglia.combined, resolution = 0.5)

saveRDS(microglia.combined, file = "/external/rprshnas01/kcni/yjiang/patchseq_scRNA_human.rds")

# if want to read it
microglia.combined <- readRDS('patchseq_scRNA_human.rds') 

# Visualization
p1 <- DimPlot(microglia.combined, reduction = "umap", group.by = "source")
p2 <- DimPlot(microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(microglia.combined, reduction = "umap", split.by = "source")

## IDENTIFY CONSERVED CELL TYPE MARKERS ----
DefaultAssay(microglia.combined) <- "RNA"
all.markers <- FindAllMarkers(microglia.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
kable(all.markers[1:10,])

setwd('/external/rprshnas01/kcni/yjiang')
write.csv(all.markers, file = 'patchseq_scRNA_clusters_markers_human.csv')

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DefaultAssay(microglia.combined) <- "integrated"
DoHeatmap(microglia.combined, features = top10$gene, size = 4)

FeaturePlot(microglia.combined, features = c("C1QC", "C1QA", "AIF1", "APOE", "TREM2", "ITGAX"), min.cutoff = "q9")

## IDENTIFY DIFERENTIAL EXPRESSED GENES ACROSS CONDITIONS ----

library(ggplot2)
library(cowplot)
library(ggrepel)
theme_set(theme_cowplot())

avg.microglia.combined <- as.data.frame(log1p(AverageExpression(microglia.combined, verbose = FALSE)$RNA))
avg.microglia.combined$gene <- rownames(avg.microglia.combined)

write.csv(avg.microglia.combined, file = "clusters_expression.csv")

cluster_0_vs_1 <- FindMarkers(microglia.combined, ident.1 = "0", ident.2 = "1", verbose = FALSE)
cluster_0_vs_2 <- FindMarkers(microglia.combined, ident.1 = "0", ident.2 = "2", verbose = FALSE)
cluster_1_vs_2 <- FindMarkers(microglia.combined, ident.1 = "1", ident.2 = "2", verbose = FALSE)
cluster_2_vs_1 <- FindMarkers(microglia.combined, ident.1 = "2", ident.2 = "1", verbose = FALSE)
cluster_1_vs_rest <- FindMarkers(microglia.combined, ident.1 = "1", verbose = FALSE)
cluster_2_vs_rest <- FindMarkers(microglia.combined, ident.1 = "2", verbose = FALSE)

# scRNA cells only (for GO analysis later)
microglia.combined$celltype.source <- paste(Idents(microglia.combined), microglia.combined$source, sep = "_")
microglia.combined$celltype <- Idents(microglia.combined)
Idents(microglia.combined) <- "celltype.source"
scRNA_1_vs_2 <- FindMarkers(microglia.combined, ident.1 = "1_scRNA", ident.2 = "2_scRNA", verbose = FALSE)
View(scRNA_1_vs_2)


write.xlsx(cluster_0_vs_1, file = "human_cluster_0_vs_1_DE.xlsx", sheetName = "human_cluster_0_vs_1_DE")
write.xlsx(cluster_0_vs_2, file = "human_cluster_0_vs_2_DE.xlsx", sheetName = "human_cluster_0_vs_2_DE")
write.xlsx(cluster_1_vs_2, file = "human_scluster_1_vs_2_DE.xlsx", sheetName = "human_cluster_1_vs_2_DE")
write.xlsx(cluster_2_vs_1, file = "human_cluster_2_vs_1_DE.xlsx", sheetName = "human_cluster_2_vs_1_DE")
write.xlsx(cluster_1_vs_rest, file = "human_cluster_1_vs_rest_DE.xlsx", sheetName = "human_cluster_1_vs_rest_DE")
write.xlsx(cluster_2_vs_rest, file = "human_cluster_2_vs_rest_DE.xlsx", sheetName = "human_cluster_2_vs_rest_DE")
write.xlsx(scRNA_1_vs_2, file = 'human_scRNA_1_vs_2_DE.xlsx', sheetName = 'human_scRNA_1_vs_2_DE')

# string to title for MGEnrichment
str_0_1 = tibble::rownames_to_column(cluster_0_vs_1)
str_0_1$rowname = stringr::str_to_title(str_0_1$rowname)

str_0_2 = tibble::rownames_to_column(cluster_0_vs_2)
str_0_2$rowname = stringr::str_to_title(str_0_2$rowname)

str_1_2 = tibble::rownames_to_column(cluster_1_vs_2)
str_1_2$rowname = stringr::str_to_title(str_1_2$rowname)

write.xlsx(str_0_1, file = "str_0_1.xlsx", sheetName = "str_0_1")
write.xlsx(str_0_2, file = "str_0_2.xlsx", sheetName = "str_0_2")
write.xlsx(str_1_2, file = "str_1_2.xlsx", sheetName = "str_1_2")


options(ggrepel.max.overlaps = Inf)

# scatterplots
avg.microglia.combined <- avg.microglia.combined %>% dplyr::rename("cluster0" = "0") %>% dplyr::rename("cluster1" = "1") %>% dplyr::rename("cluster2" = "2")
p_0_1 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("cluster 0 vs 1")
p_0_2 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("cluster 0 vs 2")
p_1_2 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("cluster 1 vs 2")

# scatterplots with OLAH subclass markers
m1_01 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_1") + geom_label_repel(aes( label = ifelse((cluster0/cluster1 > 2 | cluster1/cluster0 > 2) & gene %in% m1$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m1_02 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_2") + geom_label_repel(aes( label = ifelse((cluster0/cluster2 > 2 | cluster2/cluster0 > 2) & gene %in% m1$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m1_12 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_1_vs_2") + geom_label_repel(aes( label = ifelse((cluster1/cluster2 > 2 | cluster2/cluster1 > 2) & gene %in% m1$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")

m2_01 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_1") + geom_label_repel(aes( label = ifelse((cluster0/cluster1 > 2 | cluster1/cluster0 > 2) & gene %in% m2$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m2_02 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_2") + geom_label_repel(aes( label = ifelse((cluster0/cluster2 > 2 | cluster2/cluster0 > 2) & gene %in% m2$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m2_12 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_1_vs_2") + geom_label_repel(aes( label = ifelse((cluster1/cluster2 > 2 | cluster2/cluster1 > 2) & gene %in% m2$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")

m3_01 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_1") + geom_label_repel(aes( label = ifelse((cluster0/cluster1 > 2 | cluster1/cluster0 > 2) & gene %in% m3$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m3_02 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_2") + geom_label_repel(aes( label = ifelse((cluster0/cluster2 > 2 | cluster2/cluster0 > 2) & gene %in% m3$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m3_12 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_1_vs_2") + geom_label_repel(aes( label = ifelse((cluster1/cluster2 > 2 | cluster2/cluster1 > 2) & gene %in% m3$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")

m4_01 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_1") + geom_label_repel(aes( label = ifelse((cluster0/cluster1 > 2 | cluster1/cluster0 > 2) & gene %in% m4$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m4_02 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_2") + geom_label_repel(aes( label = ifelse((cluster0/cluster2 > 2 | cluster2/cluster0 > 2) & gene %in% m4$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m4_12 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_1_vs_2") + geom_label_repel(aes( label = ifelse((cluster1/cluster2 > 2 | cluster2/cluster1 > 2) & gene %in% m4$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")

m5_01 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_1") + geom_label_repel(aes( label = ifelse((cluster0/cluster1 > 2 | cluster1/cluster0 > 2) & gene %in% m5$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m5_02 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_2") + geom_label_repel(aes( label = ifelse((cluster0/cluster2 > 2 | cluster2/cluster0 > 2) & gene %in% m5$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m5_12 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_1_vs_2") + geom_label_repel(aes( label = ifelse((cluster1/cluster2 > 2 | cluster2/cluster1 > 2) & gene %in% m5$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")

m6_01 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_1") + geom_label_repel(aes( label = ifelse((cluster0/cluster1 > 2 | cluster1/cluster0 > 2) & gene %in% m6$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m6_02 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_2") + geom_label_repel(aes( label = ifelse((cluster0/cluster2 > 2 | cluster2/cluster0 > 2) & gene %in% m6$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m6_12 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_1_vs_2") + geom_label_repel(aes( label = ifelse((cluster1/cluster2 > 2 | cluster2/cluster1 > 2) & gene %in% m6$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")

m7_01 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_1") + geom_label_repel(aes( label = ifelse((cluster0/cluster1 > 2 | cluster1/cluster0 > 2) & gene %in% m7$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m7_02 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_2") + geom_label_repel(aes( label = ifelse((cluster0/cluster2 > 2 | cluster2/cluster0 > 2) & gene %in% m7$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m7_12 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_1_vs_2") + geom_label_repel(aes( label = ifelse((cluster1/cluster2 > 2 | cluster2/cluster1 > 2) & gene %in% m7$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")

m8_01 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_1") + geom_label_repel(aes( label = ifelse((cluster0/cluster1 > 2 | cluster1/cluster0 > 2) & gene %in% m8$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m8_02 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_2") + geom_label_repel(aes( label = ifelse((cluster0/cluster2 > 2 | cluster2/cluster0 > 2) & gene %in% m8$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m8_12 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_1_vs_2") + geom_label_repel(aes( label = ifelse((cluster1/cluster2 > 2 | cluster2/cluster1 > 2) & gene %in% m8$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")

m9_01 <- ggplot(avg.microglia.combined, aes(cluster0, cluster1)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_1") + geom_label_repel(aes( label = ifelse((cluster0/cluster1 > 2 | cluster1/cluster0 > 2) & gene %in% m9$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m9_02 <- ggplot(avg.microglia.combined, aes(cluster0, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_0_vs_2") + geom_label_repel(aes( label = ifelse((cluster0/cluster2 > 2 | cluster2/cluster0 > 2) & gene %in% m9$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")
m9_12 <- ggplot(avg.microglia.combined, aes(cluster1, cluster2)) + geom_point() + ggtitle("Olah_subtype1_cluster_1_vs_2") + geom_label_repel(aes( label = ifelse((cluster1/cluster2 > 2 | cluster2/cluster1 > 2) & gene %in% m9$gene,as.character(gene),'')), box.padding = 0.2, point.padding = 0.2) + theme_classic() + theme(legend.position="none")


## HYPERGEOMETRIC TEST ----

# m
m1 <- read.xlsx('m1_de.xlsx', sheetName = 'm1_de') #57
m2 <- read.xlsx('m2_de.xlsx', sheetName = 'm2_de') #514
m3 <- read.xlsx('m3_de.xlsx', sheetName = 'm3_de') #208
m4 <- read.xlsx('m4_de.xlsx', sheetName = 'm4_de') #454
m5 <- read.xlsx('m5_de.xlsx', sheetName = 'm5_de') #575
m6 <- read.xlsx('m6_de.xlsx', sheetName = 'm6_de') #618
m7 <- read.xlsx('m7_de.xlsx', sheetName = 'm7_de') #318
m8 <- read.xlsx('m8_de.xlsx', sheetName = 'm8_de') #444
m9 <- read.xlsx('m9_de.xlsx', sheetName = 'm9_de') #1312

# n = 21698, from previous hypergeometric test (patch and scRNA have the same number of genes)

# k (filter fdr < 0.05 and logFC > 0)

# 0 vs 1: 299
cluster_0_vs_1 <- read.xlsx('human_cluster_0_vs_1_DE.xlsx', sheetName = 'human_cluster_0_vs_1_DE')
cluster_0_vs_1$fdr <- p.adjust(cluster_0_vs_1$p_val, "fdr")
cluster_0_vs_1 <- filter(cluster_0_vs_1, fdr < 0.05)
cluster_0_vs_1 <- filter(cluster_0_vs_1, avg_log2FC > 0)
cluster_0_vs_1 <- dplyr::rename(cluster_0_vs_1, 'gene' = 'NA.')
write.xlsx(cluster_0_vs_1, file = "human_cluster_0_vs_1_FDR.xlsx", sheetName = "human_cluster_0_vs_1_FDR", append = FALSE)

# 0 vs 2: 141
cluster_0_vs_2 <- read.xlsx('human_cluster_0_vs_2_DE.xlsx', sheetName = 'human_cluster_0_vs_2_DE')
cluster_0_vs_2$fdr <- p.adjust(cluster_0_vs_2$p_val, "fdr")
cluster_0_vs_2 <- filter(cluster_0_vs_2, fdr < 0.05)
cluster_0_vs_2 <- filter(cluster_0_vs_2, avg_log2FC > 0)
cluster_0_vs_2 <- dplyr::rename(cluster_0_vs_2, 'gene' = 'NA.')
write.xlsx(cluster_0_vs_2, file = "human_cluster_0_vs_2_FDR.xlsx", sheetName = "human_cluster_0_vs_2_FDR", append = FALSE)

# 1 vs 2: 178
cluster_1_vs_2 <- read.xlsx('human_cluster_1_vs_2_DE.xlsx', sheetName = 'human_cluster_1_vs_2_DE')
cluster_1_vs_2$fdr <- p.adjust(cluster_1_vs_2$p_val, "fdr")
cluster_1_vs_2 <- filter(cluster_1_vs_2, fdr < 0.05)
cluster_1_vs_2 <- filter(cluster_1_vs_2, avg_log2FC > 0)
cluster_1_vs_2 <- dplyr::rename(cluster_1_vs_2, 'gene' = 'NA.')
write.xlsx(cluster_1_vs_2, file = "human_cluster_1_vs_2_FDR.xlsx", sheetName = "human_cluster_1_vs_2_FDR", append = FALSE)

# q
x1 <- dplyr::inner_join(m1, cluster_0_vs_1, by = 'gene', copy = FALSE) # 0
x2 <- dplyr::inner_join(m2, cluster_0_vs_1, by = 'gene', copy = FALSE) # 2
x3 <- dplyr::inner_join(m3, cluster_0_vs_1, by = 'gene', copy = FALSE) # 1
x4 <- dplyr::inner_join(m4, cluster_0_vs_1, by = 'gene', copy = FALSE) # 9
x5 <- dplyr::inner_join(m5, cluster_0_vs_1, by = 'gene', copy = FALSE) # 3
x6 <- dplyr::inner_join(m6, cluster_0_vs_1, by = 'gene', copy = FALSE) # 0
x7 <- dplyr::inner_join(m7, cluster_0_vs_1, by = 'gene', copy = FALSE) # 0
x8 <- dplyr::inner_join(m8, cluster_0_vs_1, by = 'gene', copy = FALSE) # 0
x9 <- dplyr::inner_join(m9, cluster_0_vs_1, by = 'gene', copy = FALSE) # 3

y1 <- dplyr::inner_join(m1, cluster_0_vs_2, by = 'gene', copy = FALSE) # 0
y2 <- dplyr::inner_join(m2, cluster_0_vs_2, by = 'gene', copy = FALSE) # 1
y3 <- dplyr::inner_join(m3, cluster_0_vs_2, by = 'gene', copy = FALSE) # 2
y4 <- dplyr::inner_join(m4, cluster_0_vs_2, by = 'gene', copy = FALSE) # 5
y5 <- dplyr::inner_join(m5, cluster_0_vs_2, by = 'gene', copy = FALSE) # 5
y6 <- dplyr::inner_join(m6, cluster_0_vs_2, by = 'gene', copy = FALSE) # 1
y7 <- dplyr::inner_join(m7, cluster_0_vs_2, by = 'gene', copy = FALSE) # 0
y8 <- dplyr::inner_join(m8, cluster_0_vs_2, by = 'gene', copy = FALSE) # 1
y9 <- dplyr::inner_join(m9, cluster_0_vs_2, by = 'gene', copy = FALSE) # 0

z1 <- dplyr::inner_join(m1, cluster_1_vs_2, by = 'gene', copy = FALSE) # 0
z2 <- dplyr::inner_join(m2, cluster_1_vs_2, by = 'gene', copy = FALSE) # 4
z3 <- dplyr::inner_join(m3, cluster_1_vs_2, by = 'gene', copy = FALSE) # 3
z4 <- dplyr::inner_join(m4, cluster_1_vs_2, by = 'gene', copy = FALSE) # 5
z5 <- dplyr::inner_join(m5, cluster_1_vs_2, by = 'gene', copy = FALSE) # 7
z6 <- dplyr::inner_join(m6, cluster_1_vs_2, by = 'gene', copy = FALSE) # 7
z7 <- dplyr::inner_join(m7, cluster_1_vs_2, by = 'gene', copy = FALSE) # 2
z8 <- dplyr::inner_join(m8, cluster_1_vs_2, by = 'gene', copy = FALSE) # 2
z9 <- dplyr::inner_join(m9, cluster_1_vs_2, by = 'gene', copy = FALSE) # 6

# hypergeometric test 0 vs 1
phyper(0, 57, 21698-57, 299, lower.tail = FALSE, log.p = FALSE) # 0.547042
phyper(2, 500, 21698-500, 299, lower.tail = FALSE, log.p = FALSE) # 0.9699549
phyper(1, 208, 21698-208, 299, lower.tail = FALSE, log.p = FALSE) # 0.7835308
phyper(9, 454, 21698-454, 299, lower.tail = FALSE, log.p = FALSE) # 0.09894951
phyper(3, 500, 21698-500, 299, lower.tail = FALSE, log.p = FALSE) # 0.9162837
phyper(0, 500, 21698-500, 299, lower.tail = FALSE, log.p = FALSE) # 0.9991058
phyper(0, 318, 21698-318, 299, lower.tail = FALSE, log.p = FALSE) # 0.9882666
phyper(0, 444, 21698-444, 299, lower.tail = FALSE, log.p = FALSE) # 0.9980209
phyper(3, 500, 21698-500, 299, lower.tail = FALSE, log.p = FALSE) # 0.9162837

# hypergeometric test 0 vs 2
phyper(0, 57, 21698-57, 141, lower.tail = FALSE, log.p = FALSE) # 0.3107102
phyper(1, 500, 21698-500, 141, lower.tail = FALSE, log.p = FALSE) # 0.8393059
phyper(2, 208, 21698-208, 141, lower.tail = FALSE, log.p = FALSE) # 0.1536896
phyper(5, 454, 21698-454, 141, lower.tail = FALSE, log.p = FALSE) # 0.07621712
phyper(5, 500, 21698-500, 141, lower.tail = FALSE, log.p = FALSE) # 0.1080116
phyper(1, 500, 21698-500, 141, lower.tail = FALSE, log.p = FALSE) # 0.8393059
phyper(0, 318, 21698-318, 141, lower.tail = FALSE, log.p = FALSE) # 0.8761328
phyper(1, 444, 21698-444, 141, lower.tail = FALSE, log.p = FALSE) # 0.7871565
phyper(0, 500, 21698-500, 141, lower.tail = FALSE, log.p = FALSE) # 0.963041

# hypergeometric test 1 vs 2
phyper(0, 57, 21698-57, 178, lower.tail = FALSE, log.p = FALSE) # 0.3750844
phyper(4, 514, 21698-514, 178, lower.tail = FALSE, log.p = FALSE) # 0.4140551
phyper(3, 208, 21698-208, 178, lower.tail = FALSE, log.p = FALSE) # 0.09239406
phyper(5, 454, 21698-454, 178, lower.tail = FALSE, log.p = FALSE) # 0.1706818
phyper(7, 575, 21698-575, 178, lower.tail = FALSE, log.p = FALSE) # 0.1018262
phyper(7, 618, 21698-618, 178, lower.tail = FALSE, log.p = FALSE) # 0.1368832
phyper(2, 318, 21698-318, 178, lower.tail = FALSE, log.p = FALSE) # 0.4855007
phyper(2, 444, 21698-444, 178, lower.tail = FALSE, log.p = FALSE) # 0.7088
phyper(6, 1312, 21698-1312, 178, lower.tail = FALSE, log.p = FALSE) # 0.9186179
