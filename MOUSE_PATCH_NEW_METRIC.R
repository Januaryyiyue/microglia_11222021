# 3 SEURAT INTEGRATION + NEW QCMETRIC ----

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
mouse_data <- read.csv("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/mouse_patch_seq/data/20200513_Mouse_PatchSeq_Release_count.v2.csv")

mouse_data = mouse_data %>% distinct(X, .keep_all = T)
transcript = mouse_data$X
rownames(mouse_data) = rownames(mouse_data) %>% make.names(unique=T)
mouse_data_trans = mouse_data[-1] %>% t() %>% as.data.frame()
colnames(mouse_data_trans) = transcript
mouse_data_trans # now we have the count matrix rotated 90 degrees

# open qcmetrics output
qc <- read.csv("https://raw.githubusercontent.com/keon-arbabi/patch-seq-microglia/main/output/qcMetrics.csv")
qc$sample_id = sub("-",".", qc$sample_id)
qc$sample_id = sub("-",".", qc$sample_id)
microglia_ranked <- qc %>% select(c(sample_id, Macrophage)) %>% arrange(desc(Macrophage))

top20 <- slice_max(microglia_ranked, order_by = Macrophage, prop = 0.20) # microglianess cutoff is 0.3556218
bottom20 <- slice_min(microglia_ranked, order_by = Macrophage, prop = 0.20) # microglianess cutoff is 0.01389722

microglia_ranked$group <- "low"
microglia_ranked$group[microglia_ranked$Macrophage > 0.01389722 & microglia_ranked$Macrophage < 0.3556218] <- "med"
microglia_ranked$group[microglia_ranked$Macrophage >= 0.3556218] <- "high"

write.xlsx(microglia_ranked, file = "./output/m_QC_microglia_ranked.xlsx",
           sheetName = "m_QC_microglia_ranked", append = FALSE)

top_20 <- filter(microglia_ranked, group == 'high')
bottom_20 <- filter(microglia_ranked, group == 'low')

top_20 <- dplyr::rename(top_20, transcriptomics_sample_id = sample_id)
bottom_20 <- dplyr::rename(bottom_20, transcriptomics_sample_id = sample_id)

write.xlsx(top_20, file = "./output/m_QC_top_20_microglia.xlsx",
           sheetName = "m_QC_top_20_microglia", append = FALSE)
write.xlsx(bottom_20, file = "./output/m_QC_bottom_20_microglia.xlsx",
           sheetName = "m_QC_bottom_20_microglia", append = FALSE)

# inner_join the transcriptomics_sample_id of top_20 and bottom_20 with mouse_data_trans, respectively (save as csv)
top_20_id <- dplyr::select(top_20, transcriptomics_sample_id)
bottom_20_id <- dplyr::select(bottom_20, transcriptomics_sample_id)

mouse_data_trans <- mouse_data_trans %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = 'transcriptomics_sample_id')

micro_top20 <- inner_join(mouse_data_trans, top_20_id, by = 'transcriptomics_sample_id', copy = T)
micro_bottom20 <- inner_join(mouse_data_trans, bottom_20_id, by = 'transcriptomics_sample_id', copy = T)

write.csv(micro_top20, file = "./output/m_QC_top_20_micro_count_trans.csv")
write.csv(micro_bottom20, file = "./output/m_QC_bottom_20_micro_count_trans.csv")

# reshape micro_top20 and micro_bottom20 and save them as xlsx files (these are the files to work with for the integration)
micro_top20 <- micro_top20 %>% as.data.frame()
micro_top20 = micro_top20 %>% distinct(transcriptomics_sample_id, .keep_all = T)
transcriptomics_sample_id = micro_top20$transcriptomics_sample_id
rownames(micro_top20) = rownames(micro_top20) %>% make.names(unique=T)
mouse_t = micro_top20[-1] %>% t() %>% as.data.frame()
colnames(mouse_t) = transcriptomics_sample_id

micro_bottom20 <- micro_bottom20 %>% as.data.frame()
micro_bottom20 = micro_bottom20 %>% distinct(transcriptomics_sample_id, .keep_all = T)
transcriptomics_sample_id = micro_bottom20$transcriptomics_sample_id
rownames(micro_bottom20) = rownames(micro_bottom20) %>% make.names(unique=T)
mouse_b = micro_bottom20[-1] %>% t() %>% as.data.frame()
colnames(mouse_b) = transcriptomics_sample_id

write.csv(mouse_t, file = "./output/m_QC_micro_top20_count.csv")
write.csv(mouse_b, file = "./output/m_QC_micro_bottom20_count.csv")


## JOIN METADATA + TTYPE WITH TOP20 AND BOTTOM20 ----
mouse_top <- select(micro_top20, transcriptomics_sample_id)
mouse_bottom <- select(micro_bottom20, transcriptomics_sample_id)

# open metadata and add ttype
meta <- read.csv("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/mouse_patch_seq/data/20200625_patchseq_metadata_mouse.csv")

meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)

# inner_join top20 and bottom20 with meta; save as csv
top20sub <- inner_join(mouse_top, meta, by = 'transcriptomics_sample_id', copy = T)
bottom20sub <- inner_join(mouse_bottom, meta, by = 'transcriptomics_sample_id', copy = T)

write.csv(top20sub, file = "./output/m_QC_metadata_top20.csv")
write.csv(bottom20sub, file = "./output/m_QC_metadata_bottom20.csv")

# remove everything to clear up space: rm(list = ls())


## SET UP HIGH MICROGLIANESS SEURAT OBJECT AND JOIN WITH METADATA ----
# open m_micro_top20_count.csv 
top <- read.csv('./output/m_QC_micro_top20_count.csv')
rownames(top) <- top$X

top_seurat <- CreateSeuratObject(counts = top, project = "high", min.cells = 3, min.features = 200)

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
bottom <- read.csv('./output/m_QC_micro_bottom20_count.csv')
rownames(bottom) <- bottom$X

bottom_seurat <- CreateSeuratObject(counts = bottom, project = "low", min.cells = 3, min.features = 200)

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

saveRDS(microglia.combined, file = "./output/m_QC_microglia_combined.rds")

# if want to read it
microglia.combined <- readRDS('m_QC_microglia_combined.rds')

# visualization
p1 <- DimPlot(microglia.combined, reduction = "umap", group.by = "microglianess")
p2 <- DimPlot(microglia.combined, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(microglia.combined, reduction = "umap", split.by = "microglianess")
p1 # will included in the manuscript

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
# use FindMarkers to find DE genes using Wilcoxon Rank Sum Test from limma
DefaultAssay(microglia.combined) <- "RNA"
Idents(microglia.combined) <- "microglianess" # assay should be RNA
hi_vs_lo <- FindMarkers(microglia.combined, ident.1 = 'high', ident.2 = 'low', verbose = FALSE)

write.xlsx(hi_vs_lo, file= "./output/m_QC_de.xlsx", sheetName = "m_QC_de", append = FALSE)

# visualize DE using feature plot
FeaturePlot(microglia.combined, features = c("C1qc", "C1qa", "Aif1"), split.by = "microglianess", max.cutoff = 3, cols = c("grey", "red"))

# visualize DE using violin plot
DefaultAssay(microglia.combined) <- "RNA"
microglia.combined$celltype.microglianess <- paste(Idents(microglia.combined), microglia.combined$microglianess, sep = "_")
microglia.combined$celltype <- Idents(microglia.combined)
Idents(microglia.combined) <- "celltype.microglianess"
violin <- VlnPlot(microglia.combined, features = c("C1qc", "C1qa", "Aif1"), split.by = "microglianess", group.by = "celltype", pt.size = 0, combine = FALSE)

# top vs bottom heatmap
DefaultAssay(microglia.combined) <- "integrated"
Idents(microglia.combined) <- "microglianess" # assay should be integrated here
mouse.markers <- FindAllMarkers(microglia.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- mouse.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
heatmap <- DoHeatmap(microglia.combined, features = top10$gene, group.bar = TRUE, size = 5, angle = 0) + scale_fill_gradientn(colors = c("lightsteelblue3", "white", "salmon2")) + theme(text = element_text(size = 12))

((p1 / heatmap) + plot_layout(ncol=1, heights=c(3, 4)) | wrap_plots(plots = violin, ncol = 1) ) + plot_layout(ncol=2, widths=c(3, 2))


# 5 REGRESSION ----
library(dplyr)
library(ggplot2)
library(xlsx)
library(stats)
library(tidyverse)
library(cowplot)

## OPEN METADATA AND ADD TTYPE + MICROGLIANESS
meta <- read.csv("20200625_patchseq_metadata_mouse.csv")
meta$subclass <- gsub( "\\s.*", "", meta$corresponding_AIT2.3.1_alias) # take everything before first space

meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)

write.csv(meta, file = "m_QC_metadata_ttype.csv")

micro <- read.xlsx('m_QC_microglia_ranked.xlsx', sheetName = 'm_QC_microglia_ranked')
micro <- dplyr::rename(micro, transcriptomics_sample_id = sample_id)
data <- inner_join(micro, meta, by = "transcriptomics_sample_id")

# remove subtypes with very few cells
data <- data %>% filter(subclass != 'L2/3') %>% filter(subclass != 'Meis2')
data$subclass <- as.factor(data$subclass)

## REGRESSION ----
# subclass
p1 <- ggplot(
  data,
  aes(x = subclass, y = Macrophage, fill = subclass)) + geom_point() +
  ggtitle('Neuronal subclass vs. microglia score\nin the mouse patch-seq dataset') +
  ylab('Microglia score') + xlab ('Neuronal subclass') +
  theme_bw() + geom_boxplot() + scale_fill_brewer(palette="Spectral")

microglianess <- data$Macrophage
subclass <- data$subclass
m1 <- lm(microglianess ~ subclass)
summary(m1) # Multiple R-squared:  0.009114,	Adjusted R-squared:  0.007771 

# cell soma normalized depth
p2 <- ggplot(
  data,
  aes(x = cell_soma_normalized_depth, y = Macrophage)) +
  geom_smooth(method = "lm", se = F, color = 'red4') + 
  geom_point(color = "skyblue4") + 
  ggtitle('Cell soma normalized depth vs. microglia score\nin the mouse patch-seq dataset') + ylim(0, 1) + 
  ylab('Microglia score') + xlab ('Cell soma normalized depth') +
  theme_bw()

cell_soma_norm_depth <- data$cell_soma_normalized_depth
m2 <- lm(microglianess ~ cell_soma_norm_depth)
summary(m2) # Multiple R-squared:  0.0005893,	Adjusted R-squared:  -0.001105

p1 + p2 + plot_layout(ncol = 2, widths=c(4, 3))

# if want to add more colors to graph
p3 <- ggplot(
  data,
  aes(x = cell_soma_normalized_depth, y = Macrophage)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point(aes(color = cell_soma_normalized_depth)) +
  ggtitle('Cell soma normalized depth\nvs. microglia score in the mouse\npatch-seq dataset') + ylim(0, 1) + 
  ylab('Microglia score') + xlab ('Cell soma normalized depth') +
  theme_bw() + scale_color_distiller(palette = "Spectral")

# 6 DUMMY VOLCANO ----
# https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

library(ggrepel)

de <- read.xlsx('/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/figures/files_too_large_for_git/m_QC_de.xlsx', sheetName = 'm_QC_de')

v1 <- ggplot(data = de, aes(x = avg_log2FC, y = -log(p_val))) + geom_point() + theme_minimal()

v2 <- v1 + geom_vline(xintercept = c(-2.5, 2.5), col = "red4") + geom_hline(yintercept = -log10(0.05), col = "red4")

# add a column of NAs
de$diffexpressed <- "NO"
# if avg_log2FC > 2.5 and p_val < 0.05, set as "UP" 
de$diffexpressed[de$avg_log2FC > 2.5 & de$p_val < 0.05] <- "UP"
# if avg_log2FC < -2.5 and p_val < 0.05, set as "DOWN"
de$diffexpressed[de$avg_log2FC < -2.5 & de$p_val < 0.05] <- "DOWN"

v1 <- ggplot(data = de, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed)) + geom_point() + theme_minimal()

v2 <- v1 + geom_vline(xintercept = c(-2.5, 2.5), col = "red4") +
  geom_hline(yintercept = -log10(0.05), col = "red4")

mycolors <- c("skyblue4", "black", "red4")
names(mycolors) <- c("DOWN", "UP", "NO")
v3 <- v2 + scale_colour_manual(values = mycolors)

# label
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$"NA."[de$diffexpressed != "NO"]

v4 <- ggplot(data = de, aes(x = avg_log2FC, y=-log10(p_val), col = diffexpressed, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

# plot adding up all layers we have seen so far
volcano <- ggplot(data = de, aes(x=avg_log2FC, y=-log10(p_val), col = diffexpressed, label = delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("skyblue4", "red4", "black")) +
  geom_vline(xintercept = c(-2.5, 2.5), col = "black") +
  geom_hline(yintercept = -log10(0.05), col="black")
