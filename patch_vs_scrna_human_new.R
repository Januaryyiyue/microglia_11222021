# SETUP ----

library(Seurat)
library(dplyr)
library(xlsx)
library(ggplot2)
library(cowplot)
library(ggrepel)
theme_set(theme_cowplot())


# READ PATCH + scRNA SEURAT OBJECT ----
## open and find markers ----
microglia.combined <- readRDS('patchseq_scRNA_human.rds')

DefaultAssay(microglia.combined) <- "integrated"
Idents(microglia.combined) <- "source"

patch_vs_scrna <- FindMarkers(microglia.combined, ident.1 = "patchseq", ident.2 = "scRNA", verbose = FALSE)
write.xlsx(patch_vs_scrna, file = "patch_vs_scrna_de_new.xlsx", sheetName = "patch_vs_scrna_de_new")


## visualize ----
DefaultAssay(microglia.combined) <- "RNA"
FeaturePlot(microglia.combined, features = c("TLR4", "CTSS", "ICAM1"), split.by = "source", max.cutoff = 3, cols = c("grey", "red"))

DimPlot(microglia.combined, reduction = "umap", group.by = "source")
DimPlot(microglia.combined, reduction = "umap", split.by = "source")

avg.microglia.combined <- as.data.frame(log1p(AverageExpression(microglia.combined, verbose = FALSE)$RNA))
avg.microglia.combined$gene <- rownames(avg.microglia.combined)

ggplot(avg.microglia.combined, aes(patchseq, scRNA)) + geom_point() + ggtitle("patch vs scRNA microglia") + geom_text(aes(label = ifelse(patchseq/scRNA > 2, as.character(gene),'')),hjust=0,vjust=0)


# HYPERGEOMETRIC TEST ----
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

# k (filter fdr < 0.05 and logFC > 0) got 253
patch_vs_scrna <- read.xlsx('patch_vs_scrna_de_new.xlsx', sheetName = 'patch_vs_scrna_de_new')
patch_vs_scrna$fdr <- p.adjust(patch_vs_scrna$p_val, "fdr")
patch_vs_scrna <- filter(patch_vs_scrna, fdr < 0.05)
patch_vs_scrna <- filter(patch_vs_scrna, avg_log2FC > 0)
patch_vs_scrna <- dplyr::rename(patch_vs_scrna, 'gene' = 'NA.')
write.xlsx(patch_vs_scrna, file = "patch_vs_scrna_de_new_FDR.xlsx", sheetName = "patch_vs_scrna_de_new_FDR", append = FALSE)

# q
x1 <- dplyr::inner_join(m1, patch_vs_scrna, by = 'gene', copy = FALSE) # 0
x2 <- dplyr::inner_join(m2, patch_vs_scrna, by = 'gene', copy = FALSE) # 0
x3 <- dplyr::inner_join(m3, patch_vs_scrna, by = 'gene', copy = FALSE) # 3
x4 <- dplyr::inner_join(m4, patch_vs_scrna, by = 'gene', copy = FALSE) # 3
x5 <- dplyr::inner_join(m5, patch_vs_scrna, by = 'gene', copy = FALSE) # 7
x6 <- dplyr::inner_join(m6, patch_vs_scrna, by = 'gene', copy = FALSE) # 2
x7 <- dplyr::inner_join(m7, patch_vs_scrna, by = 'gene', copy = FALSE) # 0
x8 <- dplyr::inner_join(m8, patch_vs_scrna, by = 'gene', copy = FALSE) # 2
x9 <- dplyr::inner_join(m9, patch_vs_scrna, by = 'gene', copy = FALSE) # 6

# hypergeometric test 0 vs 1
phyper(0, 57, 21698-57, 253, lower.tail = FALSE, log.p = FALSE) # 0.4879832
phyper(0, 514, 21698-514, 253, lower.tail = FALSE, log.p = FALSE) # 0.9977601
phyper(3, 208, 21698-208, 253, lower.tail = FALSE, log.p = FALSE) # 0.2252669
phyper(3, 454, 21698-454, 253, lower.tail = FALSE, log.p = FALSE) # 0.7784734
phyper(7, 575, 21698-575, 253, lower.tail = FALSE, log.p = FALSE) # 0.3564681
phyper(2, 618, 21698-618, 253, lower.tail = FALSE, log.p = FALSE) # 0.9766212
phyper(0, 318, 21698-318, 253, lower.tail = FALSE, log.p = FALSE) # 0.9766548
phyper(2, 444, 21698-444, 253, lower.tail = FALSE, log.p = FALSE) # 0.8933493
phyper(6, 1312, 21698-1312, 253, lower.tail = FALSE, log.p = FALSE) # 0.9949639
