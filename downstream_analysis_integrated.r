# ------ Analysis of scRNA-seq data: GSE229413 -------
#### Author: Aida Ortiz Garcia

## Libraries
library(openxlsx)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(data.table)
library(magrittr)
library(viridis)
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# Data: Integrated data
pancreas.combined <- readRDS("integrated_data_pancreas.rds")

DefaultAssay(pancreas.combined) <- "integrated"

#scaling data
pancreas.combined <- ScaleData(pancreas.combined, verbose = FALSE)

# We can now do PCA, which is a common way of linear dimensionality reduction. By default we use 2000 most variable genes.
pancreas.combined <- RunPCA(pancreas.combined, features = VariableFeatures(object = pancreas.combined))

# checking for the righ number of pcs
ElbowPlot(pancreas.combined, ndims = 50)

pancreas.combined <- RunUMAP(pancreas.combined, reduction = "pca", dims = 1:40)

pancreas.combined <- FindNeighbors(pancreas.combined, reduction = "pca", dims = 1:40)

#(pancreas.combined, "pancreas_integrated2.rds")

pancreas.combined <- FindClusters(pancreas.combined, resolution = 0.1)

options(repr.plot.width=10, repr.plot.height=7)
DimPlot(pancreas.combined, reduction = "umap", label = TRUE, repel = TRUE)
#ggsave("pancreas.combined_UMAP.pdf", height=12, width=18, units="cm")

options(repr.plot.width=15, repr.plot.height=7)
DimPlot(pancreas.combined, reduction = "umap", split.by = "sample")


options(repr.plot.width=10, repr.plot.height=7)
DimPlot(pancreas.combined, reduction = "umap", group.by = "sample")

pancreas.combined

# Data: Integrated data
saveRDS(pancreas.combined, "data/integrated_data_pancreas2.rds")

# Data: Integrated data
pancreas.combined <- readRDS("../data/integrated_data_pancreas2.rds")

markers_all <- FindAllMarkers(pancreas.combined, only.pos=TRUE)

head(markers_all,5)

write.xlsx(markers_all, "allmarkers_clusters.xlsx", col.names = TRUE, row.names = TRUE)

# Data: Integrated data
pancreas.combined <- readRDS("data/integrated_data_pancreas2.rds")

genes <- c("SPINK1", "PRSS3",
              "EPCAM","KRT19",
              "CDH5", "VWF",
              "PDGFRB", "TAGLN",
              "HLA-DRA", "CD14", "FCGR3A",
              "CD79A",'MS4A1',
              "CD3E",
              "MKI67",
              "INS", "CPE",
              "TPSAB1", "CPA3")

options(repr.plot.width=15, repr.plot.height=7)
DotPlot(pancreas.combined, features = genes, scale = TRUE) + scale_y_discrete(position = "left") + RotatedAxis()

options(repr.plot.width=15, repr.plot.height=10)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "SPINK1", pt.size = 1, min.cutoff = 'q10') + NoAxes() + scale_colour_viridis())
suppressMessages(p2<- FeaturePlot(pancreas.combined, "PRSS3", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
CombinePlots(plots = list(p1,p2),legend="bottom")

options(repr.plot.width=15, repr.plot.height=10)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "EPCAM", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
suppressMessages(p2<- FeaturePlot(pancreas.combined, "KRT19", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
CombinePlots(plots = list(p1,p2),legend="bottom")

options(repr.plot.width=15, repr.plot.height=10)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "CDH5", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
suppressMessages(p2<- FeaturePlot(pancreas.combined, "VWF", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
CombinePlots(plots = list(p1,p2),legend="bottom")

options(repr.plot.width=15, repr.plot.height=10)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "PDGFRB", pt.size = 1, min.cutoff = 'q10') + NoAxes() + scale_colour_viridis())
suppressMessages(p2 <- FeaturePlot(pancreas.combined, "ACTA2", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
CombinePlots(plots = list(p1,p2),legend="bottom")

options(repr.plot.width=10, repr.plot.height=7)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "CD3E", pt.size = 1, min.cutoff = 'q10') + NoAxes() + scale_colour_viridis())
CombinePlots(plots = list(p1),legend="bottom")

options(repr.plot.width=15, repr.plot.height=10)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "CD79A", pt.size = 1, min.cutoff = 'q10') + NoAxes() + scale_colour_viridis())
suppressMessages(p2<- FeaturePlot(pancreas.combined, "MS4A1", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
CombinePlots(plots = list(p1,p2),legend="bottom")

options(repr.plot.width=10, repr.plot.height=7)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "MKI67", pt.size = 1, min.cutoff = 'q10') + NoAxes() + scale_colour_viridis())
CombinePlots(plots = list(p1),legend="bottom")

options(repr.plot.width=15, repr.plot.height=10)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "INS", pt.size = 1, min.cutoff = 'q10') + NoAxes() + scale_colour_viridis())
suppressMessages(p2 <- FeaturePlot(pancreas.combined, "CPE", pt.size = 1, min.cutoff = 'q10') + NoAxes() + scale_colour_viridis())
CombinePlots(plots = list(p1,p2),legend="bottom")

options(repr.plot.width=15, repr.plot.height=10)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "HLA-DRA", pt.size = 1, min.cutoff = 'q10') + NoAxes() + scale_colour_viridis())
suppressMessages(p2<- FeaturePlot(pancreas.combined, "CD14", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
suppressMessages(p3<- FeaturePlot(pancreas.combined, "FCGR3A", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
CombinePlots(plots = list(p1,p2,p3),legend="bottom")

options(repr.plot.width=15, repr.plot.height=10)
suppressMessages(p1 <- FeaturePlot(pancreas.combined, "TPSAB1", pt.size = 1, min.cutoff = 'q10') + NoAxes() + scale_colour_viridis())
suppressMessages(p2<- FeaturePlot(pancreas.combined, "CPA3", pt.size = 1, min.cutoff = 'q10') + NoAxes() +  scale_colour_viridis())
CombinePlots(plots = list(p1,p2),legend="bottom")

options(repr.plot.width=10, repr.plot.height=7)
DimPlot(pancreas.combined, reduction = "umap", label = TRUE, repel = TRUE)

known_markers <- read.table("../data/PanglaoDB_markers_27_Mar_2020.tsv.txt", header=T, sep="\t")

head(known_markers,5)

known_markers[known_markers$official.gene.symbol=="HBM",]

# Rename all identities
pancreas.combined_renamed <- RenameIdents(object = pancreas.combined, 
                               "0" = "Macrophagues",
                               "1" = "Acinar cells",
                               "2" = "T cells",
                               "3" = "Acinar cells",
                               "4" = "Ductal cells",
                         "5" = "Fibroblasts",
                         "6" = "Ductal cells",
                              "7" = "Probably immune",
                              "8" = "Mast cells",
                              "9" = "Endothelial cells",
                              "10" = "B cells",
                              "11" = "T cells",
                              "12" = "Macrophagues",
                              "13" = "Endocrine cells",
                              "14" = "CESCs"
                  
                        )

options(repr.plot.width=10, repr.plot.height=7)
DimPlot(pancreas.combined_renamed, reduction = "umap", label = TRUE, repel = TRUE)


options(repr.plot.width=15, repr.plot.height=7)
DimPlot(pancreas.combined_renamed, reduction = "umap", split.by = "sample")


saveRDS(pancreas.combined_renamed, "../data/renamed_integrated_data_pancreas2.rds")

# Data: Integrated data
pancreas.combined_renamed <- readRDS("../data/renamed_integrated_data_pancreas2.rds")

options(repr.plot.width=15, repr.plot.height=7)
DotPlot(pancreas.combined_renamed, features = genes, scale = TRUE) + scale_y_discrete(position = "left") + RotatedAxis()

m <- SplitObject(pancreas.combined_renamed)
markers_list <- lapply(m,function(x) FindMarkers(x,ident.1="Healthy", ident.2="Adj_Normal", group.by="sample", method='MAST', logfc.threshold=0.2,min.diff.pct=0.2,only.pos=TRUE))
#markers_list <- readRDS("markers_list.rds")

head(markers_list$Fibroblasts[order(markers_list$Fibroblasts$avg_log2FC, decreasing = TRUE), ],5)

write.xlsx(markers_list, "DE_genes.xlsx", col.names = TRUE, row.names = TRUE)


