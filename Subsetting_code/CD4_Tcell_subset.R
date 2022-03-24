## ##Subsetting: CD4+ T cells####
####Set working directory####
setwd("/Volumes/GoogleDrive/Shared drives/Okkengroup/Experiments/Julius/Experiments/CITE-Sequencing/CITE-Seq (1)/CITE-Seq (1) Results/JCB_CITE-Seq_Analysis/Entire_dataset/CITE-Seq_1_code_GIT/Subsetting")

####Load requried packages####
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Matrix)
library(RColorBrewer)
library(writexl)
library(ggridges)
library(clustree)
library(scRepertoire)
library(future)
library(alakazam)
library(immunarch)
library(airr)
library(biomaRt)
library(SeuratDisk)
library(SeuratData)
library(stringr)
library(viridis)


####Define colour palettes####
col_con = viridis(50)
col = scale_colour_brewer(palette = "Spectral")

####Load Seurat object (lower case gene names - mouse)####
experiment.lc <- LoadH5Seurat("experiment.lc.H5Seurat")
head(experiment.lc[[]])

DimPlot(experiment.lc, reduction = "wnn.umap")
FeaturePlot(experiment.lc, feature = "Cd4", reduction = "wnn.umap")

####Subset by ADT feature expression of CD4####
DefaultAssay(experiment.lc) <- "ADT"
Cd4_Tcells <- subset(experiment.lc, subset = Cd4 > 3)
DimPlot(Cd4_Tcells, reduction = "wnn.umap")
DefaultAssay(Cd4_Tcells) <- "SCT"
Cd4_Tcells
head(Cd4_Tcells[[]])


####Normalise and recluster subset####
###RNA
DefaultAssay(experiment) <- "SCT"

##RNA normalisation
Cd4_Tcells = SCTransform(Cd4_Tcells, verbose = TRUE)

##RNA PCA
Cd4_Tcells <- RunPCA(Cd4_Tcells, verbose = FALSE, features = VariableFeatures(object = Cd4_Tcells))
pca_variance <- Cd4_Tcells@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #36

#RNA clustering
DefaultAssay(Cd4_Tcells) <- "SCT"
Cd4_Tcells <- FindNeighbors(Cd4_Tcells, dims = 1:36)
Cd4_Tcells <- FindClusters(Cd4_Tcells, resolution = 0.1, verbose = FALSE) #0.1 for the resolution
clustree(Cd4_Tcells, prefix = "SCT_snn_res.") + theme(legend.position="bottom")
Cd4_Tcells <-RunUMAP(Cd4_Tcells, dims = 1:36, assay = 'SCT', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
p1 <- DimPlot(Cd4_Tcells, label = TRUE, reduction = "rna.umap") +  ggtitle("RNA Clustering")

###ADT
DefaultAssay(Cd4_Tcells) <- "ADT"

#ADT normalisation
DefaultAssay(Cd4_Tcells) <- "ADT"
VariableFeatures(Cd4_Tcells) <- rownames(Cd4_Tcells[["ADT"]])
Cd4_Tcells <- NormalizeData(Cd4_Tcells, normalization.method = "CLR", margin = 2)
Cd4_Tcells <- ScaleData(Cd4_Tcells)
Cd4_Tcells <- RunPCA(Cd4_Tcells,reduction.name = 'apca')

#ADT PCA
apca_variance <- Cd4_Tcells@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #26

#ADT clustering
Cd4_Tcells <- FindNeighbors(Cd4_Tcells, dims = 1:2, reduction = "apca")
Cd4_Tcells <- FindClusters(Cd4_Tcells, resolution = 0.2, verbose = FALSE)#0.2 for the resolution
clustree(Cd4_Tcells, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
Cd4_Tcells <- RunUMAP(Cd4_Tcells, reduction = 'apca', dims = 1:26, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
p2 <- DimPlot(Cd4_Tcells, label = TRUE, reduction = "adt.umap") +  ggtitle("ADT Clustering")

###WNN
#Combine into wnn plot####
Cd4_Tcells <- FindMultiModalNeighbors(
  Cd4_Tcells, reduction.list = list("pca", "apca"), 
  dims.list = list(1:36, 1:26), modality.weight.name = "RNA.weight")

#WNN clustering
Cd4_Tcells <- RunUMAP(Cd4_Tcells, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Cd4_Tcells <- FindClusters(Cd4_Tcells, graph.name = "wsnn", algorithm = 3, resolution = 0.1, verbose = TRUE)#0.5
clustree(Cd4_Tcells, prefix = "wsnn_res.") + theme(legend.position="bottom")
p3 <- DimPlot(Cd4_Tcells, label = TRUE, reduction = "wnn.umap", label.size = 2.5)

head(Cd4_Tcells[[]])
####Subset analysis####
DefaultAssay(Cd4_Tcells) <- "ADT"
DefaultAssay(Cd4_Tcells) <- "SCT"

#By genotype
DimPlot(Cd4_Tcells, label = FALSE,reduction = "adt.umap", label.size = 2.5, group.by = "orig.ident")
DimPlot(Cd4_Tcells, label = FALSE,reduction = "adt.umap", label.size = 2.5, group.by = "Phase")


#By cell cycle genes


#Cluster identification
Idents(object = Cd4_Tcells) <- Cd4_Tcells$wsnn_res.0.1
Cd4_Tcells$seurat_clusters <- Cd4_Tcells$wsnn_res.0.1

FeaturePlot(Cd4_Tcells, features = "Cd38", reduction = "adt.umap")
Cd62L.adt <- FeaturePlot(Cd4_Tcells, features = "Cd62l", reduction = "adt.umap", cols = col_con, pt.size = 0.5) + theme_void()
FeaturePlot(Cd4_Tcells, features = "Sell", reduction = "adt.umap", cols = col_con, pt.size = 0.5) + theme_void()
Cd44.adt <- FeaturePlot(Cd4_Tcells, features = "Cd44", reduction = "adt.umap", cols = col_con, pt.size = 0.5) + theme_void()
Foxp3 <- FeaturePlot(Cd4_Tcells, features = "Foxp3", reduction = "adt.umap", cols = col_con, pt.size = 0.5) + theme_void()
Bcl6 <- FeaturePlot(Cd4_Tcells, features = "Bcl6", reduction = "adt.umap", cols = col_con, pt.size = 0.5) + theme_void()
FeaturePlot(Cd4_Tcells, features = "Tox", reduction = "adt.umap")
FeaturePlot(Cd4_Tcells, features = "Tcf7", reduction = "adt.umap", cols = col_con, pt.size = 0.5) + theme_void()
FeaturePlot(Cd4_Tcells, features = "Plac8", reduction = "adt.umap")
Tbx21 <- FeaturePlot(Cd4_Tcells, features = "Tbx21", reduction = "adt.umap", cols = col_con, pt.size = 0.5) + theme_void()
Gata3 <- FeaturePlot(Cd4_Tcells, features = "Gata3", reduction = "adt.umap", cols = col_con, pt.size = 0.5) + theme_void()



Cd62L.adt + Cd44.adt
Foxp3 + Bcl6
Tbx21 + Gata3



c3_RNA <- FindMarkers(Cd4_Tcells, ident.1 = 6, assay = "SCT")
c3_ADT <- FindMarkers(Cd4_Tcells, ident.1 = 2, assay = "ADT")
?FindMarkers


