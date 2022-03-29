## ##CITE-Seq 1: Experiment.lc.plain batch correction####
setwd("/Volumes/GoogleDrive/Shared drives/Okkengroup/Experiments/Julius/Experiments/CITE-Sequencing/CITE-Seq (1)/CITE-Seq (1) Results/JCB_CITE-Seq_Analysis/Entire_dataset/CITE-Seq_1_code_GIT/Subsetting")

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
library(escape)
library(dittoSeq)
library(SingleCellExperiment)
library(ggsci)
library(pals)
library(harmony)
library(gridExtra)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(org.Mm.eg.db)

#Load Seurat object - experiment.lc.plain
experiment.lc.plain <- LoadH5Seurat("experiment.lc.plain.h5seurat")
head(experiment.lc.plain[[]])
experiment <- experiment.lc.plain

###RNA####
###Normalise subset###
DefaultAssay(experiment) <- "RNA" #For log normalisation and SCTransform

##RNA normalisation
experiment <- NormalizeData(experiment, verbose = TRUE)
experiment <- FindVariableFeatures(experiment, nfeatures = 3000)
experiment <- ScaleData(experiment)

#Or
experiment <-  SCTransform(experiment, verbose = TRUE)
experiment[["SCT"]]

##Visualisation
top20 <-  head(VariableFeatures(experiment), 20)
plot1.1 <-  VariableFeaturePlot(experiment)
top20_plot <-  LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)

##RNA PCA
experiment <- RunPCA(experiment, verbose = FALSE, features = VariableFeatures(object = experiment))
pca_variance <- experiment@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #22

##RNA clustering
DefaultAssay(experiment) <- "RNA" #For log normalisation
DefaultAssay(experiment) <- "SCT" #For SCTransform

experiment <- FindNeighbors(experiment, dims = 1:22)
experiment <- FindClusters(experiment, resolution = 1.5, verbose = FALSE) #1.5 for the resolution
clustree(experiment, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
experiment <- RunUMAP(experiment, dims = 1:22, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
experiment_p1 <- DimPlot(experiment, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("RNA Clustering") + theme_void() + NoLegend()
experiment_p1 <- experiment_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####ADT####
DefaultAssay(experiment) <- "ADT"

##ADT normalisation
VariableFeatures(experiment) <- rownames(experiment[["ADT"]])
experiment <- NormalizeData(experiment, normalization.method = "CLR", margin = 2)
experiment <- ScaleData(experiment)

##ADT PCA
experiment <- RunPCA(experiment, reduction.name = 'apca', approx = FALSE)
apca_variance <- experiment@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #24

##ADT clustering
experiment <- FindNeighbors(experiment, dims = 1:24, reduction = "apca")
experiment <- FindClusters(experiment, resolution = 1.5, verbose = FALSE) #1.5 for the resolution
clustree(experiment, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
experiment <- RunUMAP(experiment, reduction = 'apca', dims = 1:24, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
experiment_p2 <- DimPlot(experiment, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("ADT Clustering") + theme_void() + NoLegend()
experiment_p2 <- experiment_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####WNN####
DefaultAssay(experiment) <- "RNA" #For log normalisation
DefaultAssay(experiment) <- "SCT" #For SCTransform

##Combine into wnn plot
experiment <- FindMultiModalNeighbors(
  experiment, reduction.list = list("pca", "apca"), 
  dims.list = list(1:22, 1:24), modality.weight.name = "RNA.weight")

##WNN clustering
experiment <- FindClusters(experiment, graph.name = "wsnn", algorithm = 3, resolution = 1.5, verbose = TRUE) #1.5 for the resolution
clustree(experiment, prefix = "wsnn_res.") + theme(legend.position="bottom")
experiment <- RunUMAP(experiment, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
experiment_p3 <- DimPlot(experiment, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("WNN Clustering") + theme_void() + NoLegend()
experiment_p3 <- experiment_p3 + theme(plot.title = element_text(color="black", size=25, face="bold"))

##Change Idents
Idents(experiment)
experiment[["old.ident"]] <- Idents(experiment)
Idents(experiment) <- experiment[["orig.ident"]]
experiment <- RenameIdents(experiment, `c1` = "Mb1 Cyp11a1 KO 1", `c2` = "Mb1 Cyp11a1 KO 2", `d1` = "Mb1 E1020K Cyp11a1 KO 1", `d2` = "Mb1 E1020K Cyp11a1 KO 2", `a` = "WT 1", `b` = "WT 2", `f` = "E1020K")
experiment[["orig.ident"]] <- Idents(experiment)
Idents(experiment) <- experiment[["old.ident"]]
Idents(experiment)

####Harmony reductions####
Idents(experiment) <- experiment[["orig.ident"]]
experiment[["Genotype"]] <- Idents(experiment)
experiment <- RenameIdents(experiment, `Mb1 Cyp11a1 KO 1` = "Batch 1", `Mb1 Cyp11a1 KO 2` = "Batch 1", `Mb1 E1020K Cyp11a1 KO 1` = "Batch 1", `Mb1 E1020K Cyp11a1 KO 2` = "Batch 1", `WT 1` = "Batch 2", `WT 2` = "Batch 2", `E1020K` = "Batch 2")
experiment[["Batch"]] <- Idents(experiment)
Idents(experiment) <- experiment[["old.ident"]]

##Batch effect for RNA
DefaultAssay(experiment) <- "RNA"
experiment = RunHarmony(experiment, "Batch", plot_convergence = TRUE, reduction = "pca", assay.use = "RNA", reduction.save = "harmony.rna")
harmony_rna_embeddings = Embeddings(experiment, 'harmony.rna')
experiment = experiment %>%
  RunUMAP(reduction = "harmony.rna", dims = 1:22, assay = 'RNA', reduction.name = 'harmony.rna.umap', reduction.key = 'harmony.rnaUMAP_') %>%
  FindNeighbors(reduction = "harmony.rna", dims = 1:22) %>%
  FindClusters(resolution = 1.5, verbose = FALSE)
experiment_p4 <- DimPlot(experiment, label = TRUE, reduction = "harmony.rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Harmony RNA Clustering") + theme_void() + NoLegend()
experiment_p4 <- experiment_p4 + theme(plot.title = element_text(color="black", size=25, face="bold"))


##Batch effect for ADT
DefaultAssay(experiment) <- "ADT"
experiment = RunHarmony(experiment, "Batch", plot_convergence = TRUE, reduction = "apca", assay.use = "ADT", reduction.save = "harmony.adt")
harmony_adt_embeddings = Embeddings(experiment, 'harmony.adt')
experiment = experiment %>%
  RunUMAP(reduction = "harmony.adt", dims = 1:24, assay = 'ADT', reduction.name = 'harmony.adt.umap', reduction.key = 'harmony.adtUMAP_') %>%
  FindNeighbors(reduction = "harmony.adt", dims = 1:24) %>%
  FindClusters(resolution = 1.5, verbose = FALSE)
experiment_p5 <- DimPlot(experiment, label = TRUE, reduction = "harmony.adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Harmony ADT Clustering") + theme_void() + NoLegend()
experiment_p5 <- experiment_p5 + theme(plot.title = element_text(color="black", size=25, face="bold"))


##WNN clustering with harmony reductions
DefaultAssay(experiment) <- "RNA"
experiment <- FindMultiModalNeighbors(
  experiment, reduction.list = list("harmony.rna", "harmony.adt"), 
  dims.list = list(1:22, 1:24), modality.weight.name = "harmony.weight", weighted.nn.name = "harmony.weighted.nn", snn.graph.name = "harmony.wsnn")
experiment <- FindClusters(experiment, graph.name = "harmony.wsnn", algorithm = 3, resolution = 1.5, verbose = TRUE) #1.5 for the resolution
experiment <- RunUMAP(experiment, nn.name = "harmony.weighted.nn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_")
experiment_p6 <- DimPlot(experiment, label = TRUE, reduction = "harmony.wnn.umap", pt.size = 1.2, label.size = 3, label.box = TRUE) +  ggtitle("Highlighted by cluster") + theme_bw() + NoLegend()
experiment_p6 <- experiment_p6 + theme(plot.title = element_text(color="black", size=15, face="bold")) + xlab("UMAP1") + ylab("UMAP2")

Batch_rna.umap <- DimPlot(experiment, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("RNA UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_adt.umap <- DimPlot(experiment, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("ADT UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_wnn.umap <- DimPlot(experiment, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("WNN UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))

Batch_harmony.rna.umap <- DimPlot(experiment, label = TRUE, reduction = "harmony.rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("Harmony RNA UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_harmony.adt.umap <- DimPlot(experiment, label = TRUE, reduction = "harmony.adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("Harmony ADT UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_harmony.wnn.umap <- DimPlot(experiment, label = TRUE, reduction = "harmony.wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("Harmony WNN UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))

Batch_plot <- Batch_rna.umap + Batch_adt.umap + Batch_wnn.umap + Batch_harmony.rna.umap + Batch_harmony.adt.umap + Batch_harmony.wnn.umap

plot01 <- FeaturePlot(experiment, features = "Cd19", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.3) + theme_bw() + ggtitle("CD19 TotalSeqC-binding") +
  xlab("UMAP1") + ylab("UMAP2") + theme(plot.title = element_text(color="black", size=15, face="bold"))

grid.arrange(experiment_p6, plot01, ncol = 3)

experiment
head(experiment[[]])
SaveH5Seurat(experiment, filename = "experiment.lc.batch", overwrite = TRUE)
