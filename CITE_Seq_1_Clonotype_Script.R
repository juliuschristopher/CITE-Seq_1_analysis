## ##CITE-Seq 1 Script - Clonotype Analysis####
####Setup####
###Load required packages
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

###Load Seurat object
experiment <- LoadH5Seurat("SeuratProject.h5Seurat")

###Individual plots
DefaultAssay(experiment) <- "RNA"

p1=DimPlot(experiment, label = TRUE,cols=colbig,reduction = "rna.umap", label.size = 2.5) + NoLegend()
p2=DimPlot(experiment, label = TRUE,cols=colbig,reduction = "adt.umap", label.size = 2.5) + NoLegend()
p3=DimPlot(experiment, label = TRUE,cols=colbig, reduction = "wnn.umap", label.size = 2.5) + NoLegend()

p1
p2
p3

DefaultAssay(experiment) <- "RNA"
DefaultAssay(experiment) <- "ADT"

####Simple plotting####
###Umap-wnn by mouse
plot_mouse <- DimPlot(experiment, label = TRUE,reduction = "wnn.umap", label.size = 2.5, group.by = "orig.ident") + ggtitle("Coloured by mouse")

###Umap-wnn by sample
DimPlot(experiment, label = TRUE,cols=colbig, reduction = "wnn.umap", label.size = 2.5, split.by = "orig.ident", ncol = 2) + NoLegend()

###Umap-wnn by cell cycle stage
DimPlot(experiment, label = TRUE,reduction = "wnn.umap", label.size = 2.5, group.by = "Phase") + ggtitle("Coloured by cell cycle stage")

###Feature and violin plot
FeaturePlot(experiment, features = c("CXCR5"), reduction = "wnn.umap")
VlnPlot(experiment, feature = "IgM")

####Clonotype analysis####
###Data visualization
###Percent/total number of unique clonotypes 
quantContig(combined, cloneCall = "gene+nt", scale = T) #percent of unique clonotypes of total size of the size of clonotyeps
quantContig(combined, cloneCall = "gene+nt", scale = F) #number of uniqe clonotypes

quantContig(combined, cloneCall = "gene+nt", scale = T, chain = "IGH") + ggtitle("IGH")#by IGH
quantContig(combined, cloneCall = "gene+nt", scale = F, chain = "IGH") + ggtitle("IGH")#by IGH
quantContig(combined, cloneCall = "gene+nt", scale = T, chain = "IGL") + ggtitle("IGL")#by IGL
quantContig(combined, cloneCall = "gene+nt", scale = F, chain = "IGL") + ggtitle("IGL")#by IGL

###Abundance of clonotypes
Abundance_clonotypes <- abundanceContig(combined, cloneCall = "gene", scale = F, exportTable = T)
Abundance_clonotypes <- Abundance_clonotypes %>%
  arrange(desc(Abundance))
Abundance_clonotypes

abundanceContig(combined, cloneCall = "gene", scale = T)

###Length of clonotypes
lengthContig(combined, cloneCall = "aa")
lengthContig(combined, cloneCall = "nt")

###Compare clonotypes
compareClonotypes(combined, samples = c("a", "b"), cloneCall = "aa", graph = "alluvial") #Computationally intense

###Visualise Gene Usage
vizGenes(combined, gene = "V", chain = "IGH", plot = "bar", order = "variance", scale = TRUE)
vizGenes(combined, gene = "V", chain = "IGL", plot = "bar", order = "variance", scale = TRUE)

vizGenes(combined, gene = "V", chain = "IGL", plot = "heatmap", scale = TRUE, order = "gene")

###Clonal overlap
clonalOverlap(combined, cloneCall = "gene+nt", 
              method = "morisita")

###Clonotype proportion
clonalProportion(combined, cloneCall = "gene")
clonalProportion(combined, cloneCall = "nt")

###Clonal Homeostasis
clonalHomeostasis(combined, cloneCall = "gene")
clonalHomeostasis(combined, cloneCall = "nt")
