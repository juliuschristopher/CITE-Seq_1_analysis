## ##CITE-Seq 1 Script - Seurat object generation - humanised upper case genes ####
####Setup####
#Load required packages
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
library(Matrix)
library(stringr)


#Alter working capacity
plan()

plan("multiprocess", workers = 4)
options(future.globals.maxSize= 10097152000) # 10Gb


####Load the 10X Cell Ranger output####
#Read the 10x Cell Ranger Output
c1_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Sample_GE_out/C1_GE/outs/filtered_feature_bc_matrix")
c2_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Sample_GE_out/C2_GE/outs/filtered_feature_bc_matrix")
d1_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Sample_GE_out/D1_GE/outs/filtered_feature_bc_matrix")
d2_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Sample_GE_out/D2_GE/outs/filtered_feature_bc_matrix")
a_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/A_WT_GE/outs/filtered_feature_bc_matrix")
b_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/B_WT_GE/outs/filtered_feature_bc_matrix")
f_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/F_E1020K_GE/outs/filtered_feature_bc_matrix")

#Add the sample to the cell names, consistent with antibody data below
colnames(c1_ge.data)=gsub("-1","_c1",colnames(c1_ge.data))
colnames(c2_ge.data)=gsub("-1","_c2",colnames(c2_ge.data))
colnames(d1_ge.data)=gsub("-1","_d1",colnames(d1_ge.data))
colnames(d2_ge.data)=gsub("-1","_d2",colnames(d2_ge.data))
colnames(a_ge.data)=gsub("-1","_a",colnames(a_ge.data))
colnames(b_ge.data)=gsub("-1","_b",colnames(b_ge.data))
colnames(f_ge.data)=gsub("-1","_f",colnames(f_ge.data))

#Capitalise with lower case letters (genes) - turn genes into correct nomenclature
rownames(c1_ge.data)=str_to_title(rownames(c1_ge.data))
rownames(c2_ge.data)=str_to_title(rownames(c2_ge.data))
rownames(d1_ge.data)=str_to_title(rownames(d1_ge.data))
rownames(d2_ge.data)=str_to_title(rownames(d2_ge.data))
rownames(a_ge.data)=str_to_title(rownames(a_ge.data))
rownames(b_ge.data)=str_to_title(rownames(b_ge.data))
rownames(f_ge.data)=str_to_title(rownames(f_ge.data))
head(f_ge.data)


#Basic functions to convert mouse to human gene names (biomaRt)
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(x), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F)
  #Index and replace
  i1 <- with(x, rownames(x) %in% genesV2$MGI.symbol)
  rownames(x)[i1] <- genesV2$HGNC.symbol
  #Uppercase all remaining gene names
  rownames(x)=toupper(rownames(x))
  return(x)
}

#Convert mouse genes to human genenames with convertMouseGeneList function
c1_ge.data <- convertMouseGeneList(c1_ge.data)
c2_ge.data <- convertMouseGeneList(c2_ge.data)
d1_ge.data <- convertMouseGeneList(d1_ge.data)
d2_ge.data <- convertMouseGeneList(d2_ge.data)
a_ge.data <- convertMouseGeneList(a_ge.data)
b_ge.data <- convertMouseGeneList(b_ge.data)
f_ge.data <- convertMouseGeneList(f_ge.data)

####Load 10X Antibody data####
#Read the 10x Antibody output
c1_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Antibody_fraction/C1_SP_out_2/umi_count",gene.column=1)
c2_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Antibody_fraction/C2_SP_out_2/umi_count",gene.column=1)
d1_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Antibody_fraction/D1_SP_out_2/umi_count",gene.column=1)
d2_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Antibody_fraction/D2_SP_out_2/umi_count",gene.column=1)
a_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/A_WT/umi_count",gene.column=1)
b_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/B_WT/umi_count",gene.column=1)
f_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/F_E1020K/umi_count",gene.column=1)

#Tidy up the rownames from the data
rownames(c1_ab.data)=gsub("-[^-]+$","",rownames(c1_ab.data),perl=TRUE)
rownames(c2_ab.data)=gsub("-[^-]+$","",rownames(c2_ab.data),perl=TRUE)
rownames(d1_ab.data)=gsub("-[^-]+$","",rownames(d1_ab.data),perl=TRUE)
rownames(d2_ab.data)=gsub("-[^-]+$","",rownames(d2_ab.data),perl=TRUE)
rownames(a_ab.data)=gsub("-[^-]+$","",rownames(a_ab.data),perl=TRUE)
rownames(b_ab.data)=gsub("-[^-]+$","",rownames(b_ab.data),perl=TRUE)
rownames(f_ab.data)=gsub("-[^-]+$","",rownames(f_ab.data),perl=TRUE)

#Add the Sample to the cell names in each sample
colnames(c1_ab.data)=paste(colnames(c1_ab.data),"_c1",sep="")
colnames(c2_ab.data)=paste(colnames(c2_ab.data),"_c2",sep="")
colnames(d1_ab.data)=paste(colnames(d1_ab.data),"_d1",sep="")
colnames(d2_ab.data)=paste(colnames(d2_ab.data),"_d2",sep="")
colnames(a_ab.data)=paste(colnames(a_ab.data),"_a",sep="")
colnames(b_ab.data)=paste(colnames(b_ab.data),"_b",sep="")
colnames(f_ab.data)=paste(colnames(f_ab.data),"_f",sep="")
head(f_ab.data)

####Combine 10X Cell Ranger and Antibody Data into a Suerat Object####

m <- Matrix(nrow = nrow(c1_ab.data), ncol = ncol(c1_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(c1_ab.data)
colnames(m)=colnames(c1_ge.data)
common=intersect(colnames(c1_ge.data),colnames(c1_ab.data))
m[,common]=c1_ab.data[,common]
c1 = CreateSeuratObject(counts = c1_ge.data,project="c1", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
c1[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(c2_ab.data), ncol = ncol(c2_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(c2_ab.data)
colnames(m)=colnames(c2_ge.data)
common=intersect(colnames(c2_ge.data),colnames(c2_ab.data))
m[,common]=c2_ab.data[,common]
c2 = CreateSeuratObject(counts = c2_ge.data,project="c2", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
c2[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(d1_ab.data), ncol = ncol(d1_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(d1_ab.data)
colnames(m)=colnames(d1_ge.data)
common=intersect(colnames(d1_ge.data),colnames(d1_ab.data))
m[,common]=d1_ab.data[,common]
d1 = CreateSeuratObject(counts = d1_ge.data,project="d1", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
d1[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(d2_ab.data), ncol = ncol(d2_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(d2_ab.data)
colnames(m)=colnames(d2_ge.data)
common=intersect(colnames(d2_ge.data),colnames(d2_ab.data))
m[,common]=d2_ab.data[,common]
d2 = CreateSeuratObject(counts = d2_ge.data,project="d2", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
d2[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(a_ab.data), ncol = ncol(a_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(a_ab.data)
colnames(m)=colnames(a_ge.data)
common=intersect(colnames(a_ge.data),colnames(a_ab.data))
m[,common]=a_ab.data[,common]
a = CreateSeuratObject(counts = a_ge.data,project="a", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
a[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(b_ab.data), ncol = ncol(b_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(b_ab.data)
colnames(m)=colnames(b_ge.data)
common=intersect(colnames(b_ge.data),colnames(b_ab.data))
m[,common]=b_ab.data[,common]
b = CreateSeuratObject(counts = b_ge.data,project="b", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
b[["ADT"]] <- adt_assay

m <- Matrix(nrow = nrow(f_ab.data), ncol = ncol(f_ge.data), data = 0, sparse = TRUE)
rownames(m)=rownames(f_ab.data)
colnames(m)=colnames(f_ge.data)
common=intersect(colnames(f_ge.data),colnames(f_ab.data))
m[,common]=f_ab.data[,common]
f = CreateSeuratObject(counts = f_ge.data,project="f", min.cells = 3)
adt_assay <- CreateAssayObject(counts = m)
f[["ADT"]] <- adt_assay
head(f[[]])

####Incoperate VDJ data####
#Load contig file
c1_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/C1_VDJ/outs/filtered_contig_annotations.csv")
c2_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/C2_VDJ/outs/filtered_contig_annotations.csv")
d1_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/D1_VDJ/outs/filtered_contig_annotations.csv")
d2_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/D2_VDJ/outs/filtered_contig_annotations.csv")
a_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/A_WT_VDJ/outs/filtered_contig_annotations.csv")
b_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/B_WT_VDJ/outs/filtered_contig_annotations.csv")
f_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/F_E1020K_VDJ/outs/filtered_contig_annotations.csv")

#match barcode names with GE and ADT data
c1_cl.data$barcode=gsub("-1","_c1",c1_cl.data$barcode)
c2_cl.data$barcode=gsub("-1","_c2",c2_cl.data$barcode)
d1_cl.data$barcode=gsub("-1","_d1",d1_cl.data$barcode)
d2_cl.data$barcode=gsub("-1","_d2",d2_cl.data$barcode)
a_cl.data$barcode=gsub("-1","_a",a_cl.data$barcode)
b_cl.data$barcode=gsub("-1","_b",b_cl.data$barcode)
f_cl.data$barcode=gsub("-1","_f",f_cl.data$barcode)

contig_list <- list(c1_cl.data, c2_cl.data, d1_cl.data, d2_cl.data, a_cl.data, b_cl.data, f_cl.data)
head(contig_list[[1]])

#Generate combined object
combined <- combineBCR(contig_list, samples = c("c1", "c2", "d1", "d2", "a", "b", "f"))
combined[[7]]

str(combined)
head(combined[[7]])

#Make sure barcodes are identical to GE and ADT data
combined$c1$barcode=gsub("c1_","",combined$c1$barcode)
combined$c2$barcode=gsub("c2_","",combined$c2$barcode)
combined$d1$barcode=gsub("d1_","",combined$d1$barcode)
combined$d2$barcode=gsub("d2_","",combined$d2$barcode)
combined$a$barcode=gsub("a_","",combined$a$barcode)
combined$b$barcode=gsub("b_","",combined$b$barcode)
combined$f$barcode=gsub("f_","",combined$f$barcode)
head(combined$f)

#####Process samples as one####
experiments=c(c1,c2,d1,d2,a,b,f)
experiment_names=c("c1","c2","d1","d2","a","b","f")

experiment<-merge(x= c1, y=c(c2,d1,d2,a,b,f))
head(experiment[[]])

####Merge Seurat object with VDJ data####
experiment <- combineExpression(combined, experiment, cloneCall="gene", group.by =  "sample")

head(experiment[[]])

####Quality control, filtering and normalisation####
#Mitochondrial QC metrics
experiment[["percent.Mt"]] <- PercentageFeatureSet(experiment, pattern = "^MT-")

#Remove where nCount_ADT = 0
DefaultAssay(experiment) <- "ADT"
experiment <- subset(experiment, nCount_ADT >0)
DefaultAssay(experiment) <- "RNA"

#Visualize QC metrics as violin plot
RNA_QC <- VlnPlot(experiment, features = c("nFeature_RNA", "nCount_RNA", "percent.Mt"))
ADT_QC_1 <- VlnPlot(experiment, features = "nFeature_ADT")
ADT_QC_2 <- VlnPlot(experiment, features = "nCount_ADT", y.max = 10000)

RNA_QC
ADT_QC_1
ADT_QC_2

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
feature_count.RNA_vs_percent.Mt = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "percent.Mt") + NoLegend()  +
  ylab("% of mitochondrial genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 5)

feature_count.RNA_vs_feature.RNA = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() +
  ylab("Number of genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 200)

feature_count.RNA_vs_percent.Mt + feature_count.RNA_vs_feature.RNA

feature_count.ADT_vs_feature.ADT = FeatureScatter(experiment, feature1 = "nCount_ADT", feature2 = "nFeature_ADT") + NoLegend() +
  ylab("Number of antibodies") +
  xlab("UMI counts")

feature_count.ADT_vs_feature.ADT

#Generally aim to filter out unique feature counts over 2,500 and less than 200; and percent.mt over 5%
filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes and mitochondrial percentage.")
  seurat_object = subset(seurat_object, subset = nFeature_RNA > 200  & percent.Mt < 5 & nFeature_RNA < 2500)
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  
  
  return(seurat_object)
}

experiment <-  filter_seurat(experiment)

####Normalise dataset - SCTransform####
experiment = SCTransform(experiment, verbose = TRUE)
experiment[["SCT"]]

top20 = head(VariableFeatures(experiment), 20)

plot1.1 = VariableFeaturePlot(experiment)
top20_plot = LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
top20_plot

#Cell Cycle genes
S.genes <- cc.genes.updated.2019$s.genes
S.genes <- lapply(S.genes, toupper)
G2M.genes <-  cc.genes.updated.2019$g2m.genes
G2M.genes <- lapply(G2M.genes,toupper)

experiment <- CellCycleScoring(experiment, s.features=S.genes, g2m.features=G2M.genes, set.ident = TRUE)
Idents(object = experiment) <- "old.ident"

####Dimensionality reduction - PCA####
#Perform linear dimensional reduction (PCA)
experiment <- RunPCA(experiment, verbose = FALSE, features = VariableFeatures(object = experiment))

##Visualise PCA results
print(experiment[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(experiment, dims = 1:4, reduction = "pca", nfeatures = 15)

plot5 <- DimPlot(experiment, reduction = "pca", dims = c(1,2))
plot6 <- DimPlot(experiment, reduction = "pca", dims = c(1,3))
plot5 + plot6

DimHeatmap(experiment, dims = 1:6, cells = 500, balanced = TRUE)

####Determine dimensionality of the dataset - How many principal components should be included to capture the majority of variance?####
pca_variance <- experiment@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #25

#Cluster the cells
experiment <- FindNeighbors(experiment, dims = 1:25)
experiment <- FindClusters(experiment, resolution = 1.3, verbose = FALSE) #1.3 or 1.2 for the resolution

#Cluster Tree Analysis
clustree(experiment, prefix = "SCT_snn_res.") + theme(legend.position="bottom")

#RNA UMAP
experiment <- RunUMAP(experiment, dims = 1:25)
DimPlot(experiment, label = TRUE, cols=colbig) +  ggtitle("RNA Clustering")

####Scale antibody data####
DefaultAssay(experiment) <- "ADT"
VariableFeatures(experiment) <- rownames(experiment[["ADT"]])
experiment <- NormalizeData(experiment, normalization.method = "CLR", margin = 2)
experiment <- ScaleData(experiment)
experiment <- RunPCA(experiment,reduction.name = 'apca')

#Visualise antibody PCA
print(experiment[["apca"]], dims = 1:10, nfeatures = 5)
Plot_13 <- VizDimLoadings(experiment, dims = 1:4, reduction = "apca", nfeatures = 15)
Plot_13

plot14 <- DimPlot(experiment, reduction = "apca", dims = c(1,2), group.by = "orig.ident") + ggtitle("ADT PCA")
plot15 <- DimPlot(experiment, reduction = "apca", dims = c(1,3), group.by = "orig.ident") + ggtitle("ADT PCA")
plot14 + plot15

DimHeatmap(experiment, dims = 1:6, cells = 500, balanced = TRUE, reduction = "apca")

#Determine number of PCs for ADT assay
apca_variance <- experiment@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #25

#Number of clusters for UMAP?

####Combine into wnn plot####
experiment <- FindMultiModalNeighbors(
  experiment, reduction.list = list("pca", "apca"), 
  dims.list = list(1:25, 1:25), modality.weight.name = "RNA.weight")

#UMAP plots for RNA, ADT and WNN
experiment <- RunUMAP(experiment, reduction = 'pca', dims = 1:25, assay = 'RNA', 
                      reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
experiment<- RunUMAP(experiment, reduction = 'apca', dims = 1:25, assay = 'ADT', 
                     reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
experiment <- RunUMAP(experiment, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
experiment <- FindClusters(experiment, graph.name = "wsnn", algorithm = 3, resolution = 1.0, verbose = TRUE)

#Cluster Tree Analysis of wsnn graph
clustree(experiment, prefix = "wsnn_res.") + theme(legend.position="bottom")#1.7

DefaultAssay(experiment) <- "RNA"

p1=DimPlot(experiment, label = TRUE,cols=colbig,reduction = "rna.umap", label.size = 2.5) + NoLegend()
p2=DimPlot(experiment, label = TRUE,cols=colbig,reduction = "adt.umap", label.size = 2.5) + NoLegend()
p3.uc =DimPlot(experiment.uc, label = TRUE,cols=colbig, reduction = "wnn.umap", label.size = 2.5) + NoLegend()
p3.lc =DimPlot(experiment.lc, label = TRUE,cols=colbig, reduction = "wnn.umap", label.size = 2.5) + NoLegend()


p1
p2
p3.uc
p3.lc

experiment.lc <- experiment
head(experiment[[]])
SaveH5Seurat(experiment.lc, filename = "experiment.h.uc", overwrite = TRUE)

FeaturePlot(experiment.lc, features = "Cd19", reduction = "wnn.umap")



