## ##CITE-Seq 1 Script - Seurat object generation - capitalised lower case genes ####
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
library(stringr)
library(harmony)
library(viridis)

#Alter working capacity
plan()

plan("multiprocess", workers = 4)
options(future.globals.maxSize= 10097152000) # 10Gb

#Initial setup of colour palettes
col = colorRampPalette(brewer.pal(12, 'Set3'))(20)

colbig = colorRampPalette(brewer.pal(12, 'Set3'))(50)
col_con = viridis(50)
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
?alphabet

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

#Capitalise with lower case letters (genes)
rownames(c1_ge.data)=str_to_title(rownames(c1_ge.data))
rownames(c2_ge.data)=str_to_title(rownames(c2_ge.data))
rownames(d1_ge.data)=str_to_title(rownames(d1_ge.data))
rownames(d2_ge.data)=str_to_title(rownames(d2_ge.data))
rownames(a_ge.data)=str_to_title(rownames(a_ge.data))
rownames(b_ge.data)=str_to_title(rownames(b_ge.data))
rownames(f_ge.data)=str_to_title(rownames(f_ge.data))
head(f_ge.data)

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

#Capitalise with lower case letters (antobodies)
rownames(c1_ab.data)=str_to_title(rownames(c1_ab.data))
rownames(c2_ab.data)=str_to_title(rownames(c2_ab.data))
rownames(d1_ab.data)=str_to_title(rownames(d1_ab.data))
rownames(d2_ab.data)=str_to_title(rownames(d2_ab.data))
rownames(a_ab.data)=str_to_title(rownames(a_ab.data))
rownames(b_ab.data)=str_to_title(rownames(b_ab.data))
rownames(f_ab.data)=str_to_title(rownames(f_ab.data))
head(f_ab.data)

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

####Merge suerat objects####
experiments=c(c1,c2,d1,d2,a,b,f)
experiment_names=c("c1","c2","d1","d2","a","b","f")
experiment<-merge(x= c1, y=c(c2,d1,d2,a,b,f))
head(experiment[[]])

#Merge Seurat object with VDJ data
experiment <- combineExpression(combined, experiment, cloneCall="gene", group.by =  "sample")
head(experiment[[]])

####Quality control, filtering and normalisation####
#Mitochondrial QC metrics
experiment[["percent.Mt"]] <- PercentageFeatureSet(experiment, pattern = "^Mt-")

#Remove where nCount_ADT = 0
DefaultAssay(experiment) <- "ADT"
experiment <- subset(experiment, nCount_ADT > 0)
DefaultAssay(experiment) <- "RNA"

#Visualize QC metrics as violin plot
RNA_QC <- VlnPlot(experiment, features = c("nFeature_RNA", "nCount_RNA", "percent.Mt"))
ADT_QC_1 <- VlnPlot(experiment, features = "nFeature_ADT")
ADT_QC_2 <- VlnPlot(experiment, features = "nCount_ADT", y.max = 10000)

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
feature_count.RNA_vs_percent.Mt = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "percent.Mt") + NoLegend()  +
  ylab("% of mitochondrial genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 5) 
feature_count.RNA_vs_feature.RNA = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() +
  ylab("Number of genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 200) 
feature_count.ADT_vs_feature.ADT = FeatureScatter(experiment, feature1 = "nCount_ADT", feature2 = "nFeature_ADT", jitter = TRUE) + NoLegend() +
  ylab("Number of different antibodies") + ylim(0, 30) +
  xlab("ADT counts") + xlim(0, 40000)

#Generally aim to filter out unique feature counts over 2,500 and less than 200; and percent.mt over 5%
filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes and mitochondrial percentage.")
  seurat_object = subset(seurat_object, subset = nFeature_RNA > 200  & percent.Mt < 5 & nFeature_RNA < 2500)
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  
  
  return(seurat_object)
}

head(experiment[[]])
experiment <-  filter_seurat(experiment)
SaveH5Seurat(experiment, filename = "experiment.lc.plain", overwrite = TRUE)

####Batch-correction####
#Assign merged Seurat object
seurat_object_merged = experiment

#Normalise and calculate PCAs
seurat_object_merged = seurat_object_merged %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = seurat_object_merged@var.genes, npcs = 20, verbose = FALSE)

#Add batch information to orig.ident column in metadata
Idents(seurat_object_merged)
seurat_object_merged[["old.ident"]] <- Idents(seurat_object_merged)
Idents(seurat_object_merged)
seurat_object_merged <- RenameIdents(seurat_object_merged, `c1` = "Batch 1", `c2` = "Batch 1", `d1` = "Batch 1", `d2` = "Batch 1", `a` = "Batch 2", `b` = "Batch 2", `f` = "Batch 2")
seurat_object_merged[["orig.ident"]] <- Idents(seurat_object_merged)
Idents(seurat_object_merged) <- seurat_object_merged[["old.ident"]]
Idents(seurat_object_merged)

#Run harmony
seurat_object_merged = RunHarmony(seurat_object_merged, 
                                  "orig.ident", 
                                  plot_convergence = TRUE)

#Accessing harmony embeddings:
harmony_embeddings = Embeddings(seurat_object_merged, 'harmony')
seurat_object_merged = seurat_object_merged %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)

#Visualise batches
head(seurat_object_merged[[]])
DimPlot(seurat_object_merged, reduction = "umap", group.by = "orig.ident")
FeaturePlot(seurat_object_merged, features = "Sox4")

####RNA####
#Nornmalise dataset - log normalisation#
experiment <- NormalizeData(experiment, verbose = TRUE)
experiment <- FindVariableFeatures(experiment, selection.method = "vst", nfeatures = 3000)
experiment <- ScaleData(experiment)

#Normalise dataset - SCTransform
experiment = SCTransform(experiment, verbose = TRUE)
experiment[["SCT"]]

#Visualise top variable features
top20 = head(VariableFeatures(experiment), 20)
plot1.1 = VariableFeaturePlot(experiment)
top20_plot = LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)

#Cell Cycle genes
S.genes <- cc.genes.updated.2019$s.genes
S.genes <- lapply(S.genes, str_to_title)
G2M.genes <-  cc.genes.updated.2019$g2m.genes
G2M.genes <- lapply(G2M.genes, str_to_title)
experiment <- CellCycleScoring(experiment, s.features=S.genes, g2m.features=G2M.genes, set.ident = TRUE)
Idents(object = experiment) <- "old.ident"

#RNA - PCA
experiment <- RunPCA(experiment, verbose = FALSE, features = VariableFeatures(object = experiment))
pca_variance <- experiment@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #25

#Visualise PCA results - RNA
print(experiment[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(experiment, dims = 1:4, reduction = "pca", nfeatures = 15)
PCA_RNA_1 <- DimPlot(experiment, reduction = "pca", dims = c(1,2))
PCA_RNA_2 <- DimPlot(experiment, reduction = "pca", dims = c(1,3))
DimHeatmap(experiment, dims = 1:6, cells = 500, balanced = TRUE)

#RNA clustering
DefaultAssay(experiment) <- "RNA" #For log normalisation
DefaultAssay(experiment) <- "SCT" #For SCTransform

experiment <- FindNeighbors(experiment, dims = 1:25)
experiment <- FindClusters(experiment, resolution = 1.5, verbose = FALSE)#1.5 for the resolution
clustree(experiment, prefix = "SCT_snn_res.") + theme(legend.position="bottom")
experiment <-RunUMAP(experiment, dims = 1:25, assay = 'SCT', reduction.name = 'SCT.umap', reduction.key = 'SCTUMAP_')
experiment_p1 <- DimPlot(experiment, label = TRUE, reduction = "SCT.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("SCT Clustering") + theme_void() + NoLegend()
experiment_p1 <- experiment_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####ADT####
DefaultAssay(experiment) <- "ADT"

#ADT normalisation
VariableFeatures(experiment) <- rownames(experiment[["ADT"]])
experiment <- NormalizeData(experiment, normalization.method = "CLR", margin = 2)
experiment <- ScaleData(experiment)
experiment <- RunPCA(experiment, reduction.name = 'apca')

#ADT PCA
apca_variance <- experiment@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #25

#Visualise PCA results - ADT
print(experiment[["apca"]], dims = 1:10, nfeatures = 5)
PCA_ADT_1 <- VizDimLoadings(experiment, dims = 1:4, reduction = "apca", nfeatures = 15)
PCA_ADT_2 <- DimPlot(experiment, reduction = "apca", dims = c(1,2), group.by = "orig.ident") + ggtitle("ADT PCA")
PCA_ADT_3 <- DimPlot(experiment, reduction = "apca", dims = c(1,3), group.by = "orig.ident") + ggtitle("ADT PCA")
DimHeatmap(experiment, dims = 1:6, cells = 500, balanced = TRUE, reduction = "apca")

#ADT clustering
experiment <- FindNeighbors(experiment, dims = 1:25, reduction = "apca")
experiment <- FindClusters(experiment, resolution = 1.0, verbose = FALSE) #1 for the resolution
clustree(experiment, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
experiment <- RunUMAP(experiment, reduction = 'apca', dims = 1:25, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
experiment_p2 <- DimPlot(experiment, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("ADT Clustering") + theme_void() + NoLegend()
experiment_p2 <- experiment_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####WNN####
DefaultAssay(experiment) <- "RNA" #For log normalisation
DefaultAssay(experiment) <- "SCT" #For SCTransform

#Combine into wnn plot
experiment <- FindMultiModalNeighbors(
  experiment, reduction.list = list("pca", "apca"), 
  dims.list = list(1:25, 1:25), modality.weight.name = "RNA.weight")

#WNN clustering
experiment <- FindClusters(experiment, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = TRUE) #0.5 for the resolution
clustree(experiment, prefix = "wsnn_res.") + theme(legend.position="bottom")
experiment <- RunUMAP(experiment, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
experiment_p3 <- DimPlot(experiment, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("WNN Clustering") + theme_void() + NoLegend()
experiment_p3 <- experiment_p3 + theme(plot.title = element_text(color="black", size=25, face="bold"))


head(experiment[[]])
experiment.lc <- experiment
head(experiment.lc[[]])
FeaturePlot(experiment.lc, features = "Cd19", reduction = "adt.umap")
SaveH5Seurat(experiment.lc, filename = "experiment.lc", overwrite = TRUE)
?RunUMAP


