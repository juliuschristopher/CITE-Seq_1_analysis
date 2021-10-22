## ##CITE-Seq 1 Script####
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

#Alter working capacity
plan()

plan("multiprocess", workers = 4)
options(future.globals.maxSize= 10097152000) # 10Gb

#Initial setup of colour palettes
col = colorRampPalette(brewer.pal(12, 'Set3'))(20)
colbig = colorRampPalette(brewer.pal(12, 'Set3'))(50)

####Load the 10X Cell Ranger output####
#Read the 10x Cell Ranger Output
c1_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Sample_GE_out/C1_GE/outs/filtered_feature_bc_matrix")
c2_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Sample_GE_out/C2_GE/outs/filtered_feature_bc_matrix")
d1_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Sample_GE_out/D1_GE/outs/filtered_feature_bc_matrix")
d2_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Sample_GE_out/D2_GE/outs/filtered_feature_bc_matrix")
a_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/A_WT_GE/outs/filtered_feature_bc_matrix")
b_ge.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/GE/B_WT_GE/outs/filtered_feature_bc_matrix")


#Add the sample to the cell names, consistent with antibody data below
colnames(c1_ge.data)=gsub("-1","_c1",colnames(c1_ge.data))
colnames(c2_ge.data)=gsub("-1","_c2",colnames(c2_ge.data))
colnames(d1_ge.data)=gsub("-1","_d1",colnames(d1_ge.data))
colnames(d2_ge.data)=gsub("-1","_d2",colnames(d2_ge.data))
colnames(a_ge.data)=gsub("-1","_a",colnames(a_ge.data))
colnames(b_ge.data)=gsub("-1","_b",colnames(b_ge.data))


#Uppercase the gene names for easier matching later
rownames(c1_ge.data)=toupper(rownames(c1_ge.data))
rownames(c2_ge.data)=toupper(rownames(c2_ge.data))
rownames(d1_ge.data)=toupper(rownames(d1_ge.data))
rownames(d2_ge.data)=toupper(rownames(d2_ge.data))
rownames(a_ge.data)=toupper(rownames(a_ge.data))
rownames(b_ge.data)=toupper(rownames(b_ge.data))


head(b_ge.data)
####Load 10X Antibody data####
#Read the 10x Antibody output
c1_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Antibody_fraction/C1_SP_out_2/umi_count",gene.column=1)
c2_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Antibody_fraction/C2_SP_out_2/umi_count",gene.column=1)
d1_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Antibody_fraction/D1_SP_out_2/umi_count",gene.column=1)
d2_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/Antibody_fraction/D2_SP_out_2/umi_count",gene.column=1)
a_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/A_WT/umi_count",gene.column=1)
b_ab.data <- Read10X(data.dir = "~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/second_batch_data_CP/B_WT/umi_count",gene.column=1)


#Tidy up the rownames from the data
rownames(c1_ab.data)=gsub("-[^-]+$","",rownames(c1_ab.data),perl=TRUE)
rownames(c2_ab.data)=gsub("-[^-]+$","",rownames(c2_ab.data),perl=TRUE)
rownames(d1_ab.data)=gsub("-[^-]+$","",rownames(d1_ab.data),perl=TRUE)
rownames(d2_ab.data)=gsub("-[^-]+$","",rownames(d2_ab.data),perl=TRUE)
rownames(a_ab.data)=gsub("-[^-]+$","",rownames(a_ab.data),perl=TRUE)
rownames(b_ab.data)=gsub("-[^-]+$","",rownames(b_ab.data),perl=TRUE)

#Add the Sample to the cell names in each sample
colnames(c1_ab.data)=paste(colnames(c1_ab.data),"_c1",sep="")
colnames(c2_ab.data)=paste(colnames(c2_ab.data),"_c2",sep="")
colnames(d1_ab.data)=paste(colnames(d1_ab.data),"_d1",sep="")
colnames(d2_ab.data)=paste(colnames(d2_ab.data),"_d2",sep="")
colnames(a_ab.data)=paste(colnames(a_ab.data),"_a",sep="")
colnames(b_ab.data)=paste(colnames(b_ab.data),"_b",sep="")

head(b_ab.data)
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


head(b[[]])

####Incoperate VDJ data####
#Load contig file
c1_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/C1_VDJ/outs/filtered_contig_annotations.csv")
c2_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/C2_VDJ/outs/filtered_contig_annotations.csv")
d1_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/D1_VDJ/outs/filtered_contig_annotations.csv")
d2_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/D2_VDJ/outs/filtered_contig_annotations.csv")
a_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/A_WT_VDJ/outs/filtered_contig_annotations.csv")
b_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/B_WT_VDJ/outs/filtered_contig_annotations.csv")

#match barcode names with GE and ADT data
c1_cl.data$barcode=gsub("-1","_c1",c1_cl.data$barcode)
c2_cl.data$barcode=gsub("-1","_c2",c2_cl.data$barcode)
d1_cl.data$barcode=gsub("-1","_d1",d1_cl.data$barcode)
d2_cl.data$barcode=gsub("-1","_d2",d2_cl.data$barcode)
a_cl.data$barcode=gsub("-1","_a",a_cl.data$barcode)
b_cl.data$barcode=gsub("-1","_b",b_cl.data$barcode)

contig_list <- list(c1_cl.data, c2_cl.data, d1_cl.data, d2_cl.data, a_cl.data, b_cl.data)
head(contig_list[[1]])

#Generate combined object
combined <- combineBCR(contig_list, samples = c("c1", "c2", "d1", "d2", "a", "b"))
combined[[1]]

str(combined)
head(combined[[1]])

#Make sure barcodes are identical to GE and ADT data
combined$c1$barcode=gsub("c1_","",combined$c1$barcode)
combined$c2$barcode=gsub("c2_","",combined$c2$barcode)
combined$d1$barcode=gsub("d1_","",combined$d1$barcode)
combined$d2$barcode=gsub("d2_","",combined$d2$barcode)
combined$a$barcode=gsub("a_","",combined$a$barcode)
combined$b$barcode=gsub("b_","",combined$b$barcode)

head(combined$c1)
#####Process samples as one####
experiments=c(c1,c2,d1,d2,a,b)
experiment_names=c("c1","c2","d1","d2","a","b")

experiment<-merge(x= c1, y=c(c2,d1,d2,a,b))

experiment
str(experiment)
head(experiment[[]])

####Merge seurat object with VDJ data####
experiment <- combineExpression(combined,
                                experiment,
                                cloneCall="gene", group.by = "sample")

head(experiment[[]])

####Quality control, filtering, normalisation and scaling####
#Mitochondrial QC metrics
experiment[["percent.mt"]] <- PercentageFeatureSet(experiment, pattern = "^MT-")

#Remove where nCount_ADT = 0
DefaultAssay(experiment) <- "ADT"
experiment <- subset(experiment, nCount_ADT >0)
DefaultAssay(experiment) <- "RNA"

#Visualize QC metrics as violin plot
RNA_QC <- VlnPlot(experiment, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ADT_QC <- VlnPlot(experiment, features = c("nFeature_ADT", "nCount_ADT", "percent.mt"))

RNA_QC
ADT_QC
#FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
feature_count.RNA_vs_percent.mt = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()  +
  ylab("% of mitochondrial genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 5) 

feature_count.RNA_vs_feature.RNA = FeatureScatter(experiment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() +
  ylab("Number of genes") +
  xlab("UMI counts") + 
  geom_hline(yintercept = 200) 

feature_count.RNA_vs_percent.mt + feature_count.RNA_vs_feature.RNA

#Generally aim to filter out unique feature counts over 2,500 and less than 200; and percent.mt over 5%
filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes and mitochondrial percentage.")
  seurat_object = subset(seurat_object, subset = nFeature_RNA > 200  & percent.mt < 5 & nFeature_RNA < 2500)
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  
  
  return(seurat_object)
}

experiment <-  filter_seurat(experiment)

####Normalise dataset - SCTransform####
experiment = SCTransform(experiment, verbose = TRUE)
experiment[["SCT"]]

top20 = head(VariableFeatures(experiment), 30)

plot1.1 = VariableFeaturePlot(experiment)
plot2.1 = LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
plot2.1

#Find variable features
experiment <- FindVariableFeatures(experiment, selection.method = "vst")
top20 = head(VariableFeatures(experiment), 20)

plot3 = VariableFeaturePlot(experiment)
plot4 = LabelPoints(plot = plot3, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
plot4

#Cell Cycle genes
S.genes = cc.genes.updated.2019$s.genes
G2M.genes = cc.genes.updated.2019$g2m.genes

experiment = CellCycleScoring(experiment, s.features=S.genes, g2m.features=G2M.genes, set.ident = TRUE)

####Dimensionality reduction -PCA####
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
experiment <- FindNeighbors(experiment, dims = 1:30)
experiment <- FindClusters(experiment, resolution = 1.0, verbose = FALSE)
experiment <- RunUMAP(experiment, dims = 1:30)
DimPlot(experiment, label = TRUE, cols=colbig) +  ggtitle("RNA Clustering")

####Scale antibody data####
DefaultAssay(experiment) <- "ADT"
VariableFeatures(experiment) <- rownames(experiment[["ADT"]])
experiment <- NormalizeData(experiment, normalization.method = "CLR", margin = 2)
experiment <- ScaleData(experiment)
experiment <- RunPCA(experiment,reduction.name = 'apca')

#Visualise antobody PCA
print(experiment[["apca"]], dims = 1:10, nfeatures = 5)
Plot_13 <- VizDimLoadings(experiment, dims = 1:4, reduction = "apca", nfeatures = 15)
Plot_13

####Combine into wnn plot####
experiment <- FindMultiModalNeighbors(
  experiment, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight")

#Determing resolution (cluster number) for wnn plot?
clustree(experiment, prefix = "wsnn_res.") +
  theme(legend.position="bottom")


#UMPA plots for RNA, ADT and WNN
experiment <- RunUMAP(experiment, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                      reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
experiment<- RunUMAP(experiment, reduction = 'apca', dims = 1:18, assay = 'ADT', 
                     reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
experiment <- RunUMAP(experiment, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
experiment <- FindClusters(experiment, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = TRUE)


DefaultAssay(experiment) <- "RNA"

p1=DimPlot(experiment, label = TRUE,cols=colbig,reduction = "rna.umap", label.size = 2.5) + NoLegend()
p2=DimPlot(experiment, label = TRUE,cols=colbig,reduction = "adt.umap", label.size = 2.5) + NoLegend()
p3=DimPlot(experiment, label = TRUE,cols=colbig, reduction = "wnn.umap", label.size = 2.5) + NoLegend()

p1|p2|p3

p1
p2
p3

###Umap-wnn by mouse
plot_mouse <- DimPlot(experiment, label = TRUE,reduction = "wnn.umap", label.size = 2.5, group.by = "orig.ident")
plot_mouse

###Umap-wnn by sample
DimPlot(experiment, label = TRUE,cols=colbig, reduction = "wnn.umap", label.size = 2.5, split.by = "orig.ident", ncol = 2) + NoLegend()

###Umap-wnn by cell cycle stage
DimPlot(experiment, label = TRUE,reduction = "wnn.umap", label.size = 2.5, group.by = "Phase")

head(experiment[[]])
###Match the RNA Names to the Antibodies, this should be checked
list1=c(rownames(a_ab.data))
list2=c("PTPRC","FAS","CD19","IGHM","CR2","FCER2A","CD93","CD83","CD86","IGHD","CD8A","SELL","CD44","CD4","CXCR5","PDCD1","IL2RA","CD274","PDCD1LG2","CTLA4","CD80","CD40","CD69","ICOS","CD38","TNFRSF18")

####Analysis of clusters####
##Different plotting options
DefaultAssay(experiment) <- "RNA"
DefaultAssay(experiment) <- "ADT"

FeaturePlot(experiment, features = c("CD19", "CD4", "CD8A", "PRDM1", "PPBP", "NKG7", "CST3", "FOXP3", "B220"), reduction = "wnn.umap")

RidgePlot(experiment, features = c("CD19", "CYP11A1"), ncol = 2)

FeaturePlot(experiment, features = c("IGHV1-53", "IGKV3-4", "IGHD1-1"), reduction = "wnn.umap")
FeaturePlot(experiment, feature = "IGHG", reduction = "wnn.umap")

FeaturePlot(experiment, features = c("CD44"), reduction = "wnn.umap")
VlnPlot(experiment, feature = "CTLA4")
?VlnPlot
p3

Idents(object = experiment) <- "old.ident"
head(experiment[[]])
##Finding all the markers
experiment.markers <- FindAllMarkers(experiment, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
experiment.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(experiment, features = top10$gene) + NoLegend()

##DE genes of individual clusters
Cluster_11 <- FindMarkers(experiment, ident.1 = 11, assay = "RNA")
Cluster_11_adt <- FindMarkers(experiment, ident.1 = 11, assay = "ADT")

Cluster_4 <- FindMarkers(experiment, ident.1 = 4, assay = "RNA")
Cluster_4_adt <- FindMarkers(experiment, ident.1 = 4, assay = "ADT")

##Subsetting unknown cluster
Unknown_cells <- subset(experiment, idents = c(2, 12, 27, 33))
Unknown_cells <- FindClusters(Unknown_cells, resolution = 0.8, verbose = FALSE, graph.name = "wsnn")
Unknown_cells <- RunUMAP(Unknown_cells, dims = 1:30, reduction.name = "unknown.umap")
DimPlot(Unknown_cells, label = TRUE, cols=colbig, reduction = "unknown.umap", label.size = 2.5) + NoLegend()
FeaturePlot(Unknown_cells, "SOX4", reduction = "unknown.umap")


####Clonotype analysis####
##Data visualization
#Percent/total number of unique clonotypes 
quantContig(combined, cloneCall = "gene+nt", scale = T) #percent of unique clonotypes of total size of the size of clonotyeps
quantContig(combined, cloneCall = "gene+nt", scale = F) #number of uniqe clonotypes

quantContig(combined, cloneCall = "gene+nt", scale = T, chain = "IGL")#by chain

#Abundance of clonotypes
Abundance_clonotypes <- abundanceContig(combined, cloneCall = "gene", scale = F, exportTable = T)
Abundance_clonotypes <- Abundance_clonotypes %>%
  arrange(desc(Abundance))
Abundance_clonotypes

abundanceContig(combined, cloneCall = "gene", scale = T)

#Length of clonotypes
lengthContig(combined, cloneCall = "aa")
lengthContig(combined, cloneCall = "nt")

#Compare clonotypes
compareClonotypes(combined, samples = c("a", "b"), cloneCall = "aa", graph = "alluvial") #Computationally intese

#Visualise Gene Usage
vizGenes(combined, gene = "D", chain = "IGH", plot = "bar", order = "variance", scale = TRUE)
vizGenes(combined, gene = "V", chain = "IGL", plot = "bar", order = "variance", scale = TRUE)

vizGenes(combined, gene = "V", chain = "IGL", plot = "heatmap", scale = TRUE, order = "gene")

#Clonal overlap
clonalOverlap(combined, cloneCall = "gene+nt", 
              method = "morisita")

#Clonotype proportion
clonalProportion(combined, cloneCall = "gene")
clonalProportion(combined, cloneCall = "nt")

#Clonal Homeostasis
clonalHomeostasis(combined, cloneCall = "gene")
clonalHomeostasis(combined, cloneCall = "nt")


