## ##CITE-Seq 1 Script - Cluster Analysis####
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

#Initial setup of colour palettes
col = colorRampPalette(brewer.pal(12, 'Set3'))(20)
colbig = colorRampPalette(brewer.pal(12, 'Set3'))(50)

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
FeaturePlot(experiment, features = c("APOE"), reduction = "wnn.umap")
VlnPlot(experiment, feature = "APOE")


####Finding all the markers####
experiment.markers <- FindAllMarkers(experiment, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
experiment.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(experiment, features = top10$gene) + NoLegend()

####DEG of each cluster####
###Find DEGs
Cluster_0 <- FindMarkers(experiment, ident.1 = 0, assay = "RNA")
Cluster_0_adt <- FindMarkers(experiment, ident.1 = 0, assay = "ADT")

Cluster_1 <- FindMarkers(experiment, ident.1 = 1, assay = "RNA")
Cluster_1_adt <- FindMarkers(experiment, ident.1 = 1, assay = "ADT")

Cluster_2 <- FindMarkers(experiment, ident.1 = 2, assay = "RNA")
Cluster_2_adt <- FindMarkers(experiment, ident.1 = 2, assay = "ADT")

Cluster_3 <- FindMarkers(experiment, ident.1 = 3, assay = "RNA")
Cluster_3_adt <- FindMarkers(experiment, ident.1 = 3, assay = "ADT")

Cluster_4 <- FindMarkers(experiment, ident.1 = 4, assay = "RNA")
Cluster_4_adt <- FindMarkers(experiment, ident.1 = 4, assay = "ADT")

Cluster_5 <- FindMarkers(experiment, ident.1 = 5, assay = "RNA")
Cluster_5_adt <- FindMarkers(experiment, ident.1 = 5, assay = "ADT")

Cluster_6 <- FindMarkers(experiment, ident.1 = 6, assay = "RNA")
Cluster_6_adt <- FindMarkers(experiment, ident.1 = 6, assay = "ADT")

Cluster_7 <- FindMarkers(experiment, ident.1 = 7, assay = "RNA")
Cluster_7_adt <- FindMarkers(experiment, ident.1 = 7, assay = "ADT")

Cluster_8 <- FindMarkers(experiment, ident.1 = 8, assay = "RNA")
Cluster_8_adt <- FindMarkers(experiment, ident.1 = 8, assay = "ADT")

Cluster_9 <- FindMarkers(experiment, ident.1 = 9, assay = "RNA")
Cluster_9_adt <- FindMarkers(experiment, ident.1 = 9, assay = "ADT")

Cluster_10 <- FindMarkers(experiment, ident.1 = 10, assay = "RNA")
Cluster_10_adt <- FindMarkers(experiment, ident.1 = 10, assay = "ADT")

Cluster_11 <- FindMarkers(experiment, ident.1 = 11, assay = "RNA")
Cluster_11_adt <- FindMarkers(experiment, ident.1 = 11, assay = "ADT")

Cluster_12 <- FindMarkers(experiment, ident.1 = 12, assay = "RNA")
Cluster_12_adt <- FindMarkers(experiment, ident.1 = 12, assay = "ADT")

Cluster_13 <- FindMarkers(experiment, ident.1 = 13, assay = "RNA")
Cluster_13_adt <- FindMarkers(experiment, ident.1 = 13, assay = "ADT")

Cluster_14 <- FindMarkers(experiment, ident.1 = 14, assay = "RNA")
Cluster_14_adt <- FindMarkers(experiment, ident.1 = 14, assay = "ADT")

Cluster_15 <- FindMarkers(experiment, ident.1 = 15, assay = "RNA")
Cluster_15_adt <- FindMarkers(experiment, ident.1 = 15, assay = "ADT")

Cluster_16 <- FindMarkers(experiment, ident.1 = 16, assay = "RNA")
Cluster_16_adt <- FindMarkers(experiment, ident.1 = 16, assay = "ADT")

Cluster_17 <- FindMarkers(experiment, ident.1 = 17, assay = "RNA")
Cluster_17_adt <- FindMarkers(experiment, ident.1 = 17, assay = "ADT")

Cluster_18 <- FindMarkers(experiment, ident.1 = 18, assay = "RNA")
Cluster_18_adt <- FindMarkers(experiment, ident.1 = 18, assay = "ADT")

Cluster_19 <- FindMarkers(experiment, ident.1 = 19, assay = "RNA")
Cluster_19_adt <- FindMarkers(experiment, ident.1 = 19, assay = "ADT")

Cluster_20 <- FindMarkers(experiment, ident.1 = 20, assay = "RNA")
Cluster_20_adt <- FindMarkers(experiment, ident.1 = 20, assay = "ADT")

Cluster_21 <- FindMarkers(experiment, ident.1 = 21, assay = "RNA")
Cluster_21_adt <- FindMarkers(experiment, ident.1 = 21, assay = "ADT")

Cluster_22 <- FindMarkers(experiment, ident.1 = 22, assay = "RNA")
Cluster_22_adt <- FindMarkers(experiment, ident.1 = 22, assay = "ADT")

Cluster_23 <- FindMarkers(experiment, ident.1 = 23, assay = "RNA")
Cluster_23_adt <- FindMarkers(experiment, ident.1 = 23, assay = "ADT")

Cluster_24 <- FindMarkers(experiment, ident.1 = 24, assay = "RNA")
Cluster_24_adt <- FindMarkers(experiment, ident.1 = 24, assay = "ADT")

Cluster_25 <- FindMarkers(experiment, ident.1 = 25, assay = "RNA")
Cluster_25_adt <- FindMarkers(experiment, ident.1 = 25, assay = "ADT")

Cluster_26 <- FindMarkers(experiment, ident.1 = 26, assay = "RNA")
Cluster_26_adt <- FindMarkers(experiment, ident.1 = 26, assay = "ADT")

Cluster_27 <- FindMarkers(experiment, ident.1 = 27, assay = "RNA")
Cluster_27_adt <- FindMarkers(experiment, ident.1 = 27, assay = "ADT")

Cluster_28 <- FindMarkers(experiment, ident.1 = 28, assay = "RNA")
Cluster_28_adt <- FindMarkers(experiment, ident.1 = 28, assay = "ADT")

Cluster_29 <- FindMarkers(experiment, ident.1 = 29, assay = "RNA")
Cluster_29_adt <- FindMarkers(experiment, ident.1 = 29, assay = "ADT")

Cluster_30 <- FindMarkers(experiment, ident.1 = 30, assay = "RNA")
Cluster_30_adt <- FindMarkers(experiment, ident.1 = 30, assay = "ADT")

Cluster_31 <- FindMarkers(experiment, ident.1 = 31, assay = "RNA")
Cluster_31_adt <- FindMarkers(experiment, ident.1 = 31, assay = "ADT")

Cluster_32 <- FindMarkers(experiment, ident.1 = 32, assay = "RNA")
Cluster_32_adt <- FindMarkers(experiment, ident.1 = 32, assay = "ADT")

Cluster_33 <- FindMarkers(experiment, ident.1 = 33, assay = "RNA")
Cluster_33_adt <- FindMarkers(experiment, ident.1 = 33, assay = "ADT")

Cluster_34 <- FindMarkers(experiment, ident.1 = 34, assay = "RNA")
Cluster_34_adt <- FindMarkers(experiment, ident.1 = 34, assay = "ADT")

Cluster_35 <- FindMarkers(experiment, ident.1 = 35, assay = "RNA")
Cluster_35_adt <- FindMarkers(experiment, ident.1 = 35, assay = "ADT")

Cluster_36 <- FindMarkers(experiment, ident.1 = 36, assay = "RNA")
Cluster_36_adt <- FindMarkers(experiment, ident.1 = 36, assay = "ADT")

Cluster_37 <- FindMarkers(experiment, ident.1 = 37, assay = "RNA")
Cluster_37_adt <- FindMarkers(experiment, ident.1 = 37, assay = "ADT")

Cluster_38 <- FindMarkers(experiment, ident.1 = 38, assay = "RNA")
Cluster_38_adt <- FindMarkers(experiment, ident.1 = 38, assay = "ADT")

Cluster_39 <- FindMarkers(experiment, ident.1 = 39, assay = "RNA")
Cluster_39_adt <- FindMarkers(experiment, ident.1 = 39, assay = "ADT")

Cluster_40 <- FindMarkers(experiment, ident.1 = 40, assay = "RNA")
Cluster_40_adt <- FindMarkers(experiment, ident.1 = 40, assay = "ADT")

###Export Cluster-associated gene lists
Cluster_4_exp <- tibble::rownames_to_column(Cluster_4, "Genes")
write_xlsx(Cluster_4_exp, "Cluster_gene_signatures/Cluster_4/Cluster_4.xlsx")

####TFIDF####
###Define the function
tfidf = function(data,target,universe){
  if(!all(target %in% universe))
    stop('Target must be a subset of universe')
  nObs = Matrix::rowSums(data[,target,drop=FALSE]>0)
  nTot = Matrix::rowSums(data[,universe,drop=FALSE]>0)
  tf = nObs/length(target)
  idf = log(length(universe)/nTot)
  score = tf*idf
  #Calculate p-value for significance based on using a hypergeometric distribution to simulate the results of infinite random sampling
  pvals = phyper(nObs-1,nTot,length(universe)-nTot,length(target),lower.tail=FALSE)
  qvals = p.adjust(pvals,method='BH')
  ntf = (exp(-idf)*length(universe)-tf*length(target))/(length(universe)-length(target))
  return(data.frame(geneFrequency=tf,
                    geneFrequencyOutsideCluster=ntf,
                    geneFrequencyGlobal=exp(-idf),
                    geneExpression=Matrix::rowMeans(data[,target,drop=FALSE]),
                    geneExpressionOutsideCluster = Matrix::rowMeans(data[,universe[!(universe%in%target)],drop=FALSE]),
                    geneExpressionGlobal = Matrix::rowMeans(data),
                    idf=idf,
                    tfidf=score,
                    qval=qvals)[order(score,decreasing=TRUE),])
}

###Select cluster
TFIDF.c0 <- WhichCells(object = experiment, ident = 0)
TFIDF.c1 <- WhichCells(object = experiment, ident = 1)
TFIDF.c2 <- WhichCells(object = experiment, ident = 2)
TFIDF.c3 <- WhichCells(object = experiment, ident = 3)
TFIDF.c4 <- WhichCells(object = experiment, ident = 4)
TFIDF.c5 <- WhichCells(object = experiment, ident = 5)
TFIDF.c6 <- WhichCells(object = experiment, ident = 6)
TFIDF.c7 <- WhichCells(object = experiment, ident = 7)
TFIDF.c8 <- WhichCells(object = experiment, ident = 8)
TFIDF.c9 <- WhichCells(object = experiment, ident = 9)
TFIDF.c10 <- WhichCells(object = experiment, ident = 10)
TFIDF.c11 <- WhichCells(object = experiment, ident = 11)
TFIDF.c12 <- WhichCells(object = experiment, ident = 12)
TFIDF.c13 <- WhichCells(object = experiment, ident = 13)
TFIDF.c14 <- WhichCells(object = experiment, ident = 14)
TFIDF.c15 <- WhichCells(object = experiment, ident = 15)
TFIDF.c16 <- WhichCells(object = experiment, ident = 16)
TFIDF.c17 <- WhichCells(object = experiment, ident = 17)
TFIDF.c18 <- WhichCells(object = experiment, ident = 18)
TFIDF.c19 <- WhichCells(object = experiment, ident = 19)
TFIDF.c20 <- WhichCells(object = experiment, ident = 20)
TFIDF.c21 <- WhichCells(object = experiment, ident = 21)
TFIDF.c22 <- WhichCells(object = experiment, ident = 22)
TFIDF.c23 <- WhichCells(object = experiment, ident = 23)
TFIDF.c24 <- WhichCells(object = experiment, ident = 24)
TFIDF.c25 <- WhichCells(object = experiment, ident = 25)
TFIDF.c26 <- WhichCells(object = experiment, ident = 26)
TFIDF.c27 <- WhichCells(object = experiment, ident = 27)
TFIDF.c28 <- WhichCells(object = experiment, ident = 28)
TFIDF.c29 <- WhichCells(object = experiment, ident = 29)
TFIDF.c30 <- WhichCells(object = experiment, ident = 30)
TFIDF.c31 <- WhichCells(object = experiment, ident = 31)
TFIDF.c32 <- WhichCells(object = experiment, ident = 32)
TFIDF.c33 <- WhichCells(object = experiment, ident = 33)
TFIDF.c34 <- WhichCells(object = experiment, ident = 34)
TFIDF.c35 <- WhichCells(object = experiment, ident = 35)
TFIDF.c36 <- WhichCells(object = experiment, ident = 36)
TFIDF.c37 <- WhichCells(object = experiment, ident = 37)
TFIDF.c38 <- WhichCells(object = experiment, ident = 38)
TFIDF.c39 <- WhichCells(object = experiment, ident = 39)
TFIDF.c40 <- WhichCells(object = experiment, ident = 40)

###Set to RNA assay
DefaultAssay(experiment)<-"RNA"

##Run function - Higher TFIDF score and lower q value genes are more uniquely expressed
TFIDF.c0.genes <- tfidf(GetAssayData(experiment), TFIDF.c0, colnames(experiment))
TFIDF.c1.genes <- tfidf(GetAssayData(experiment), TFIDF.c1, colnames(experiment))
TFIDF.c2.genes <- tfidf(GetAssayData(experiment), TFIDF.c2, colnames(experiment))
TFIDF.c3.genes <- tfidf(GetAssayData(experiment), TFIDF.c3, colnames(experiment))
TFIDF.c4.genes <- tfidf(GetAssayData(experiment), TFIDF.c4, colnames(experiment))
TFIDF.c5.genes <- tfidf(GetAssayData(experiment), TFIDF.c5, colnames(experiment))
TFIDF.c6.genes <- tfidf(GetAssayData(experiment), TFIDF.c6, colnames(experiment))
TFIDF.c7.genes <- tfidf(GetAssayData(experiment), TFIDF.c7, colnames(experiment))
TFIDF.c8.genes <- tfidf(GetAssayData(experiment), TFIDF.c8, colnames(experiment))
TFIDF.c9.genes <- tfidf(GetAssayData(experiment), TFIDF.c9, colnames(experiment))
TFIDF.c10.genes <- tfidf(GetAssayData(experiment), TFIDF.c10, colnames(experiment))
TFIDF.c11.genes <- tfidf(GetAssayData(experiment), TFIDF.c11, colnames(experiment))
TFIDF.c13.genes <- tfidf(GetAssayData(experiment), TFIDF.c13, colnames(experiment))
TFIDF.c14.genes <- tfidf(GetAssayData(experiment), TFIDF.c14, colnames(experiment))
TFIDF.c15.genes <- tfidf(GetAssayData(experiment), TFIDF.c15, colnames(experiment))
TFIDF.c16.genes <- tfidf(GetAssayData(experiment), TFIDF.c16, colnames(experiment))
TFIDF.c17.genes <- tfidf(GetAssayData(experiment), TFIDF.c17, colnames(experiment))
TFIDF.c18.genes <- tfidf(GetAssayData(experiment), TFIDF.c18, colnames(experiment))
TFIDF.c19.genes <- tfidf(GetAssayData(experiment), TFIDF.c19, colnames(experiment))
TFIDF.c20.genes <- tfidf(GetAssayData(experiment), TFIDF.c20, colnames(experiment))
TFIDF.c21.genes <- tfidf(GetAssayData(experiment), TFIDF.c21, colnames(experiment))
TFIDF.c22.genes <- tfidf(GetAssayData(experiment), TFIDF.c22, colnames(experiment))
TFIDF.c23.genes <- tfidf(GetAssayData(experiment), TFIDF.c23, colnames(experiment))
TFIDF.c24.genes <- tfidf(GetAssayData(experiment), TFIDF.c24, colnames(experiment))
TFIDF.c25.genes <- tfidf(GetAssayData(experiment), TFIDF.c25, colnames(experiment))
TFIDF.c26.genes <- tfidf(GetAssayData(experiment), TFIDF.c26, colnames(experiment))
TFIDF.c27.genes <- tfidf(GetAssayData(experiment), TFIDF.c27, colnames(experiment))
TFIDF.c28.genes <- tfidf(GetAssayData(experiment), TFIDF.c28, colnames(experiment))
TFIDF.c29.genes <- tfidf(GetAssayData(experiment), TFIDF.c29, colnames(experiment))
TFIDF.c30.genes <- tfidf(GetAssayData(experiment), TFIDF.c30, colnames(experiment))
TFIDF.c31.genes <- tfidf(GetAssayData(experiment), TFIDF.c31, colnames(experiment))
TFIDF.c32.genes <- tfidf(GetAssayData(experiment), TFIDF.c32, colnames(experiment))
TFIDF.c33.genes <- tfidf(GetAssayData(experiment), TFIDF.c33, colnames(experiment))
TFIDF.c34.genes <- tfidf(GetAssayData(experiment), TFIDF.c34, colnames(experiment))
TFIDF.c35.genes <- tfidf(GetAssayData(experiment), TFIDF.c35, colnames(experiment))
TFIDF.c36.genes <- tfidf(GetAssayData(experiment), TFIDF.c36, colnames(experiment))
TFIDF.c37.genes <- tfidf(GetAssayData(experiment), TFIDF.c37, colnames(experiment))
TFIDF.c38.genes <- tfidf(GetAssayData(experiment), TFIDF.c38, colnames(experiment))
TFIDF.c39.genes <- tfidf(GetAssayData(experiment), TFIDF.c39, colnames(experiment))
TFIDF.c40.genes <- tfidf(GetAssayData(experiment), TFIDF.c40, colnames(experiment))

###Export table
TFIDF_cluster_4 <- tibble::rownames_to_column(TFIDF.c4.genes, "Genes")
write_xlsx(TFIDF_cluster_4, "Cluster_gene_signatures/Cluster_4/TFIDF_cluster_4.xlsx")

####Addmodulescore####
###Load gene signature file
Hallmark_sig <- read.csv("Population_signatures/Hallmark_signatures.csv",header = T, sep = ',')

###Load individual signature lists
Apoptosis_list <- list(Hallmark_sig$Apoptosis)
PI3Ksig_list <- list(Hallmark_sig$PI3K_AKT_MTOR_Signalling)
Cholesterol_list <- list(Hallmark_sig$Cholesterol_Homeostasis)
IL2sig_list <- list(Hallmark_sig$IL2_STAT5_Signalling)
IL6sig_list <- list(Hallmark_sig$IL6_JAK_STAT3_Signalling)


###convertMouseGeneList - function
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  mousex <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(mousex)
}

###Convert human to mouse genes
Apoptosis_list <- convertHumanGeneList(Apoptosis_list)
lapply(Apoptosis_list, toupper)
str(Apoptosis_list)

###Set correct assay
DefaultAssay(experiment) <- "RNA"

##A#ddmodulescore
experiment <-AddModuleScore(experiment, features = Apoptosis_list, name = "Apoptosis_enrichment")
experiment <-AddModuleScore(experiment, features = PI3Ksig_list, name = "PI3Ksig_enrichment")
experiment <-AddModuleScore(experiment, features = Cholesterol_list, name = "Cholesterol_enrichment")
experiment <-AddModuleScore(experiment, features = IL2sig_list, name = "IL2sig_enrichment")
experiment <-AddModuleScore(experiment, features = IL6sig_list, name = "IL6sig_enrichment")

###Visualise enrichment scores
VlnPlot(experiment, c("IL6sig_enrichment1"), pt.size  = 0.1)+NoLegend()
FeaturePlot(experiment, c("PI3Ksig_enrichment1"),cols=c("Blue", "Yellow"), reduction = "wnn.umap")

####Subsetting of clusters####
B_cell_c4 <- subset(experiment, idents = 4)

B_cell_c4 = SCTransform(B_cell_c4, verbose = TRUE)
top20_c4 = head(VariableFeatures(B_cell_c4), 20)
plot_top20_c4 = VariableFeaturePlot(B_cell_c4)
plot_top20_c4 = LabelPoints(plot = plot_top20_c4, points = top20_c4, repel = TRUE, xnudge = 0, ynudge = 0)

B_cell_c4 <- RunPCA(B_cell_c4, verbose = FALSE, features = VariableFeatures(object = B_cell_c4))
print(B_cell_c4[["pca"]], dims = 1:5, nfeatures = 5)
pca_variance_c4 <- B_cell_c4@reductions$pca@stdev^2
plot(pca_variance_c4/sum(pca_variance_c4), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01)


B_cell_c4 <- FindNeighbors(B_cell_c4, dims = 1:20)
B_cell_c4 <- FindClusters(B_cell_c4, resolution = 0.1, verbose = FALSE)
B_cell_c4 <- RunUMAP(B_cell_c4, reduction = 'pca', dims = 1:23, reduction.name = "B_cell_c4.umap")
DimPlot(B_cell_c4, label = TRUE, cols=colbig) +  ggtitle("RNA Clustering")



