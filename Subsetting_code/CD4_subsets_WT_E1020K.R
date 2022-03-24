## ##Subsetting: CD4+ cells####
#Julius Christopher Baeck
#CITE-Seq (1) analysis - WT and E1020K
#Subsetting of all T cells

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
library(ggpubr)
library(data.table)

####Define functions####
###TFIDF
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

#Stacked violin plot
#remove the x-axis text and tick
#plot.margin to adjust the white space between each plot.
#... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


#main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  #Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  #change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


####Define colour palettes####
col_con = viridis(50)
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)


####Load Seurat object (lower case gene names - mouse)####
experiment.lc.plain <- LoadH5Seurat("experiment.lc.plain.h5seurat")
head(experiment.lc.plain[[]])

####Subsetting####
test <- experiment.lc.plain
head(test[[]])

##Subset based on Genoyype: WT and E1020K
DefaultAssay(test) <-  "RNA"
CD4cells <- subset(test, subset = cloneType == "NA") #16119
CD4cells <- subset(CD4cells, subset = orig.ident != "c1")
CD4cells <- subset(CD4cells, subset = orig.ident != "c2")
CD4cells <- subset(CD4cells, subset = orig.ident != "d1")
CD4cells <- subset(CD4cells, subset = orig.ident != "d2") #7351

##Subset based on x3 Cd4 ADT expression and 0x Cd19 and 0x Cd8a expression
DefaultAssay(CD4cells) <-  "ADT"
CD4_cells <- subset(CD4cells, subset = Cd4 > 3) #4739
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <- subset(CD4_cells, subset = Cd4 > 0) #1616
CD4_cells <- subset(CD4_cells, subset = Cd19 < 1) #1606
CD4_cells <- subset(CD4_cells, subset = Cd8a < 1) #1583
head(CD4_cells[[]])


###RNA####
###Normalise subset###
DefaultAssay(CD4_cells) <- "RNA" #For log normalisation
DefaultAssay(CD4_cells) <- "SCT" #For SCTransform

##RNA normalisation
CD4_cells <- NormalizeData(CD4_cells, verbose = TRUE)
CD4_cells <- FindVariableFeatures(CD4_cells, nfeatures = 3000)
CD4_cells <- ScaleData(CD4_cells)

#Or
CD4_cells <-  SCTransform(CD4_cells, verbose = TRUE)
CD4_cells[["SCT"]]

##Visualisation
top20 <-  head(VariableFeatures(CD4_cells), 20)
plot1.1 <-  VariableFeaturePlot(CD4_cells)
top20_plot <-  LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)

##RNA PCA
CD4_cells <- RunPCA(CD4_cells, verbose = FALSE, features = VariableFeatures(object = CD4_cells))
pca_variance <- CD4_cells@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #30

##RNA clustering
DefaultAssay(CD4_cells) <- "RNA" #For log normalisation
DefaultAssay(CD4_cells) <- "SCT" #For SCTransform

CD4_cells <- FindNeighbors(CD4_cells, dims = 1:30)
CD4_cells <- FindClusters(CD4_cells, resolution = 0.5, verbose = FALSE) #0.5 for the resolution
clustree(CD4_cells, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
CD4_cells <-RunUMAP(CD4_cells, dims = 1:30, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
CD4_cells_p1 <- DimPlot(CD4_cells, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("RNA Clustering") + theme_void() + NoLegend()
CD4_cells_p1 <- CD4_cells_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####ADT####
DefaultAssay(CD4_cells) <- "ADT"

##ADT normalisation
VariableFeatures(CD4_cells) <- rownames(CD4_cells[["ADT"]])
CD4_cells <- NormalizeData(CD4_cells, normalization.method = "CLR", margin = 2)
CD4_cells <- ScaleData(CD4_cells)

##ADT PCA
CD4_cells <- RunPCA(CD4_cells, reduction.name = 'apca', approx = FALSE)
apca_variance <- CD4_cells@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #26

##ADT clustering
CD4_cells <- FindNeighbors(CD4_cells, dims = 1:26, reduction = "apca")
CD4_cells <- FindClusters(CD4_cells, resolution = 1.0, verbose = FALSE) #1.0 for the resolution
clustree(CD4_cells, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
CD4_cells <- RunUMAP(CD4_cells, reduction = 'apca', dims = 1:26, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
CD4_cells_p2 <- DimPlot(CD4_cells, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("ADT Clustering") + theme_void() + NoLegend()
CD4_cells_p2 <- CD4_cells_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####WNN####
DefaultAssay(CD4_cells) <- "RNA" #For log normalisation
DefaultAssay(CD4_cells) <- "SCT" #For SCTransform

##Combine into wnn plot
CD4_cells <- FindMultiModalNeighbors(
  CD4_cells, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:26), modality.weight.name = "RNA.weight")

##WNN clustering
CD4_cells <- FindClusters(CD4_cells, graph.name = "wsnn", algorithm = 3, resolution = 0.9, verbose = TRUE) #0.9 for the resolution
clustree(CD4_cells, prefix = "wsnn_res.") + theme(legend.position="bottom")
CD4_cells <- RunUMAP(CD4_cells, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
CD4_cells_p3 <- DimPlot(CD4_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Seurat Clusters") + theme_classic() + NoLegend()
CD4_cells_p3 <- CD4_cells_p3 + theme(plot.title = element_text(color="black", size=15)) + xlab("UMAP1") + ylab("UMAP2")

##Change Idents
Idents(CD4_cells)
CD4_cells[["old.ident"]] <- Idents(CD4_cells)
Idents(CD4_cells) <- CD4_cells[["orig.ident"]]
CD4_cells<- RenameIdents(CD4_cells, `a` = "WT 1", `b` = "WT 2", `f` = "E1020K")
CD4_cells[["orig.ident"]] <- Idents(CD4_cells)
Idents(CD4_cells) <- CD4_cells[["old.ident"]]
Idents(CD4_cells)

CD4_cells
head(CD4_cells[[]])
SaveH5Seurat(CD4_cells, filename = "CD4_cells_WT_and_E1020K", overwrite = TRUE)
CD4_cells <- LoadH5Seurat("CD4_cells_WT_and_E1020K.h5seurat")
Idents(CD4_cells) <- CD4_cells$wsnn_res.0.9
CD4_cells_p3 <- DimPlot(CD4_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Seurat Clusters") + theme_classic() + NoLegend()
CD4_cells_p3 <- CD4_cells_p3 + theme(plot.title = element_text(color="black", size=15, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")



####Subset analysis####
DefaultAssay(CD4_cells) <- "ADT"
DefaultAssay(CD4_cells) <- "RNA"

##Set seurat.clusters to wnn.umap clusters
Idents(CD4_cells) <- CD4_cells$wsnn_res.0.9
CD4_cells<- RenameIdents(CD4_cells, `0` = "Cluster 0", `1` ="Cluster 1", `2` ="Cluster 2", `3` = "Cluster 3", `4` = "Cluster 4",`5` ="Cluster 5", `6` = "Cluster 6", `7` = "Cluster 7", `8` = "Cluster 8", `9` = "Cluster 9")
CD4_cells$seurat_clusters <- Idents(CD4_cells)
Idents(CD4_cells) <- CD4_cells$wsnn_res.0.9
CD4_cells_p3

##By genotype
CD4_cells_Mouse <- DimPlot(CD4_cells, label = FALSE ,reduction = "wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, label.box = FALSE, repel = FALSE, cols = c("#F8766D", "#B79F00", "#7CAE00"))
CD4_cells_Mouse <- CD4_cells_Mouse + theme_void() +
  ggtitle("Highlighted by genotype") +
  theme(plot.title = element_text(color="black", size=15, face="bold"))

show_col(hue_pal()(12))

Sample.a <- subset(CD4_cells, subset = Genotype == "WT 1")
Sample.a.plot <- DimPlot(Sample.a, label = FALSE ,reduction = "wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#F8766D")) +
  theme_void() +
  ggtitle("WT 1") +
  theme(plot.title = element_text(color="black", size=15, face="bold"))

Sample.b <- subset(CD4_cells, subset = Genotype == "WT 2")
Sample.b.plot <- DimPlot(Sample.b, label = FALSE ,reduction = "wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#B79F00")) +
  theme_void() +
  ggtitle("WT 2") +
  theme(plot.title = element_text(color="black", size=15, face="bold"))

Sample.f <- subset(CD4_cells, subset = Genotype == "E1020K")
Sample.f.plot <- DimPlot(Sample.f, label = FALSE ,reduction = "wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#7CAE00")) +
  theme_void() +
  ggtitle("E1020K") +
  theme(plot.title = element_text(color="black", size=15, face="bold"))

Sample.a.b <- subset(CD4_cells, subset = Genotype != "E1020K")
Sample.a.b.plot <- DimPlot(Sample.a.b, label = FALSE ,reduction = "wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#F8766D","#B79F00")) +
  theme_void() +
  ggtitle("WT combined") +
  theme(plot.title = element_text(color="black", size=15, face="bold"))

CD4_cells_p3 + CD4_cells_cell_cylce + Sample.a.b.plot + Sample.f.plot

##By cell cycle genes
S.genes <- cc.genes.updated.2019$s.genes
S.genes <- lapply(S.genes, str_to_title)
G2M.genes <-  cc.genes.updated.2019$g2m.genes
G2M.genes <- lapply(G2M.genes, str_to_title)
CD4_cells <- CellCycleScoring(CD4_cells, s.features=S.genes, g2m.features=G2M.genes, set.ident = FALSE)
Idents(CD4_cells)
head(CD4_cells[[]])


CD4_cells_cell_cylce <- DimPlot(CD4_cells, label = FALSE, reduction = "wnn.umap", group.by = "Phase", pt.size = 1.3) +  ggtitle("Highlighted by cell-cycle stage") + theme_void()
CD4_cells_cell_cylce <- CD4_cells_cell_cylce + theme(plot.title = element_text(color="black", size=15, face="bold"))

CD4_cells_p3 + CD4_cells_cell_cylce + CD4_cells_Mouse

##Cell numbers and percent by genotype by clusters - harmony.wnn)
cell.numbers <- table(CD4_cells@meta.data$seurat_clusters, CD4_cells@meta.data$Genotype)
cell.numbers <- as.data.frame.matrix(cell.numbers)

cell_number_heatmap <- pheatmap::pheatmap(t(cell.numbers), cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                                          cellwidth = 30,cellheight = 30, angle_col = 45)

CD4_cells_meta <- CD4_cells@meta.data
genotype.numbers <- CD4_cells_meta %>% dplyr::count(Genotype)
genotype.numbers.vector <- genotype.numbers %>% pull(n)
str(genotype.numbers.vector)
cell.percent <- sweep(cell.numbers, 2, genotype.numbers.vector, "/")
cell.percent <- cell.percent*100

cell_percent_heatmap <- pheatmap::pheatmap(t(cell.percent), cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                                           cellwidth = 30,cellheight = 30, angle_col = 45)

cell.numbers <- tibble::rownames_to_column(cell.numbers, "Cluster")
cell.percent <- tibble::rownames_to_column(cell.percent, "Population")
write_xlsx(cell.numbers, "cell.numbers.xlsx")
write_xlsx(cell.percent, "cell.percent.xlsx")

CD4_cells_p3

cell.percent.WT1 <- cell.percent[, 0:2]
colnames(cell.percent) <- c("Population", "WT1", "WT2", "E1020K")
rownames(cell.percent) <- cell.percent$Population
cell.percent <- cell.percent[, 2:4] %>%
  round(2)
cell.percent[] <- paste0(as.matrix(cell.percent), "%")
cell.percent <- tibble::rownames_to_column(cell.percent, "Population")


#Pie charts
pie_labels <- paste0(round(100 * count_2/sum(count_2), 2), "%")

pie(cell.percent$WT1, labels = cell.percent$Population, main="Pie Chart of WT 1 Populations", col = col2)
pie(cell.percent$WT2, labels = cell.percent$Population, main="Pie Chart of WT 2 Populations", col = col2)
pie(cell.percent$E1020K, labels = cell.percent$Population, main="Pie Chart of E1020K Populations", col = col2)

ggplot(cell.percent, aes(x = "", y = WT1, fill = Population)) +
  geom_col(color = "black") +
  geom_label(aes(label = WT1),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  guides(fill = guide_legend(title = "Populations")) +
  coord_polar(theta = "y") + 
  theme_void()


#Stacked bar charts
cell.percent <- as.data.frame(cell.percent)
str(cell.percent)
unstack(cell.percent)
rownames(cell.percent) <- cell.percent$Population
cell.percent <- cell.percent[, 2:4]

ggplot(cell.percent, aes(x = x, fill = rownames(cell.percent))) + 
  geom_bar()

barplot(as.matrix(cell.percent), col = turbo(10))
?barplot

library(tidyr)
cell.percent.long <- gather(cell.percent, Genotype, Percentage, WT1:E1020K, factor_key = TRUE)
ggplot(cell.percent.long, aes(fill=Population, y=Percentage, x=Genotype)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Population percentages across genotypes") +
  geom_text(aes(label = Percentage), position = position_stack(vjust = 0.5), size = 3) +
  theme_classic() +
  xlab("") +
  theme(plot.title = element_text(color="black", size=15))



##Find markers in each cluster
head(CD4_cells[[]])
Idents(CD4_cells) <- CD4_cells$wsnn_res.0.9
CD4_cells@meta.data$seurat_clusters <- CD4_cells$wsnn_res.0.9


CD4cells_RNA_c0 <- FindMarkers(CD4_cells, ident.1 = 0, assay = "RNA")
CD4cells_RNA_c0 <- tibble::rownames_to_column(CD4cells_RNA_c0, "Genes")
write_xlsx(CD4cells_RNA_c0, "CD4cells_RNA_c0.xlsx")
CD4cells_ADT_c0 <- FindMarkers(CD4_cells, ident.1 = 0, assay = "ADT")
CD4cells_ADT_c0 <- tibble::rownames_to_column(CD4cells_ADT_c0, "Genes")
write_xlsx(CD4cells_ADT_c0, "CD4cells_ADT_c0.xlsx")

CD4cells_RNA_c2 <- FindMarkers(CD4_cells, ident.1 = 2, assay = "RNA")
CD4cells_RNA_c2 <- tibble::rownames_to_column(CD4cells_RNA_c2, "Genes")
write_xlsx(CD4cells_RNA_c2, "CD4cells_RNA_c2.xlsx")
CD4cells_ADT_c2 <- FindMarkers(CD4_cells, ident.1 = 2, assay = "ADT")
CD4cells_ADT_c2 <- tibble::rownames_to_column(CD4cells_ADT_c2, "Genes")
write_xlsx(CD4cells_ADT_c2, "CD4cells_ADT_c2.xlsx")

CD4cells_RNA_c1 <- FindMarkers(CD4_cells, ident.1 = 1, assay = "RNA")
CD4cells_RNA_c1 <- tibble::rownames_to_column(CD4cells_RNA_c1, "Genes")
write_xlsx(CD4cells_RNA_c1, "CD4cells_RNA_c1.xlsx")
CD4cells_ADT_c1 <- FindMarkers(CD4_cells, ident.1 = 1, assay = "ADT")
CD4cells_ADT_c1 <- tibble::rownames_to_column(CD4cells_ADT_c1, "Genes")
write_xlsx(CD4cells_ADT_c1, "CD4cells_ADT_c1.xlsx")

CD4cells_RNA_c3 <- FindMarkers(CD4_cells, ident.1 = 3, assay = "RNA")
CD4cells_RNA_c3 <- tibble::rownames_to_column(CD4cells_RNA_c3, "Genes")
write_xlsx(CD4cells_RNA_c3, "CD4cells_RNA_c3.xlsx")
CD4cells_ADT_c3 <- FindMarkers(CD4_cells, ident.1 = 3, assay = "ADT")
CD4cells_ADT_c3 <- tibble::rownames_to_column(CD4cells_ADT_c3, "Genes")
write_xlsx(CD4cells_ADT_c3, "CD4cells_ADT_c3.xlsx")

CD4cells_RNA_c4 <- FindMarkers(CD4_cells, ident.1 = 4, assay = "RNA")
CD4cells_RNA_c4 <- tibble::rownames_to_column(CD4cells_RNA_c4, "Genes")
write_xlsx(CD4cells_RNA_c4, "CD4cells_RNA_c4.xlsx")
CD4cells_ADT_c4 <- FindMarkers(CD4_cells, ident.1 = 4, assay = "ADT")
CD4cells_ADT_c4 <- tibble::rownames_to_column(CD4cells_ADT_c4, "Genes")
write_xlsx(CD4cells_ADT_c4, "CD4cells_ADT_c4.xlsx")

CD4cells_RNA_c5 <- FindMarkers(CD4_cells, ident.1 = 5, assay = "RNA")
CD4cells_RNA_c5 <- tibble::rownames_to_column(CD4cells_RNA_c5, "Genes")
write_xlsx(CD4cells_RNA_c5, "CD4cells_RNA_c5.xlsx")
CD4cells_ADT_c5 <- FindMarkers(CD4_cells, ident.1 = 5, assay = "ADT")
CD4cells_ADT_c5 <- tibble::rownames_to_column(CD4cells_ADT_c5, "Genes")
write_xlsx(CD4cells_ADT_c5, "CD4cells_ADT_c5.xlsx")

CD4cells_RNA_c6 <- FindMarkers(CD4_cells, ident.1 = 6, assay = "RNA")
CD4cells_RNA_c6 <- tibble::rownames_to_column(CD4cells_RNA_c6, "Genes")
write_xlsx(CD4cells_RNA_c6, "CD4cells_RNA_c6.xlsx")
CD4cells_ADT_c6 <- FindMarkers(CD4_cells, ident.1 = 6, assay = "ADT")
CD4cells_ADT_c6 <- tibble::rownames_to_column(CD4cells_ADT_c6, "Genes")
write_xlsx(CD4cells_ADT_c6, "CD4cells_ADT_c6.xlsx")

CD4cells_RNA_c7 <- FindMarkers(CD4_cells, ident.1 = 7, assay = "RNA")
CD4cells_RNA_c7 <- tibble::rownames_to_column(CD4cells_RNA_c7, "Genes")
write_xlsx(CD4cells_RNA_c7, "CD4cells_RNA_c7.xlsx")
CD4cells_ADT_c7 <- FindMarkers(CD4_cells, ident.1 = 7, assay = "ADT")
CD4cells_ADT_c7 <- tibble::rownames_to_column(CD4cells_ADT_c7, "Genes")
write_xlsx(CD4cells_ADT_c7, "CD4cells_ADT_c7.xlsx")

CD4cells_RNA_c8 <- FindMarkers(CD4_cells, ident.1 = 8, assay = "RNA")
CD4cells_RNA_c8 <- tibble::rownames_to_column(CD4cells_RNA_c8, "Genes")
write_xlsx(CD4cells_RNA_c8, "CD4cells_RNA_c8.xlsx")
CD4cells_ADT_c8 <- FindMarkers(CD4_cells, ident.1 = 8, assay = "ADT")
CD4cells_ADT_c8 <- tibble::rownames_to_column(CD4cells_ADT_c8, "Genes")
write_xlsx(CD4cells_ADT_c8, "CD4cells_ADT_c8.xlsx")

CD4cells_RNA_c9 <- FindMarkers(CD4_cells, ident.1 = 9, assay = "RNA")
CD4cells_RNA_c9 <- tibble::rownames_to_column(CD4cells_RNA_c9, "Genes")
write_xlsx(CD4cells_RNA_c9, "CD4cells_RNA_c9.xlsx")
CD4cells_ADT_c9 <- FindMarkers(CD4_cells, ident.1 = 9, assay = "ADT")
CD4cells_ADT_c9 <- tibble::rownames_to_column(CD4cells_ADT_c9, "Genes")
write_xlsx(CD4cells_ADT_c9, "CD4cells_ADT_c9.xlsx")

##TFIDF
DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c0_genes <- WhichCells(object = CD4_cells, ident = "0")
CD4cells_TFIDF_c0 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c0_genes, colnames(CD4_cells))
CD4cells_TFIDF_c0 <- tibble::rownames_to_column(CD4cells_TFIDF_c0, "Genes")
write_xlsx(CD4cells_TFIDF_c0, "CD4cells_TFIDF_c0.xlsx")

DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c1_genes <- WhichCells(object = CD4_cells, ident = "1")
CD4cells_TFIDF_c1 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c1_genes, colnames(CD4_cells))
CD4cells_TFIDF_c1 <- tibble::rownames_to_column(CD4cells_TFIDF_c1, "Genes")
write_xlsx(CD4cells_TFIDF_c1, "CD4cells_TFIDF_c1.xlsx")

DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c2_genes <- WhichCells(object = CD4_cells, ident = "2")
CD4cells_TFIDF_c2 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c2_genes, colnames(CD4_cells))
CD4cells_TFIDF_c2 <- tibble::rownames_to_column(CD4cells_TFIDF_c2, "Genes")
write_xlsx(CD4cells_TFIDF_c2, "CD4cells_TFIDF_c2.xlsx")

DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c3_genes <- WhichCells(object = CD4_cells, ident = "3")
CD4cells_TFIDF_c3 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c3_genes, colnames(CD4_cells))
CD4cells_TFIDF_c3 <- tibble::rownames_to_column(CD4cells_TFIDF_c3, "Genes")
write_xlsx(CD4cells_TFIDF_c3, "CD4cells_TFIDF_c3.xlsx")

DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c4_genes <- WhichCells(object = CD4_cells, ident = "4")
CD4cells_TFIDF_c4 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c4_genes, colnames(CD4_cells))
CD4cells_TFIDF_c4 <- tibble::rownames_to_column(CD4cells_TFIDF_c4, "Genes")
write_xlsx(CD4cells_TFIDF_c4, "CD4cells_TFIDF_c4.xlsx")

DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c5_genes <- WhichCells(object = CD4_cells, ident = "5")
CD4cells_TFIDF_c5 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c5_genes, colnames(CD4_cells))
CD4cells_TFIDF_c5 <- tibble::rownames_to_column(CD4cells_TFIDF_c5, "Genes")
write_xlsx(CD4cells_TFIDF_c5, "CD4cells_TFIDF_c5.xlsx")

DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c6_genes <- WhichCells(object = CD4_cells, ident = "6")
CD4cells_TFIDF_c6 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c6_genes, colnames(CD4_cells))
CD4cells_TFIDF_c6 <- tibble::rownames_to_column(CD4cells_TFIDF_c6, "Genes")
write_xlsx(CD4cells_TFIDF_c6, "CD4cells_TFIDF_c6.xlsx")

DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c7_genes <- WhichCells(object = CD4_cells, ident = "7")
CD4cells_TFIDF_c7 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c7_genes, colnames(CD4_cells))
CD4cells_TFIDF_c7 <- tibble::rownames_to_column(CD4cells_TFIDF_c7, "Genes")
write_xlsx(CD4cells_TFIDF_c7, "CD4cells_TFIDF_c7.xlsx")

DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c8_genes <- WhichCells(object = CD4_cells, ident = "8")
CD4cells_TFIDF_c8 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c8_genes, colnames(CD4_cells))
CD4cells_TFIDF_c8 <- tibble::rownames_to_column(CD4cells_TFIDF_c8, "Genes")
write_xlsx(CD4cells_TFIDF_c8, "CD4cells_TFIDF_c8.xlsx")

DefaultAssay(CD4_cells)<-"RNA"
CD4cells_TFIDF_c9_genes <- WhichCells(object = CD4_cells, ident = "9")
CD4cells_TFIDF_c9 <- tfidf(GetAssayData(CD4_cells), CD4cells_TFIDF_c9_genes, colnames(CD4_cells))
CD4cells_TFIDF_c9 <- tibble::rownames_to_column(CD4cells_TFIDF_c9, "Genes")
write_xlsx(CD4cells_TFIDF_c9, "CD4cells_TFIDF_c9.xlsx")

###Individual markers
##Abs
DefaultAssay(CD4_cells)<-"ADT"
CD4_cells_p3
CD44 <- FeaturePlot(CD4_cells, features = "Cd44", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CD44") + theme(plot.title = element_text(color="black", size=15))
CD62L <- FeaturePlot(CD4_cells, features = "Cd62l", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CD62L")+ theme(plot.title = element_text(color="black", size=15))
PD_1 <- FeaturePlot(CD4_cells, features = "Pd-1", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("PD-1")+ theme(plot.title = element_text(color="black", size=15))
CXCR5 <- FeaturePlot(CD4_cells, features = "Cxcr5", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CXCR5")+ theme(plot.title = element_text(color="black", size=15))
CD25 <- FeaturePlot(CD4_cells, features = "Cd-25", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CD25")+ theme(plot.title = element_text(color="black", size=15))
GITR <- FeaturePlot(CD4_cells, features = "Gitr", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("GITR")+ theme(plot.title = element_text(color="black", size=15))
CD69 <- FeaturePlot(CD4_cells, features = "Cd69", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CD69")+ theme(plot.title = element_text(color="black", size=15))
PD_L1 <- FeaturePlot(CD4_cells, features = "Pd-L1", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("PD-L1")+ theme(plot.title = element_text(color="black", size=15))
PD_L2 <- FeaturePlot(CD4_cells, features = "Pd-L2", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("PD-L2")+ theme(plot.title = element_text(color="black", size=15))
CD86 <- FeaturePlot(CD4_cells, features = "Cd86", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CD86")+ theme(plot.title = element_text(color="black", size=15))
CD80 <- FeaturePlot(CD4_cells, features = "Cd80", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CD80")+ theme(plot.title = element_text(color="black", size=15))
CD38 <- FeaturePlot(CD4_cells, features = "Cd38", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CD38")+ theme(plot.title = element_text(color="black", size=15))
CD95 <- FeaturePlot(CD4_cells, features = "Cd95", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CD95")+ theme(plot.title = element_text(color="black", size=15))
ICOS <- FeaturePlot(CD4_cells, features = "Icos", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("ICOS")+ theme(plot.title = element_text(color="black", size=15))
CTLA4 <- FeaturePlot(CD4_cells, features = "Ctla4", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("CTLA4")+ theme(plot.title = element_text(color="black", size=15))
VlnPlot(CD4_cells, features = "Nt5e")


CD4_cells_3 <- DimPlot(CD4_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 2.5, label.box = TRUE) +  ggtitle("Seurat Clusters") + theme_void() + NoLegend() + theme(plot.title = element_text(color="black", size=15))
#Central and Effector Memory cells
Memory_cells <- grid.arrange(CD44, Cd44, CD62L, SELL, ncol = 2)
Central_memory_cells <- FeaturePlot(CD4_cells, features = c("Cd44", "Cd62l"), blend = TRUE, reduction = "wnn.umap", pt.size = 1.8) + theme_void()

#T follicular helper cells
Tfh_cells <- grid.arrange(CD4_cells_3, BCL6, PD_1, CXCR5, ncol = 2)
CD4_cells_3 <- DimPlot(CD4_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 2.5, label.box = TRUE) +  ggtitle("Seurat Clusters") + theme_void() + NoLegend() + theme(plot.title = element_text(color="black", size=15))

#Tregs
Tregs <- grid.arrange(CD4_cells_3, FOXP3, CD25, Il2ra, GITR, Tnfrsf18, ncol = 3)

CD44 + CD62L
grid.arrange(CD25, FOXP3, BCL6, CD4_cells_p3, GZMB, ncol = 2)
CD80 + CD86
CD69

CD25 + CD4_cells_p3

FeaturePlot(CD4_cells, features = c("Foxp3", "Bcl6"), blend = TRUE, reduction = "wnn.umap", pt.size = 1.8) + theme_void()

##Canonical RNA markers
DefaultAssay(CD4_cells)<-"RNA"
FOXP3 <- FeaturePlot(CD4_cells, features = "Foxp3", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("Foxp3")+ theme(plot.title = element_text(color="black", size=15))
BCL6 <- FeaturePlot(CD4_cells, features = "Bcl6", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("Bcl6")+ theme(plot.title = element_text(color="black", size=15))
Cd44 <- FeaturePlot(CD4_cells, features = "Cd44", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("Cd44")+ theme(plot.title = element_text(color="black", size=15))
SELL <- FeaturePlot(CD4_cells, features = "Sell", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("Sell")+ theme(plot.title = element_text(color="black", size=15))
Il2ra <- FeaturePlot(CD4_cells, features = "Il2ra", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("Il2ra")+ theme(plot.title = element_text(color="black", size=15))
Tnfrsf18 <- FeaturePlot(CD4_cells, features = "Tnfrsf18", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + ggtitle("Tnfrsf18")+ theme(plot.title = element_text(color="black", size=15))


##Cytotoxic markers
GZMB <- FeaturePlot(CD4_cells, features = "Gzmb", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
IFNG <- FeaturePlot(CD4_cells, features = "Ifng", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
TNFA <- FeaturePlot(CD4_cells, features = "Tnf", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
Il2 <- FeaturePlot(CD4_cells, features = "Il2", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))

grid.arrange(CD4_cells_3, Cd44, GZMB, GZMK, EOMES, NKG7, KLRK1, TBX21, ncol= 4)

#Other markers
SELL <- FeaturePlot(CD4_cells, features = "Sell", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
Cd44 <- FeaturePlot(CD4_cells, features = "Cd44", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))



CCL5 <- FeaturePlot(CD4_cells, features = "Ccl5", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
NKG7 <- FeaturePlot(CD4_cells, features = "Nkg7", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
HOPX <- FeaturePlot(CD4_cells, features = "Hopx", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))

#Cytotoxic markers Grace assessed with flow cytometry
DefaultAssay(CD4_cells)<-"RNA"

EOMES <- FeaturePlot(CD4_cells, features = "Eomes", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Cd44 <- FeaturePlot(CD4_cells, features = "Cd44", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Eomes <- VlnPlot(CD4_cells, features = "Eomes") + 
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  NoLegend()

GZMB <- FeaturePlot(CD4_cells, features = "Gzmb", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Gzmb <- VlnPlot(CD4_cells, features = "Gzmb") + 
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  NoLegend()

GZMK <- FeaturePlot(CD4_cells, features = "Gzmk", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Gzmk <- VlnPlot(CD4_cells, features = "Gzmk") + 
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  NoLegend()

NKG7 <- FeaturePlot(CD4_cells, features = "Nkg7", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Nkg7 <- VlnPlot(CD4_cells, features = "Nkg7") + 
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  NoLegend()

TBX21 <- FeaturePlot(CD4_cells, features = "Tbx21", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Tbx21 <- VlnPlot(CD4_cells, features = "Tbx21") + 
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  NoLegend()

KLRK1 <- FeaturePlot(CD4_cells, features = "Klrk1", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Klrk1 <- VlnPlot(CD4_cells, features = "Klrk1") + 
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  NoLegend()

DefaultAssay(CD4_cells)<-"ADT"
CD80 <- FeaturePlot(CD4_cells, features = "Cd80", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Cd80 <- VlnPlot(CD4_cells, features = "Cd80") + 
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  NoLegend()

CD86 <- FeaturePlot(CD4_cells, features = "Cd86", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Cd86 <- VlnPlot(CD4_cells, features = "Cd86") + 
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  NoLegend()

Th1_plot <- grid.arrange(EOMES, GZMB, GZMK, NKG7, TBX21, KLRK1,  ncol = 3)

Th1_plot_Vln <- VlnPlot(CD4_cells, features = CD4_Th1_genes, stack = TRUE, flip = TRUE, pt.size = 0.5) + 
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  NoLegend()

CD4_Th1_Dot <- DotPlot(CD4_cells, features = CD4_Th1_genes, cluster.idents = FALSE) + RotatedAxis() + scale_colour_viridis() + theme_classic()


grid.arrange(CD4_cells_p3, Th1_plot, Th1_plot_Vln, CD4_Th1_Dot, ncol = 2)
grid.arrange(CD4_cells_p3, Th1_plot, ncol = 2)

#Markers for cluster 3
IFNG <- FeaturePlot(CD4_cells, features = "Ifng", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

LAG3 <- FeaturePlot(CD4_cells, features = "Lag3", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

Il21 <- FeaturePlot(CD4_cells, features = "Il21", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

DefaultAssay(CD4_cells) <- "ADT"
PD_1 <- FeaturePlot(CD4_cells, features = "Pd-1", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

DefaultAssay(CD4_cells) <- "ADT"
PD_L1 <- FeaturePlot(CD4_cells, features = "Pd-L1", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

DefaultAssay(CD4_cells) <- "ADT"
CD38 <- FeaturePlot(CD4_cells, features = "Cd38", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

grid.arrange(CD4_cells_3, CD44, PD_1, PD_L1, CD38, IFNG, LAG3, Il21, ncol = 4)

#Markers for cluster 6
DefaultAssay(CD4_cells) <- "ADT"
GITR <- FeaturePlot(CD4_cells, features = "Gitr", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

DefaultAssay(CD4_cells) <- "ADT"
CD86 <- FeaturePlot(CD4_cells, features = "Cd86", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

DefaultAssay(CD4_cells) <- "ADT"
CD40 <- FeaturePlot(CD4_cells, features = "Cd40", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

FeaturePlot(CD4_cells, features = "Cd40", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

FeaturePlot(CD4_cells, features = "Il1rl1", reduction = "wnn.umap", cols = col_con, pt.size = 3) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

DefaultAssay(CD4_cells) <- "RNA"
CD81 <- FeaturePlot(CD4_cells, features = "Cd81", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

DefaultAssay(CD4_cells) <- "RNA"
CCR8 <- FeaturePlot(CD4_cells, features = "Ccr8", reduction = "wnn.umap", cols = col_con, pt.size = 0.8) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

DefaultAssay(CD4_cells) <- "RNA"
CD154 <- FeaturePlot(CD4_cells, features = "Cd40lg", reduction = "wnn.umap", cols = col_con, pt.size = 1.5) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

FeaturePlot(CD4_cells, features = "Il10ra", reduction = "wnn.umap", cols = col_con, pt.size = 1.5) +
  theme_light() +
  theme(plot.title = element_text(color="black", size=15)) +
  xlab("UMAP1") +
  ylab("UMAP2")

VlnPlot(CD4_cells, features = "Cd40lg")

grid.arrange(CD4_cells_3, PD_1, GITR, CD86, CD38, CD40, CD81, CCR8, ncol = 4)
##GSEA
Book1 <- read.csv("Book1.csv",header = T, sep = ',')

HALLMARK_INTERFERON_ALPHA_RESPONSE <- list(Book1$HALLMARK_INTERFERON_ALPHA_RESPONSE)
HALLMARK_INTERFERON_ALPHA_RESPONSE <- lapply(HALLMARK_INTERFERON_ALPHA_RESPONSE, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = HALLMARK_INTERFERON_ALPHA_RESPONSE, name = "HALLMARK_INTERFERON_ALPHA_RESPONSE", n = 22, search = FALSE)
VlnPlot(CD4_cells, c("HALLMARK_INTERFERON_ALPHA_RESPONSE1"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(CD4_cells, c("HALLMARK_INTERFERON_ALPHA_RESPONSE1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

Book2 <- read.csv(file.choose("CD4_Tcells_Signatures.csv"),header = T, sep = ',')


#Naive and Memory T cells
#Naive T cells
Tn <- list(Book2$Tn)
Tn <- lapply(Tn, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tn, name = "Tn", search = FALSE)
Cluster1 <-AddModuleScore(Cluster1, features = Tn, name = "Tn", search = FALSE)
VlnPlot(Cluster1, c("Tn1"), pt.size  = 0.5, cols = turbo(10))
VlnPlot(CD4_cells, c("Tn1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tn1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Central memory T cells
Tcm <- list(Book2$Tcm)
Tcm <- lapply(Tcm, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tcm, name = "Tcm", search = FALSE)
Cluster1 <-AddModuleScore(Cluster1, features = Tcm, name = "Tcm", search = FALSE)
VlnPlot(Cluster1, c("Tcm1"), pt.size  = 0.5, cols = turbo(10))
VlnPlot(CD4_cells, c("Tcm1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tcm1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Effector memory T cells
Tem <- list(Book2$Tem)
Tem <- lapply(Tem, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tem, name = "Tem", search = FALSE)
Cluster1 <-AddModuleScore(Cluster1, features = Tem, name = "Tem", search = FALSE)
VlnPlot(Cluster1, c("Tem1"), pt.size  = 0.5, cols = turbo(10))
VlnPlot(CD4_cells, c("Tem1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(Cluster1, c("Tem1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Effector memory T cells re-expressing CD45RA
Temra <- list(Book2$Temra)
Temra <- lapply(Temra, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Temra, name = "Temra", search = FALSE)
Cluster1 <-AddModuleScore(Cluster1, features = Temra, name = "Temra", search = FALSE)
VlnPlot(CD4_cells, c("Temra1"), pt.size  = 0.5, cols = turbo(10))
VlnPlot(Cluster1, c("Temra1"), pt.size  = 0.5, cols = turbo(5))
FeaturePlot(CD4_cells, c("Temra1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Natural Tregs
nTreg <- list(Book2$nTreg)
nTreg <- lapply(nTreg, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = nTreg, name = "nTreg", search = FALSE)
VlnPlot(CD4_cells, c("nTreg1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("nTreg1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Naive and memory T cells - heatmap
Data<-CD4_cells@meta.data
Forest <-Data %>% 
  group_by(`Genotype`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$Genotype
Forest<-Forest[ , c(56:60)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)


Data<-CD4_cells@meta.data
Forest <-Data %>% 
  group_by(`seurat_clusters`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$seurat_clusters
Forest<-Forest[ , c(56:60)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)


#Naive or memory T cells stimulated with Th-subset-associated cytokines
#Th1 cells
Tn_Th1 <- list(Book2$Tn_Th1)
Tn_Th1 <- lapply(Tn_Th1, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tn_Th1, name = "Tn_Th1", search = FALSE)
VlnPlot(CD4_cells, c("Tn_Th11"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tn_Th11"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

Tm_Th1 <- list(Book2$Tm_Th1)
Tm_Th1 <- lapply(Tm_Th1, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tm_Th1, name = "Tm_Th1", search = FALSE)
VlnPlot(CD4_cells, c("Tm_Th11"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tm_Th11"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Th2
Tn_Th2 <- list(Book2$Tn_Th2)
Tn_Th2 <- lapply(Tn_Th2, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tn_Th2, name = "Tn_Th2", search = FALSE)
VlnPlot(CD4_cells, c("Tn_Th21"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tn_Th21"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Th17 and iTreg
Tn_Th17 <- list(Book2$Tn_Th17)
Tn_Th17 <- lapply(Tn_Th17, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tn_Th17, name = "Tn_Th17", search = FALSE)
Cluster1 <-AddModuleScore(Cluster1, features = Tn_Th17, name = "Tn_Th17", search = FALSE)
VlnPlot(Cluster1, c("Tn_Th171"), pt.size  = 0.5, cols = turbo(10))
VlnPlot(CD4_cells, c("Tn_Th171"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tn_Th171"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

Tn_iTreg <- list(Book2$Tn_iTreg)
Tn_iTreg <- lapply(Tn_iTreg, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tn_iTreg, name = "Tn_iTreg", search = FALSE)
VlnPlot(CD4_cells, c("Tn_iTreg1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tn_iTreg1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

Tm_Th17_iTreg <- list(Book2$Tm_Th17_iTreg)
Tm_Th17_iTreg <- lapply(Tm_Th17_iTreg, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tm_Th17_iTreg, name = "Tm_Th17_iTreg", search = FALSE)
Cluster1 <-AddModuleScore(Cluster1, features = Tm_Th17_iTreg, name = "Tm_Th17_iTreg", search = FALSE)
VlnPlot(Cluster1, c("Tm_Th17_iTreg1"), pt.size  = 0.5, cols = turbo(10))
VlnPlot(CD4_cells, c("Tm_Th17_iTreg1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tm_Th17_iTreg1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#IFNbeta response
IFNbeta_resp <- list(Book2$IFNbeta_resp)
IFNbeta_resp <- lapply(IFNbeta_resp, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = IFNbeta_resp, name = "IFNbeta_resp", search = FALSE)
VlnPlot(CD4_cells, c("IFNbeta_resp1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("IFNbeta_resp1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

Tm_IFNbeta_resp <- list(Book2$Tm_IFNbeta_resp)
Tm_IFNbeta_resp <- lapply(Tm_IFNbeta_resp, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tm_IFNbeta_resp, name = "Tm_IFNbeta_resp", search = FALSE)
VlnPlot(CD4_cells, c("Tm_IFNbeta_resp1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tm_IFNbeta_resp1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Naive and memory T cell stumuated with Th-assocaited cytokines - heatmap
Data<-CD4_cells@meta.data
Forest <-Data %>% 
  group_by(`Genotype`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$Genotype
Forest<-Forest[ , c(48:55)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)


Data<-CD4_cells@meta.data
Forest <-Data %>% 
  group_by(`seurat_clusters`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$seurat_clusters
Forest<-Forest[ , c(48:55)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)

#Summary naive and memory T cells stimulated with Th-subsets
#Activated (TCR and CD28 stimulation) Temra
TEMRA_Th0 <- list(Book2$TEMRA_Th0)
TEMRA_Th0 <- lapply(TEMRA_Th0, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TEMRA_Th0, name = "TEMRA_Th0", search = FALSE)
VlnPlot(CD4_cells, c("TEMRA_Th01"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TEMRA_Th01"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#iTreg/Th17-stimulated (cytokines inducing iTreg/Th17 differentiation) Temra cells
TEMRA_iTreg_Th17 <- list(Book2$TEMRA_iTreg_Th17)
TEMRA_iTreg_Th17 <- lapply(TEMRA_iTreg_Th17, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TEMRA_iTreg_Th17, name = "TEMRA_iTreg_Th17", search = TRUE)
VlnPlot(CD4_cells, c("TEMRA_iTreg_Th171"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TEMRA_iTreg_Th171"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Activated Tem cells
TEM_Th0 <- list(Book2$TEM_Th0)
TEM_Th0 <- lapply(TEM_Th0, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TEM_Th0, name = "TEM_Th0", search = FALSE)
VlnPlot(CD4_cells, c("TEM_Th01"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TEM_Th01"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#iTreg/Th17-stimulated Tem cells
TEM_iTreg_Th17 <- list(Book2$TEM_iTreg_Th17)
TEM_iTreg_Th17 <- lapply(TEM_iTreg_Th17, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TEM_iTreg_Th17, name = "TEM_iTreg_Th17", search = FALSE)
VlnPlot(CD4_cells, c("TEM_iTreg_Th171"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TEM_iTreg_Th171"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Activated Tcm1 cells
TCM1_Th0 <- list(Book2$TCM1_Th0)
TCM1_Th0 <- lapply(TCM1_Th0, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TCM1_Th0, name = "TCM1_Th0", search = FALSE)
VlnPlot(CD4_cells, c("TCM1_Th01"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TCM1_Th01"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#iTreg/Th17-stimulated Tcm1 cells
TCM1_iTreg_Th17 <- list(Book2$TCM1_iTreg_Th17)
TCM1_iTreg_Th17 <- lapply(TCM1_iTreg_Th17, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TCM1_iTreg_Th17, name = "TCM1_iTreg_Th17", search = FALSE)
VlnPlot(CD4_cells, c("TCM1_iTreg_Th171"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TCM1_iTreg_Th171"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Activated Tcm2 cells
TCM2_Th0 <- list(Book2$TCM2_Th0)
TCM2_Th0 <- lapply(TCM2_Th0, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TCM2_Th0, name = "TCM2_Th0", search = FALSE)
VlnPlot(CD4_cells, c("TCM2_Th01"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TCM2_Th01"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#iTreg/Th17-stimulated Tcm2 cells
TCM2_iTreg_Th17 <- list(Book2$TCM2_iTreg_Th17)
TCM2_iTreg_Th17 <- lapply(TCM2_iTreg_Th17, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TCM2_iTreg_Th17, name = "TCM2_iTreg_Th17", search = FALSE)
VlnPlot(CD4_cells, c("TCM2_iTreg_Th171"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TCM2_iTreg_Th171"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#iTregstimulated naive T cells
TN_iTreg <- list(Book2$TN_iTreg)
TN_iTreg <- lapply(TN_iTreg, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TN_iTreg, name = "TN_iTreg", search = FALSE)
VlnPlot(CD4_cells, c("TN_iTreg1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TN_iTreg1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Activated naive T cells
TN_Th0 <- list(Book2$TN_Th0)
TN_Th0 <- lapply(TN_Th0, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TN_Th0, name = "TN_Th0", search = FALSE)
VlnPlot(CD4_cells, c("TN_Th01"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TN_Th01"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Th2-stimulated naive T cells
TN_Th2 <- list(Book2$TN_Th2)
TN_Th2 <- lapply(TN_Th2, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TN_Th2, name = "TN_Th2", search = FALSE)
VlnPlot(CD4_cells, c("TN_Th21"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TN_Th21"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Th17-stimulated T cells
TN_Th17 <- list(Book2$TN_Th17)
TN_Th17 <- lapply(TN_Th17, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TN_Th17, name = "TN_Th17", search = FALSE)
VlnPlot(CD4_cells, c("TN_Th171"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TN_Th171"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#iTreg/Th17-stimulated T cells
TN_iTreg.Th17 <- list(Book2$TN_iTreg.Th17)
TN_iTreg.Th17 <- lapply(TN_iTreg.Th17, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TN_iTreg.Th17, name = "TN_iTreg.Th17", search = FALSE)
VlnPlot(CD4_cells, c("TN_iTreg.Th171"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TN_iTreg.Th171"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Natural Treg cells
nTregs_combined <- list(Book2$nTregs_combined)
nTregs_combined <- lapply(nTregs_combined, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = nTregs_combined, name = "nTregs_combined", search = FALSE)
VlnPlot(CD4_cells, c("nTregs_combined1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("nTregs_combined1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Summary T cell subsets - heatmap
Data<-CD4_cells@meta.data
Forest <-Data %>% 
  group_by(`Genotype`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$Genotype
Forest<-Forest[ , c(62:76)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)


Data<-CD4_cells@meta.data
Forest <-Data %>% 
  group_by(`seurat_clusters`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$seurat_clusters
Forest<-Forest[ , c(62:76)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)



#High for IFN-response genes
IFN_responsive_CD4Tcells <- list(Book2$IFN_responsive_CD4Tcells)
IFN_responsive_CD4Tcells <- lapply(IFN_responsive_CD4Tcells, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = IFN_responsive_CD4Tcells, name = "IFN_responsive_CD4Tcells", search = FALSE)
VlnPlot(CD4_cells, c("IFN_responsive_CD4Tcells1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("IFN_responsive_CD4Tcells1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#High mitotic genes
Mitotic_CD4Tcells <- list(Book2$Mitotic_CD4Tcells)
Mitotic_CD4Tcells <- lapply(Mitotic_CD4Tcells, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Mitotic_CD4Tcells, name = "Mitotic_CD4Tcells", search = FALSE)
VlnPlot(CD4_cells, c("Mitotic_CD4Tcells1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Mitotic_CD4Tcells1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#High for heatshock proteins
HeatShockhi_CD4Tcells <- list(Book2$HeatShockhi_CD4Tcells)
HeatShockhi_CD4Tcells <- lapply(HeatShockhi_CD4Tcells, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = HeatShockhi_CD4Tcells, name = "HeatShockhi_CD4Tcells", search = FALSE)
VlnPlot(CD4_cells, c("HeatShockhi_CD4Tcells1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("HeatShockhi_CD4Tcells1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#High for effector T cells
High_effector_CD4Tcells <- list(Book2$High_effector_CD4Tcells)
High_effector_CD4Tcells <- lapply(High_effector_CD4Tcells, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = High_effector_CD4Tcells, name = "High_effector_CD4Tcells", search = FALSE)
VlnPlot(CD4_cells, c("High_effector_CD4Tcells1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("High_effector_CD4Tcells1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)


#"Healthy" and intra-tumoural CD4 T cells
#CD4 Tcm cells
CD4cm <- list(Book2$CD4cm)
CD4cm <- lapply(CD4cm, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4cm, name = "CD4cm", search = FALSE)
VlnPlot(CD4_cells, c("CD4cm1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4cm1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#IL2RA-hi (CD25-hi) CD4 T cells
CD4IL2RAHI <- list(Book2$CD4IL2RAHI)
CD4IL2RAHI <- lapply(CD4IL2RAHI, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4IL2RAHI, name = "CD4IL2RAHI", search = FALSE)
VlnPlot(CD4_cells, c("CD4IL2RAHI1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4IL2RAHI1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#IL2RA-lo (CD25-lo) CD4 T cells
CD4IL2RALO <- list(Book2$CD4IL2RALO)
CD4IL2RALO <- lapply(CD4IL2RALO, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4IL2RALO, name = "CD4IL2RALO", search = FALSE)
VlnPlot(CD4_cells, c("CD4IL2RALO1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4IL2RALO1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Activated CD4 T cells
CD4ACTIVATED <- list(Book2$CD4ACTIVATED)
CD4ACTIVATED <- lapply(CD4ACTIVATED, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4ACTIVATED, name = "CD4ACTIVATED", search = FALSE)
VlnPlot(CD4_cells, c("CD4ACTIVATED1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4ACTIVATED1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#GZMB CD4 T cells
CD4GZMB <- list(Book2$CD4GZMB)
CD4GZMB <- lapply(CD4GZMB, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4GZMB, name = "CD4GZMB", search = FALSE)
VlnPlot(CD4_cells, c("CD4GZMB1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4GZMB1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#GZMK CD4 T cells
CD4GZMK <- list(Book2$CD4GZMK)
CD4GZMK <- lapply(CD4GZMK, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4GZMK, name = "CD4GZMK", search = FALSE)
VlnPlot(CD4_cells, c("CD4GZMK1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4GZMK1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Th17 CD4 T cells
CD4TH17 <- list(Book2$CD4TH17)
CD4TH17 <- lapply(CD4TH17, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4TH17, name = "CD4TH17", search = FALSE)
VlnPlot(CD4_cells, c("CD4TH171"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4TH171"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#CXCL13 CD4 T cells
CD4CXCL13 <- list(Book2$CD4CXCL13)
CD4CXCL13 <- lapply(CD4CXCL13, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4CXCL13, name = "CD4CXCL13", search = FALSE)
VlnPlot(CD4_cells, c("CD4CXCL131"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4CXCL131"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#High MT-genes CD4 T cells
CD4MITO <- list(Book2$CD4MITO)
CD4MITO <- lapply(CD4MITO, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4MITO, name = "CD4MITO", search = FALSE)
VlnPlot(CD4_cells, c("CD4MITO1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4MITO1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#High HSP CD4 T cells
CD4HSP <- list(Book2$CD4HSP)
CD4HSP <- lapply(CD4HSP, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4HSP, name = "CD4HSP", search = FALSE)
VlnPlot(CD4_cells, c("CD4HSP1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4HSP1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Proliferating CD4 T cells
CD4PROLIF <- list(Book2$CD4PROLIF)
CD4PROLIF <- lapply(CD4PROLIF, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4PROLIF, name = "CD4PROLIF", search = FALSE)
VlnPlot(CD4_cells, c("CD4PROLIF1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4PROLIF1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Healthy and tumour CD4 T cell - heatmap
Data<-CD4_cells@meta.data
Forest <-Data %>% 
  group_by(`Genotype`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$Genotype
Forest<-Forest[ , c(80:90)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)


Data<-CD4_cells@meta.data
Forest <-Data %>% 
  group_by(`seurat_clusters`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$seurat_clusters
Forest<-Forest[ , c(80:90)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)



#IL2RA-hi CD4 T cells - tumour vs healthy
CD4IL2RAHI.tumor.v.normal <- list(Book2$CD4IL2RAHI.tumor.v.normal)
CD4IL2RAHI.tumor.v.normal <- lapply(CD4IL2RAHI.tumor.v.normal, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4IL2RAHI.tumor.v.normal, name = "CD4IL2RAHI.tumor.v.normal", search = FALSE)
VlnPlot(CD4_cells, c("CD4IL2RAHI.tumor.v.normal1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4IL2RAHI.tumor.v.normal1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#IL2RA-lo CD4 T cells - tumour vs healthy
CD4IL2RALO.tumor.v.normal <- list(Book2$CD4IL2RALO.tumor.v.normal)
CD4IL2RALO.tumor.v.normal <- lapply(CD4IL2RAHI.tumor.v.normal, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4IL2RALO.tumor.v.normal, name = "CD4IL2RALO.tumor.v.normal", search = FALSE)
VlnPlot(CD4_cells, c("CD4IL2RALO.tumor.v.normal1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4IL2RALO.tumor.v.normal1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#GZMB CD4 T cells - tumour vs healthy
CD4GZMB.tumor.v.normal <- list(Book2$CD4GZMB.tumor.v.normal)
CD4GZMB.tumor.v.normal <- lapply(CD4GZMB.tumor.v.normal, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4GZMB.tumor.v.normal, name = "CD4GZMB.tumor.v.normal", search = FALSE)
VlnPlot(CD4_cells, c("CD4GZMB.tumor.v.normal1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4GZMB.tumor.v.normal1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#GZMK CD4 T cells - tumour vs healthy
CD4GZMK.tumor.v.normal <- list(Book2$CD4GZMK.tumor.v.normal)
CD4GZMK.tumor.v.normal <- lapply(CD4GZMK.tumor.v.normal, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = CD4GZMK.tumor.v.normal, name = "CD4GZMK.tumor.v.normal", search = FALSE)
VlnPlot(CD4_cells, c("CD4GZMK.tumor.v.normal1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("CD4GZMK.tumor.v.normal1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Cytotoxic (TC2) and non-cytotoxic (TC1) T cells from supercentenarians
#Non-cytotoxic CD4 T cells
TC_1 <- list(Book2$TC_1)
TC_1 <- lapply(TC_1, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TC_1, name = "TC_1", search = FALSE)
VlnPlot(CD4_cells, c("TC_11"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("TC_11"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Cytotoxic CD4 T cells
TC_2 <- list(Book2$TC_2)
TC_2 <- lapply(TC_2, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = TC_2, name = "TC_2", search = FALSE)
plot40 <- VlnPlot(CD4_cells, c("TC_21"), pt.size  = 0.5) + ggtitle("TC2") + xlab("") + NoLegend()
plot41 <- FeaturePlot(CD4_cells, c("TC_21"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("TC2")

grid.arrange(plot41, plot40, ncol = 2)

#IL10+ CD44-hi CD4+ T cells
#Th1-like
Th1_sig <- list(Book2$Th1_sig)
Th1_sig <- lapply(Th1_sig, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Th1_sig, name = "Th1_sig", search = FALSE)
VlnPlot(CD4_cells, c("Th1_sig1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Th1_sig1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Tfh-like
Tfh_sig <- list(Book2$Tfh_sig)
Tfh_sig <- lapply(Tfh_sig, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Tfh_sig, name = "Tfh_sig", search = FALSE)
VlnPlot(CD4_cells, c("Tfh_sig1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Tfh_sig1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)


#Looking for C1 clusters upon all CD4 cells
C1_cluster0 <- list(Book2$C1_Cluster0)
C1_cluster0  <- lapply(C1_cluster0, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = C1_cluster0 , name = "C1_cluster0", search = FALSE)
VlnPlot(CD4_cells, c("C1_cluster01"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("C1_cluster01"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

C1_cluster3 <- list(Book2$Cluster3)
C1_cluster3  <- lapply(C1_cluster3, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = C1_cluster3 , name = "C1_cluster3", search = FALSE)
VlnPlot(CD4_cells, c("C1_cluster31"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("C1_cluster31"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

C1_cluster4 <- list(Book2$C1_Cluster4)
C1_cluster4  <- lapply(C1_cluster4, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = C1_cluster4 , name = "C1_cluster4", search = FALSE)
VlnPlot(CD4_cells, c("C1_cluster41"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("C1_cluster41"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)



Book2 <- read.csv(file.choose("CD4_Tcells_Signatures.csv"),header = T, sep = ',')
CD4_cells_p3


grid.arrange(CD44, Cd44, CD62L, SELL, CD4_cells_p3, ncol = 2)
FeaturePlot(CD4_cells, features = "Ncr1", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))


CD4_cells_p3
CD4_cells_Mouse

##Comparing individual cluster with each other
CD25_hi_abs <- FindMarkers(CD4_cells, ident.1 = 5, ident.2 = 4, assay = "ADT")
CD25_hi_rna  <- FindMarkers(CD4_cells, ident.1 = 5, ident.2 = 4, assay = "RNA")
CD25_lo_rna  <- FindMarkers(CD4_cells, ident.1 = 4, ident.2 = 5, assay = "RNA")

Cluster3_vs_Cluster1_ADT <- FindMarkers(CD4_cells, ident.1 = 3, ident.2 = 1, assay = "ADT")
Cluster3_vs_Cluster1_RNA <- FindMarkers(CD4_cells, ident.1 = 3, ident.2 = 1, assay = "RNA")

backup <- CD4_cells


##Looking for idnvidual markers
DefaultAssay(CD4_cells) <- "RNA"
FeaturePlot(CD4_cells, features = "Ifng", reduction = "wnn.umap", cols = col_con, pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
VlnPlot(CD4_cells, features = "Gzmb")

#GSEA analysis via GO and KEGG
Idents(CD4_cells) <- CD4_cells$seurat_clusters
CD4.markers <- FindAllMarkers(CD4_cells, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
CD4.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top100_all <- CD4.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% arrange(desc(avg_log2FC))
top100_pval <- subset(top100_all, rowSums(top100_all[5] < 0.05) > 0)
dim(top100_pval)
df <- top100_pval[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)

dfsample$`Cluster 0` = bitr(dfsample$`Cluster 0`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 1` = bitr(dfsample$`Cluster 1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 2` = bitr(dfsample$`Cluster 2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 3` = bitr(dfsample$`Cluster 3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 4` = bitr(dfsample$`Cluster 4`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 5` = bitr(dfsample$`Cluster 5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 6` = bitr(dfsample$`Cluster 6`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 7` = bitr(dfsample$`Cluster 7`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 8` = bitr(dfsample$`Cluster 8`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 9` = bitr(dfsample$`Cluster 9`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)

genelist <- list("0" = dfsample$`Cluster 0`$ENTREZID,"1" = dfsample$`Cluster 1`$ENTREZID, "3" = dfsample$`Cluster 3`$ENTREZID, "4" = dfsample$`Cluster 4`$ENTREZID,
                 "5" = dfsample$`Cluster 5`$ENTREZID, "6" = dfsample$`Cluster 6`$ENTREZID, "7" = dfsample$`Cluster 7`$ENTREZID, "8" = dfsample$`Cluster 8`$ENTREZID, "9" = dfsample$`Cluster 9`$ENTREZID)

GOclusterplot <- compareCluster(geneClusters = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
dotplot(GOclusterplot, font.size = 8, label_format = 100)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", organism = "mmu")
dotplot(KEGGclusterplot,  font.size = 8, label_format = 100)

##Subset Cluster 1
Idents(CD4_cells)
Cluster1 <- subset(CD4_cells, idents = 1)
DefaultAssay(Cluster1) <- "RNA" #For log normalisation
Cluster1 <- NormalizeData(Cluster1, verbose = TRUE)
Cluster1 <- FindVariableFeatures(Cluster1, nfeatures = 3000)
Cluster1 <- ScaleData(Cluster1)
top20_plot <-  head(VariableFeatures(Cluster1), 20)
plot2.1 <-  VariableFeaturePlot(Cluster1)
top20_plot_plot <-  LabelPoints(plot = plot2.1, points = top20_plot, repel = TRUE, xnudge = 0, ynudge = 0) #CCL5 as a diffining marker?
Cluster1 <- RunPCA(Cluster1, verbose = FALSE, features = VariableFeatures(object = Cluster1))
pca_variance <- Cluster1@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #30
DefaultAssay(Cluster1) <- "RNA" #For log normalisation
Cluster1 <- FindNeighbors(Cluster1, dims = 1:30)
Cluster1 <- FindClusters(Cluster1, resolution = 1, verbose = FALSE) #0.6 for the resolution
clustree(Cluster1, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
Cluster1 <-RunUMAP(Cluster1, dims = 1:30, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
Cluster1_p1 <- DimPlot(Cluster1, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("RNA Clustering") + theme_void() + NoLegend()
Cluster1_p1 <- Cluster1_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))
DefaultAssay(Cluster1) <- "ADT"
VariableFeatures(Cluster1) <- rownames(Cluster1[["ADT"]])
Cluster1 <- NormalizeData(Cluster1, normalization.method = "CLR", margin = 2)
Cluster1 <- ScaleData(Cluster1)
Cluster1 <- RunPCA(Cluster1, reduction.name = 'apca', approx = FALSE)
apca_variance <- Cluster1@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #26
Cluster1 <- FindNeighbors(Cluster1, dims = 1:26, reduction = "apca")
Cluster1 <- FindClusters(Cluster1, resolution = 1.0, verbose = FALSE) #1.0 for the resolution
clustree(Cluster1, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
Cluster1 <- RunUMAP(Cluster1, reduction = 'apca', dims = 1:26, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
Cluster1_p2 <- DimPlot(Cluster1, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("ADT Clustering") + theme_void() + NoLegend()
Cluster1_p2 <- Cluster1_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))
DefaultAssay(Cluster1) <- "RNA" #For log normalisation
Cluster1 <- FindMultiModalNeighbors(
  Cluster1, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:26), modality.weight.name = "RNA.weight")
Cluster1 <- FindClusters(Cluster1, graph.name = "wsnn", algorithm = 3, resolution = 1.0, verbose = TRUE) #0.9 for the resolution
clustree(Cluster1, prefix = "wsnn_res.") + theme(legend.position="bottom")
Cluster1 <- RunUMAP(Cluster1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Cluster1_p3 <- DimPlot(Cluster1, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Seurat Clusters") + theme_classic() + NoLegend()
Cluster1_p3 <- Cluster1_p3 + theme(plot.title = element_text(color="black", size=15)) + xlab("UMAP1") + ylab("UMAP2")


##Analysis of Cluster 1
#Cell numbers and percent per cluster

cell.numbers.2 <- table(Cluster1@meta.data$seurat_clusters, Cluster1@meta.data$Genotype)
cell.numbers.2 <- as.data.frame.matrix(cell.numbers.2)

cell_number_heatmap_2 <- pheatmap::pheatmap(t(cell.numbers.2), cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                                            cellwidth = 30,cellheight = 30, angle_col = 45, col = col_con)

Cluster1_meta <- Cluster1@meta.data
genotype.numbers.2 <- Cluster1_meta %>% dplyr::count(Genotype)
genotype.numbers.vector.2 <- genotype.numbers.2 %>% pull(n)
str(genotype.numbers.vector.2)
cell.percent.2 <- sweep(cell.numbers.2, 2, genotype.numbers.vector.2, "/")
cell.percent.2 <- cell.percent.2*100

cell_percent_heatmap_2 <- pheatmap::pheatmap(t(cell.percent.2), cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                                             cellwidth = 30,cellheight = 30, angle_col = 45, col = col_con)

#Cluster-associated markers
Cluster1_RNA_c0 <- FindMarkers(Cluster1, ident.1 = 0, assay = "RNA")
Cluster1_RNA_c1 <- FindMarkers(Cluster1, ident.1 = 1, assay = "RNA")
Cluster1_RNA_c2 <- FindMarkers(Cluster1, ident.1 = 2, assay = "RNA")
Cluster1_RNA_c3 <- FindMarkers(Cluster1, ident.1 = 3, assay = "RNA")
Cluster1_RNA_c4 <- FindMarkers(Cluster1, ident.1 = 4, assay = "RNA")

Cluster1_ADT_c0 <- FindMarkers(Cluster1, ident.1 = 0, assay = "ADT")
Cluster1_ADT_c1 <- FindMarkers(Cluster1, ident.1 = 1, assay = "ADT")
Cluster1_ADT_c2 <- FindMarkers(Cluster1, ident.1 = 2, assay = "ADT")
Cluster1_ADT_c3 <- FindMarkers(Cluster1, ident.1 = 3, assay = "ADT")
Cluster1_ADT_c4 <- FindMarkers(Cluster1, ident.1 = 4, assay = "ADT")

DefaultAssay(Cluster1)<-"RNA"
Idents(Cluster1)
Cluster1_TFIDF_c0_genes <- WhichCells(object = Cluster1, ident = "0")
Cluster1_TFIDF_c0 <- tfidf(GetAssayData(Cluster1), Cluster1_TFIDF_c0_genes, colnames(Cluster1))

Cluster1_TFIDF_c1_genes <- WhichCells(object = Cluster1, ident = "1")
Cluster1_TFIDF_c1 <- tfidf(GetAssayData(Cluster1), Cluster1_TFIDF_c1_genes, colnames(Cluster1))

Cluster1_TFIDF_c2_genes <- WhichCells(object = Cluster1, ident = "2")
Cluster1_TFIDF_c2 <- tfidf(GetAssayData(Cluster1), Cluster1_TFIDF_c2_genes, colnames(Cluster1))

Cluster1_TFIDF_c3_genes <- WhichCells(object = Cluster1, ident = "3")
Cluster1_TFIDF_c3 <- tfidf(GetAssayData(Cluster1), Cluster1_TFIDF_c3_genes, colnames(Cluster1))

Cluster1_TFIDF_c4_genes <- WhichCells(object = Cluster1, ident = "4")
Cluster1_TFIDF_c4 <- tfidf(GetAssayData(Cluster1), Cluster1_TFIDF_c4_genes, colnames(Cluster1))

#No obvious indication that the E1020K animal only has Temra cells
Cluster1_E1020K <- subset(Cluster1, subset = Genotype == "E1020K")
Cluster1_WT <- subset(Cluster1, subset = Genotype != "E1020K")

FeaturePlot(Cluster1, features = c("Cd44", "Ptprc"), blend = TRUE, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
FeaturePlot(Cluster1_E1020K, features = c("Cd44", "Ptprc"), blend = TRUE, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
FeaturePlot(Cluster1_WT, features = c("Cd44", "Ptprc"), blend = TRUE, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))

#Loking at idnividual markers
FeaturePlot(Cluster1, features = "Bap1", cols = col_con, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
VlnPlot(CD4_cells, features = "Gzmb")
DimPlot(Cluster1, reduction = "wnn.umap", split.by = "Genotype") + Cluster1_p3

#Cluster identification
#Cluster 0 - CD4 Tcm cells
Cluster1_CD44_CD62L + Cluster1_p3

#Cluster 1 - high for CD86, CD44, CD62L-lo - CD86-hi Tem cells
FeaturePlot(Cluster1, features = "Klhl6", cols = col_con, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
VlnPlot(Cluster1, features = "Klhl6")

#Cluster 2 - 
FeaturePlot(Cluster1, features = "Klhl6", cols = col_con, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
VlnPlot(Cluster1, features = "Klhl6")
C1_c2_vs_c1_abs <- FindMarkers(Cluster1, ident.1 = 2, ident.2 = 1, assay = "ADT")
C1_c2_vs_c1 <- FindMarkers(Cluster1, ident.1 = 2, ident.2 = 1, assay = "RNA")

#Cluster 3 - cytotoxic cells (Th1 cells)
Cluster1 <-AddModuleScore(Cluster1, features = TC_2, name = "TC_2", search = FALSE)
VlnPlot(Cluster1, c("TC_21"), pt.size  = 0.5, cols = turbo(10))
Cluster1_TC2 <- FeaturePlot(Cluster1, c("TC_21"), cols=col_con, reduction = "wnn.umap", pt.size  = 3)
grid.arrange(Cluster1_TC2, Cluster1_p3, ncol = 2)

Cluster1 <-AddModuleScore(Cluster1, features = Tn_Th1, name = "Tn_Th1", search = FALSE)
VlnPlot(Cluster1, c("Tn_Th11"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(Cluster1, c("Tn_Th11"), cols=col_con, reduction = "wnn.umap", pt.size  = 3)

#Th1 cells are low or high for CD44 expression (some of them are more memory-like), but high for CD38
DefaultAssay(Cluster1) <- "ADT"
FeaturePlot(Cluster1, features = c("Cd44", "Gzmk"), blend = TRUE, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
FeaturePlot(Cluster1, features = c("Cd44", "Eomes"), blend = TRUE, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
FeaturePlot(Cluster1, features = c("Cd44", "Nkg7"), blend = TRUE, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))

DefaultAssay(CD4_cells) <- "RNA"
FeaturePlot(CD4_cells, features = "Il21", cols = col_con, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
CD4_cells_p3 + BCL6

#Cluster 4 - Ccr2 and Ccr10 expresing cells

CD4_cells_p3
##Subset Cluster 3 - CD4_cells
Idents(CD4_cells)
Cluster3 <- subset(CD4_cells, idents = "Cluster 3")
DefaultAssay(Cluster3) <- "RNA" #For log normalisation
Cluster3 <- NormalizeData(Cluster3, verbose = TRUE)
Cluster3 <- FindVariableFeatures(Cluster3, nfeatures = 3000)
Cluster3 <- ScaleData(Cluster3)
top20_plot <-  head(VariableFeatures(Cluster3), 20)
plot2.1 <-  VariableFeaturePlot(Cluster3)
top20_plot_plot <-  LabelPoints(plot = plot2.1, points = top20_plot, repel = TRUE, xnudge = 0, ynudge = 0) #CCL5 as a diffining marker?
Cluster3 <- RunPCA(Cluster3, verbose = FALSE, features = VariableFeatures(object = Cluster3))
pca_variance <- Cluster3@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #30
DefaultAssay(Cluster3) <- "RNA" #For log normalisation
Cluster3 <- FindNeighbors(Cluster3, dims = 1:30)
Cluster3 <- FindClusters(Cluster3, resolution = 1, verbose = FALSE) #0.6 for the resolution
clustree(Cluster3, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
Cluster3 <-RunUMAP(Cluster3, dims = 1:30, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
Cluster3_p1 <- DimPlot(Cluster3, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("RNA Clustering") + theme_void() + NoLegend()
Cluster3_p1 <- Cluster3_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))
DefaultAssay(Cluster3) <- "ADT"
VariableFeatures(Cluster3) <- rownames(Cluster3[["ADT"]])
Cluster3 <- NormalizeData(Cluster3, normalization.method = "CLR", margin = 2)
Cluster3 <- ScaleData(Cluster3)
Cluster3 <- RunPCA(Cluster3, reduction.name = 'apca', approx = FALSE)
apca_variance <- Cluster3@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #22
Cluster3 <- FindNeighbors(Cluster3, dims = 1:22, reduction = "apca")
Cluster3 <- FindClusters(Cluster3, resolution = 1.0, verbose = FALSE) #1.0 for the resolution
clustree(Cluster3, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
Cluster3 <- RunUMAP(Cluster3, reduction = 'apca', dims = 1:26, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
Cluster3_p2 <- DimPlot(Cluster3, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("ADT Clustering") + theme_void() + NoLegend()
Cluster3_p2 <- Cluster3_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))
DefaultAssay(Cluster3) <- "RNA" #For log normalisation
Cluster3 <- FindMultiModalNeighbors(
  Cluster3, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:22), modality.weight.name = "RNA.weight")

Cluster3_p2
Idents(Cluster3)
#Clusters based on antibodies - analysis
Cluster3_RNA_c0 <- FindMarkers(Cluster3, ident.1 = 0, assay = "RNA")
Cluster3_RNA_c1 <- FindMarkers(Cluster3, ident.1 = 1, assay = "RNA")
Cluster3_RNA_c2 <- FindMarkers(Cluster3, ident.1 = 2, assay = "RNA")
Cluster3_RNA_c3 <- FindMarkers(Cluster3, ident.1 = 3, assay = "RNA")

Cluster3_ADT_c0 <- FindMarkers(Cluster3, ident.1 = 0, assay = "ADT")
Cluster3_ADT_c1 <- FindMarkers(Cluster3, ident.1 = 1, assay = "ADT")
Cluster3_ADT_c2 <- FindMarkers(Cluster3, ident.1 = 2, assay = "ADT")
Cluster3_ADT_c3 <- FindMarkers(Cluster3, ident.1 = 3, assay = "ADT")

DefaultAssay(Cluster3)<-"RNA"
Idents(Cluster3)
Cluster3_TFIDF_c0_genes <- WhichCells(object = Cluster3, ident = "0")
Cluster3_TFIDF_c0 <- tfidf(GetAssayData(Cluster3), Cluster3_TFIDF_c0_genes, colnames(Cluster3))

Cluster3_TFIDF_c1_genes <- WhichCells(object = Cluster3, ident = "1")
Cluster3_TFIDF_c1 <- tfidf(GetAssayData(Cluster3), Cluster3_TFIDF_c1_genes, colnames(Cluster3))

Cluster3_TFIDF_c2_genes <- WhichCells(object = Cluster3, ident = "2")
Cluster3_TFIDF_c2 <- tfidf(GetAssayData(Cluster3), Cluster3_TFIDF_c2_genes, colnames(Cluster3))

Cluster3_TFIDF_c3_genes <- WhichCells(object = Cluster3, ident = "3")
Cluster3_TFIDF_c3 <- tfidf(GetAssayData(Cluster3), Cluster3_TFIDF_c3_genes, colnames(Cluster3))

#Subset bt Genotype
Cluster3_E1020K <- subset(Cluster3, subset = Genotype == "E1020K")
Cluster3_WT <- subset(Cluster3, subset = Genotype != "E1020K")

DimPlot(Cluster3_E1020K, reduction = "adt.umap")
DimPlot(Cluster3_WT, reduction = "adt.umap")

#General analysis of Custer 3
CD4_cells_3 + CD38 + PD_1 + PD_L1 + ICOS

DefaultAssay(CD4_cells) <- "RNA"
LAG3 <- FeaturePlot(CD4_cells, features = "Lag3", cols = col_con, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
IL21 <- FeaturePlot(CD4_cells, features = "Il21", cols = col_con, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
IFNG <- FeaturePlot(CD4_cells, features = "Ifng", cols = col_con, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
ITGB1 <- FeaturePlot(CD4_cells, features = "Itgb1", cols = col_con, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))

FeaturePlot(CD4_cells, features = "Ccr9", cols = col_con, reduction = "wnn.umap", pt.size = 1.8) + theme_void() + theme(plot.title = element_text(color="black", size=15))
CD4_cells_p3

#Il21-hi Th1 effector cells
PD1_ICOS_CD73_CD200_Tcells <- list(Book2$PD1_ICOS_CD73_CD200_Tcells)
PD1_ICOS_CD73_CD200_Tcells  <- lapply(PD1_ICOS_CD73_CD200_Tcells, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = PD1_ICOS_CD73_CD200_Tcells , name = "PD1_ICOS_CD73_CD200_Tcells", search = FALSE)
VlnPlot(CD4_cells, c("PD1_ICOS_CD73_CD200_Tcells1"), pt.size  = 0.5)
FeaturePlot(CD4_cells, c("PD1_ICOS_CD73_CD200_Tcells1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

PD1_ICOS_CD73_CD200_Tcells_Th1_like <- list(Book2$PD1_ICOS_CD73_CD200_Tcells_Th1_like)
PD1_ICOS_CD73_CD200_Tcells_Th1_like  <- lapply(PD1_ICOS_CD73_CD200_Tcells_Th1_like, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = PD1_ICOS_CD73_CD200_Tcells_Th1_like, name = "PD-1_ICOS_CD73_CD200_Tcells_Th1_like", search = FALSE)
plot43 <- VlnPlot(CD4_cells, c("PD1_ICOS_CD73_CD200_Tcells_Th1_like1"), pt.size  = 0.5) + ggtitle("") + NoLegend() + xlab("")
plot44 <- FeaturePlot(CD4_cells, c("PD1_ICOS_CD73_CD200_Tcells_Th1_like1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for Th1-like PD1hi ICOShi CD73hi CD200hi CD4+ T cells")
grid.arrange(plot44, plot43, ncol = 2)

#Cluster 3 (and other E1020K-enriched clusters have an aging signature):
Age_associated.genes <- list(Book2$Age_associated.genes)
Age_associated.genes  <- lapply(Age_associated.genes, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Age_associated.genes , name = "Age_associated.genes", search = FALSE)
VlnPlot(CD4_cells, c("Age_associated.genes1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Age_associated.genes1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Age-associated CD4 subsets
Old_aTreg <- list(Book2$Old_aTreg)
Old_aTreg  <- lapply(Old_aTreg, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Old_aTreg, name = "Old_aTreg", search = FALSE)
p31 <- VlnPlot(CD4_cells, c("Old_aTreg1"), pt.size  = 0.5) + theme_void() + ggtitle("Old aTregs")
p32 <- FeaturePlot(CD4_cells, c("Old_aTreg1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Old aTregs")

Young_aTreg <- list(Book2$Young_aTreg)
Young_aTreg <- lapply(Young_aTreg, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Young_aTreg, name = "Young_aTreg", search = FALSE)
VlnPlot(CD4_cells, c("Young_aTreg1"), pt.size  = 0.5, cols = turbo(10))
p33 <- FeaturePlot(CD4_cells, c("Young_aTreg1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Young aTregs")

grid.arrange(CD4_cells_3, p32, p33, p31, ncol = 2)

Old_Exhausted <- list(Book2$Old_Exhausted)
Old_Exhausted <- lapply(Old_Exhausted, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Old_Exhausted, name = "Old_Exhausted", search = FALSE)
p28 <- VlnPlot(CD4_cells, c("Old_Exhausted1"), pt.size  = 0.5) + theme_void() + ggtitle("Old Exhausted")
p29 <- FeaturePlot(CD4_cells, c("Old_Exhausted1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Old Exhausted")

Young_Exhausted <- list(Book2$Young_Exhausted)
Young_Exhausted <- lapply(Young_Exhausted, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Young_Exhausted, name = "Young_Exhausted", search = FALSE)
VlnPlot(CD4_cells, c("Young_Exhausted1"), pt.size  = 0.5, cols = turbo(10))
p30 <- FeaturePlot(CD4_cells, c("Young_Exhausted1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)  + theme_void() + ggtitle("Young Exhausted")

grid.arrange(CD4_cells_3, p29, p30, p28, ncol = 2)

Old_Cytotoxic <- list(Book2$Old_Cytotoxic)
Old_Cytotoxic <- lapply(Old_Cytotoxic, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Old_Cytotoxic, name = "Old_Cytotoxic", search = FALSE)
p20 <- VlnPlot(CD4_cells, c("Old_Cytotoxic1"), pt.size  = 0.5) + ggtitle("Old Cytotoxic") + xlab("") + theme_void()
p21 <- FeaturePlot(CD4_cells, c("Old_Cytotoxic1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + ggtitle("Old Cytotoxic") + theme_void()

col2 <- (hue_pal()(12))

Young_Cytotoxic <- list(Book2$Young_Cytotoxic)
Young_Cytotoxic <- lapply(Young_Cytotoxic, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Young_Cytotoxic, name = "Young_Cytotoxic", search = FALSE)
p22 <- VlnPlot(CD4_cells, c("Young_Cytotoxic1"), pt.size  = 0.5, cols = turbo(10)) + ggtitle("Young Cytotoxic") + NoLegend() + xlab("")
p23 <- FeaturePlot(CD4_cells, c("Young_Cytotoxic1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + ggtitle("Young Cytotoxic") +theme_void()

p24 <- grid.arrange(p20, p22, ncol = 1)
p25 <- grid.arrange(CD4_cells_3, p21, p23, p20, ncol = 2)


Old_EM <- list(Book2$Old_EM)
Old_EM <- lapply(Old_EM, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Old_EM, name = "Old_EM", search = FALSE)
p34 <- VlnPlot(CD4_cells, c("Old_EM1"), pt.size  = 0.5) + theme_void() + ggtitle("Old Effector Memory")
p35 <- FeaturePlot(CD4_cells, c("Old_EM1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Old Effector Memory")

Young_EM <- list(Book2$Young_EM)
Young_EM <- lapply(Young_EM, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Young_EM, name = "Young_EM", search = FALSE)
VlnPlot(CD4_cells, c("Young_EM1"), pt.size  = 0.5, cols = turbo(10))
p36 <- FeaturePlot(CD4_cells, c("Young_EM1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Young Effector Memory")

grid.arrange(CD4_cells_3, p35, p36, p34, ncol = 2)

Old_rTreg <- list(Book2$Old_rTreg)
Old_rTreg <- lapply(Old_rTreg, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Old_rTreg, name = "Old_rTreg", search = FALSE)
VlnPlot(CD4_cells, c("Old_rTreg1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Old_rTreg1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

Young_rTreg <- list(Book2$Young_rTreg)
Young_rTreg <- lapply(Young_rTreg, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Young_rTreg, name = "Young_rTreg", search = FALSE)
VlnPlot(CD4_cells, c("Young_rTreg1"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Young_rTreg1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

Old_Naive_Isg15 <- list(Book2$Old_Naive_Isg15)
Old_Naive_Isg15 <- lapply(Old_Naive_Isg15, str_to_ti, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Old_Naive_Isg15, name = "Old_Naive_Isg15", search = FALSE)
VlnPlot(CD4_cells, c("Old_Naive_Isg151"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Old_Naive_Isg151"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

Young_Naive_Isg15 <- list(Book2$Young_Naive_Isg15)
Young_Naive_Isg15 <- lapply(Young_Naive_Isg15, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Young_Naive_Isg15, name = "Young_Naive_Isg15", search = FALSE)
VlnPlot(CD4_cells, c("Young_Naive_Isg151"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("Young_Naive_Isg151"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

Old_Naive <- list(Book2$Old_Naive)
Old_Naive <- lapply(Old_Naive, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Old_Naive, name = "Old_Naive", search = FALSE)
p26 <- VlnPlot(CD4_cells, c("Old_Naive1"), pt.size  = 0.5) + theme_void() + ggtitle("Old Naive")
p27 <- FeaturePlot(CD4_cells, c("Old_Naive1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + ggtitle("Old Naive") + theme_void()

Young_Naive <- list(Book2$Young_Naive)
Young_Naive <- lapply(Young_Naive, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = Young_Naive, name = "Young_Naive", search = FALSE)
VlnPlot(CD4_cells, c("Young_Naive1"), pt.size  = 0.5, cols = turbo(10))
p28 <- FeaturePlot(CD4_cells, c("Young_Naive1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + ggtitle("Young Naive") + theme_void()

grid.arrange(CD4_cells_3, p27, p28, p26, ncol = 2)

#Heatmap young and old popualtions
Data<-CD4_cells@meta.data
Forest <-Data %>%
  group_by(`Genotype`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$Genotype
Forest<-Forest[ , c(46:58)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)


Data<-CD4_cells@meta.data
Forest <-Data %>% 
  group_by(`seurat_clusters`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$seurat_clusters
Forest<-Forest[ , c(46:58, 64)]
colnames(Forest) <- c("Old aTregs", "Young aTregs", "Old Exhausted", "Young Exhausted", "Old Cytotoxic", "Young Cytotoxic", "Old Effector Memory",
                      "Old rTregs", "Young rTregs", "Old Naive Isg15", "Young Naive Isg15", "Old Naive", "Young Naive", "Young Effector Memory")
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = T, cluster_cols = T,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45)

#Investigating similarity between Cluster 3 cells and 4PD1hi cells
X4PD1_hi <- list(Book2$X4PD1_hi)
X4PD1_hi <- lapply(X4PD1_hi, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = X4PD1_hi, name = "X4PD1_hi", search = FALSE)
plot44 <- VlnPlot(CD4_cells, c("X4PD1_hi1"), pt.size  = 0.5) + NoLegend() + ggtitle("") + xlab("")
plot45 <- FeaturePlot(CD4_cells, c("X4PD1_hi1"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("4PD1hi")
grid.arrange(plot45, plot44, ncol = 2)


#Pan-cancer T cell signatures
IFNG_TFH_or_TH1 <- list(Book2$IFNG_TFH_or_TH1)
IFNG_TFH_or_TH1 <- lapply(IFNG_TFH_or_TH1, str_to_title)
DefaultAssay(CD4_cells) <- "RNA"
CD4_cells <-AddModuleScore(CD4_cells, features = IFNG_TFH_or_TH1, name = "IFNG_TFH_or_TH1", search = FALSE)
VlnPlot(CD4_cells, c("IFNG_TFH_or_TH11"), pt.size  = 0.5, cols = turbo(10))
FeaturePlot(CD4_cells, c("IFNG_TFH_or_TH11"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)


Idents(CD4_cells)
Idents(CD4_cells) <- CD4_cells$wsnn_res.0.9
levels(CD4_cells) <- Idents(CD4_cells)
FeaturePlot(CD4_cells, c("Il21"), cols=col_con, reduction = "wnn.umap", pt.size  = 1.5)

#Rename cluster idents
Idents(CD4_cells) <- CD4_cells$wsnn_res.0.9
CD4_cells<- RenameIdents(CD4_cells, `0` = "Naive CD4+ T cells", `1` ="Efffector Memmory and Cytotoxic CD4+ T cells", `2` ="T follicular helper cells", `3` = "Exhausted, IL21+ Th1 cells", `4` = "CD25lo T regulatory cells",`5` ="CD25hi T regulatory cells", `6` = "Activated/age-associated T regulatory cells", `7` = "Th17/Th9 cells", `8` = "Cytotoxic/Innate CD4+ T cells", `9` = "Contamination")
CD4_cells$seurat_clusters <- Idents(CD4_cells)
Idents(CD4_cells) <- CD4_cells$seurat_clusters
CD4_cells_p3 <- DimPlot(CD4_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Seurat Clusters") + theme_classic() + NoLegend()
CD4_cells_p3 <- CD4_cells_p3 + theme(plot.title = element_text(color="black", size=15, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")
CD4_cells_p3


Book2 <- read.csv(file.choose("CD4_Tcells_Signatures.csv"),header = T, sep = ',')
CD4_cells_p3

