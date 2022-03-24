## ##Subsetting: B cells####
#Julius Christopher Baeck
#CITE-Seq (1) analysis
#Subsetting of all B cells

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
experiment.lc <- LoadH5Seurat("experiment.lc.h5seurat")
experiment.lc.plain <- LoadH5Seurat("experiment.lc.plain.h5seurat")
head(experiment.lc[[]])
head(experiment.lc.plain[[]])

####Subsetting####
test <- experiment.lc.plain

##Based on clonotype abd x3 CD19 ADT
DefaultAssay(test) <-  "RNA"
Bcells <- subset(test, subset = cloneType != "NA") #12794
DefaultAssay(Bcells) <-  "ADT"
Bcells <- subset(Bcells, subset = Cd19 > 2) #12636
B_cells_p6

##Based on minimum of x3 CD19 ADT and 1x CD19 RNA
DefaultAssay(test) <-  "ADT"
Bcells1 <- subset(test, subset = Cd19 > 2)
DefaultAssay(Bcells1) <-  "RNA"
Bcells1 <- subset(Bcells1, subset = Cd19 > 0) #10592


Bcells4 <- subset(test, subset = Cd19 > 2)
Bcells4 <- subset(Bcells4, subset = Cd4 == 0)

##Based on other B cells associated genes and a minimum of x3 CD19 ADT 
DefaultAssay(test) <-  "RNA"
Bcells2 <- subset(test, subset = Ms4a1 > 0)
DefaultAssay(Bcells2) <-  "ADT"
Bcells2 <- subset(Bcells2, subset = Cd19 > 2) #13665

DefaultAssay(test) <-  "RNA"
Bcells3 <- subset(test, subset = Cd79a > 0)
DefaultAssay(Bcells3) <-  "ADT"
Bcells3 <- subset(Bcells3, subset = Cd19 > 2) #14360

DefaultAssay(test) <-  "RNA"
Bcells4 <- subset(test, subset = Cd79b > 0)
DefaultAssay(Bcells4) <-  "ADT"
Bcells4 <- subset(Bcells4, subset = Cd19 > 2) #14297

##Selecteion of Bcells (clonotype + CD19 ADT > 2) for further downstream analysis
B_cells <- Bcells
head(B_cells[[]])

###RNA####
###Normalise subset###
DefaultAssay(B_cells) <- "RNA" #For log normalisation
DefaultAssay(B_cells) <- "SCT" #For SCTransform

##RNA normalisation
B_cells <- NormalizeData(B_cells, verbose = TRUE)
B_cells <- FindVariableFeatures(B_cells, nfeatures = 3000)
B_cells <- ScaleData(B_cells)

#Or
B_cells <-  SCTransform(B_cells, verbose = TRUE)
B_cells[["SCT"]]

##Visualisation
top20 <-  head(VariableFeatures(B_cells), 20)
plot1.1 <-  VariableFeaturePlot(B_cells)
top20_plot <-  LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)

##RNA PCA
B_cells <- RunPCA(B_cells, verbose = FALSE, features = VariableFeatures(object = B_cells))
pca_variance <- B_cells@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #35

##RNA clustering
DefaultAssay(B_cells) <- "RNA" #For log normalisation
DefaultAssay(B_cells) <- "SCT" #For SCTransform

B_cells <- FindNeighbors(B_cells, dims = 1:35)
B_cells <- FindClusters(B_cells, resolution = 1.0, verbose = FALSE) #1.0 for the resolution
clustree(B_cells, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
B_cells <-RunUMAP(B_cells, dims = 1:35, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
B_cells_p1 <- DimPlot(B_cells, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("RNA Clustering") + theme_void() + NoLegend()
B_cells_p1 <- B_cells_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####ADT####
DefaultAssay(B_cells) <- "ADT"

##ADT normalisation
VariableFeatures(B_cells) <- rownames(B_cells[["ADT"]])
B_cells <- NormalizeData(B_cells, normalization.method = "CLR", margin = 2)
B_cells <- ScaleData(B_cells)
B_cells <- RunPCA(B_cells, reduction.name = 'apca', approx = FALSE)

##ADT PCA
apca_variance <- B_cells@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #27

##ADT clustering
B_cells <- FindNeighbors(B_cells, dims = 1:27, reduction = "apca")
B_cells <- FindClusters(B_cells, resolution = 1.0, verbose = FALSE) #1.0 for the resolution
clustree(B_cells, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
B_cells <- RunUMAP(B_cells, reduction = 'apca', dims = 1:27, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
B_cells_p2 <- DimPlot(B_cells, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("ADT Clustering") + theme_void() + NoLegend()
B_cells_p2 <- B_cells_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))

####WNN####
DefaultAssay(B_cells) <- "RNA" #For log normalisation
DefaultAssay(B_cells) <- "SCT" #For SCTransform

##Combine into wnn plot
B_cells <- FindMultiModalNeighbors(
  B_cells, reduction.list = list("pca", "apca"), 
  dims.list = list(1:35, 1:27), modality.weight.name = "RNA.weight")

##WNN clustering
B_cells <- FindClusters(B_cells, graph.name = "wsnn", algorithm = 3, resolution = 1.0, verbose = TRUE) #1.0 for the resolution
clustree(B_cells, prefix = "wsnn_res.") + theme(legend.position="bottom")
B_cells <- RunUMAP(B_cells, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
B_cells_p3 <- DimPlot(B_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("WNN Clustering") + theme_void() + NoLegend()
B_cells_p3 <- B_cells_p3 + theme(plot.title = element_text(color="black", size=25, face="bold"))

##Change Idents
Idents(B_cells)
B_cells[["old.ident"]] <- Idents(B_cells)
Idents(B_cells) <- B_cells[["orig.ident"]]
B_cells<- RenameIdents(B_cells, `c1` = "Mb1 Cyp11a1 KO 1", `c2` = "Mb1 Cyp11a1 KO 2", `d1` = "Mb1 E1020K Cyp11a1 KO 1", `d2` = "Mb1 E1020K Cyp11a1 KO 2", `a` = "WT 1", `b` = "WT 2", `f` = "E1020K")
B_cells[["orig.ident"]] <- Idents(B_cells)
Idents(B_cells) <- B_cells[["old.ident"]]
Idents(B_cells)

####Harmony reductions####
Idents(B_cells) <- B_cells[["orig.ident"]]
B_cells[["Genotype"]] <- Idents(B_cells)
B_cells <- RenameIdents(B_cells, `Mb1 Cyp11a1 KO 1` = "Batch 1", `Mb1 Cyp11a1 KO 2` = "Batch 1", `Mb1 E1020K Cyp11a1 KO 1` = "Batch 1", `Mb1 E1020K Cyp11a1 KO 2` = "Batch 1", `WT 1` = "Batch 2", `WT 2` = "Batch 2", `E1020K` = "Batch 2")
B_cells[["Batch"]] <- Idents(B_cells)
Idents(B_cells) <- B_cells[["old.ident"]]
head(B_cells[[]])

##Batch effect for RNA
DefaultAssay(B_cells) <- "RNA"
B_cells = RunHarmony(B_cells, "Batch", plot_convergence = TRUE, reduction = "pca", assay.use = "RNA", reduction.save = "harmony.rna")
harmony_rna_embeddings = Embeddings(B_cells, 'harmony.rna')
B_cells = B_cells %>%
  RunUMAP(reduction = "harmony.rna", dims = 1:35, assay = 'RNA', reduction.name = 'harmony.rna.umap', reduction.key = 'harmony.rnaUMAP_') %>%
  FindNeighbors(reduction = "harmony.rna", dims = 1:35) %>%
  FindClusters(resolution = 1.0, verbose = FALSE)
B_cells_p4 <- DimPlot(B_cells, label = TRUE, reduction = "harmony.rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Harmony RNA Clustering") + theme_void() + NoLegend()
B_cells_p4 <- B_cells_p4 + theme(plot.title = element_text(color="black", size=25, face="bold"))


##Batch effect for ADT
DefaultAssay(B_cells) <- "ADT"
B_cells = RunHarmony(B_cells, "Batch", plot_convergence = TRUE, reduction = "apca", assay.use = "ADT", reduction.save = "harmony.adt")
harmony_adt_embeddings = Embeddings(B_cells, 'harmony.adt')
B_cells = B_cells %>%
  RunUMAP(reduction = "harmony.adt", dims = 1:27, assay = 'ADT', reduction.name = 'harmony.adt.umap', reduction.key = 'harmony.adtUMAP_') %>%
  FindNeighbors(reduction = "harmony.adt", dims = 1:27) %>%
  FindClusters(resolution = 1.0, verbose = FALSE)
B_cells_p5 <- DimPlot(B_cells, label = TRUE, reduction = "harmony.adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE) +  ggtitle("Harmony ADT Clustering") + theme_void() + NoLegend()
B_cells_p5 <- B_cells_p5 + theme(plot.title = element_text(color="black", size=25, face="bold"))


##WNN clustering with harmony reductions
DefaultAssay(B_cells) <- "RNA"
B_cells <- FindMultiModalNeighbors(
  B_cells, reduction.list = list("harmony.rna", "harmony.adt"), 
  dims.list = list(1:35, 1:27), modality.weight.name = "harmony.weight", weighted.nn.name = "harmony.weighted.nn", snn.graph.name = "harmony.wsnn")
B_cells <- FindClusters(B_cells, graph.name = "harmony.wsnn", algorithm = 3, resolution = 0.6, verbose = TRUE) #0.6 for the resolution
B_cells <- RunUMAP(B_cells, nn.name = "harmony.weighted.nn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_")
B_cells_p6 <- DimPlot(B_cells, label = TRUE, reduction = "harmony.wnn.umap", pt.size = 1.2, label.size = 3, label.box = TRUE) +  ggtitle("Highlighted by cluster") + theme_void() + NoLegend()
B_cells_p6 <- B_cells_p6 + theme(plot.title = element_text(color="black", size=20, face="bold"))

Batch_rna.umap <- DimPlot(B_cells, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("RNA UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_adt.umap <- DimPlot(B_cells, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("ADT UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_wnn.umap <- DimPlot(B_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("WNN UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))

Batch_harmony.rna.umap <- DimPlot(B_cells, label = TRUE, reduction = "harmony.rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("Harmony RNA UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_harmony.adt.umap <- DimPlot(B_cells, label = TRUE, reduction = "harmony.adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("Harmony ADT UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))
Batch_harmony.wnn.umap <- DimPlot(B_cells, label = TRUE, reduction = "harmony.wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, group.by = "Batch") + theme_bw() + ggtitle("Harmony WNN UMAP") + NoLegend() + theme(plot.title = element_text(color="black", size=20, face="bold"))

Batch_plot <- Batch_rna.umap + Batch_adt.umap + Batch_wnn.umap + Batch_harmony.rna.umap + Batch_harmony.adt.umap + Batch_harmony.wnn.umap

B_cells
head(B_cells[[]])
SaveH5Seurat(B_cells, filename = "B_cells_2", overwrite = TRUE)

#Or
B_cells <- LoadH5Seurat("B_cells_2.h5seurat")


####Subset analysis####
DefaultAssay(B_cells) <- "ADT"
DefaultAssay(B_cells) <- "RNA"

##Set identity to harmony.wnn.umap clusters
Idents(B_cells) <- B_cells$harmony.wsnn_res.0.6
B_cells<- RenameIdents(B_cells, `0` = "Cluster 0", `1` ="Cluster 1", `2` ="Cluster 2", `3` = "Cluster 3", `4` = "Cluster 4",`5` ="Cluster 5", `6` = "Cluster 6", `7` = "Cluster 7", `8` = "Cluster 8", `9` = "Cluster 9", `10` = "Cluster 10", `11` = "Cluster 11", `12` = "Cluster 12", `13` = "Cluster 13", `14` = "Cluster 14", `15` = "Cluster 15", `16` = "Cluster 16")
B_cells$seurat_clusters <- Idents(B_cells)
Idents(B_cells) <- B_cells$harmony.wsnn_res.0.6


##By genotype
B_cells_Mouse <- DimPlot(B_cells, label = FALSE ,reduction = "harmony.wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, label.box = FALSE, repel = FALSE)
B_cells_Mouse$data$Genotype = factor(x = B_cells_Mouse$data$Genotype, levels = c("WT 1", 'WT 2', "E1020K", "Mb1 Cyp11a1 KO 1", "Mb1 Cyp11a1 KO 2", "Mb1 E1020K Cyp11a1 KO 1", "Mb1 E1020K Cyp11a1 KO 2"))
B_cells_Mouse <- B_cells_Mouse + theme_void() +
  ggtitle("Highlighted by genotype") +
  theme(plot.title = element_text(color="black", size=20, face="bold"))

show_col(hue_pal()(12))

Sample.a <- subset(B_cells, subset = Genotype == "WT 1")
Sample.a.plot <- DimPlot(Sample.a, label = FALSE ,reduction = "harmony.wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#F8766D")) +
  theme_void() +
  ggtitle("")

Sample.b <- subset(B_cells, subset = Genotype == "WT 2")
Sample.b.plot <- DimPlot(Sample.b, label = FALSE ,reduction = "harmony.wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#B79F00")) +
  theme_void() +
  ggtitle("")

Sample.c <- subset(B_cells, subset = Genotype == "E1020K")
Sample.c.plot <- DimPlot(Sample.c, label = FALSE ,reduction = "harmony.wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#7CAE00")) +
  theme_void() +
  ggtitle("")

Sample.d <- subset(B_cells, subset = Genotype == "Mb1 Cyp11a1 KO 1")
Sample.d.plot <- DimPlot(Sample.d, label = FALSE ,reduction = "harmony.wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#00C08B")) +
  theme_void() +
  ggtitle("")

Sample.e <- subset(B_cells, subset = Genotype == "Mb1 Cyp11a1 KO 2")
Sample.e.plot <- DimPlot(Sample.e, label = FALSE ,reduction = "harmony.wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#00B4F0")) +
  theme_void() +
  ggtitle("")

Sample.f <- subset(B_cells, subset = Genotype == "Mb1 E1020K Cyp11a1 KO 1")
Sample.f.plot <- DimPlot(Sample.f, label = FALSE ,reduction = "harmony.wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#C77CFF")) +
  theme_void() +
  ggtitle("")

Sample.g <- subset(B_cells, subset = Genotype == "Mb1 E1020K Cyp11a1 KO 2")
Sample.g.plot <- DimPlot(Sample.g, label = FALSE ,reduction = "harmony.wnn.umap", group.by = "Genotype", pt.size = 1.2, label.size = 6, cols = c("#F564E3")) +
  theme_void() +
  ggtitle("")

B_cells_p6 + B_cells_Mouse + Sample.a.plot + Sample.b.plot + Sample.c.plot + Sample.d.plot + Sample.e.plot + Sample.f.plot + Sample.g.plot

##By cell cycle genes
S.genes <- cc.genes.updated.2019$s.genes
S.genes <- lapply(S.genes, str_to_title)
G2M.genes <-  cc.genes.updated.2019$g2m.genes
G2M.genes <- lapply(G2M.genes, str_to_title)
B_cells <- CellCycleScoring(B_cells, s.features=S.genes, g2m.features=G2M.genes, set.ident = FALSE)
Idents(B_cells)

B_cells_cell_cylce <- DimPlot(B_cells, label = FALSE, reduction = "harmony.wnn.umap", group.by = "Phase", pt.size = 1.3) +  ggtitle("Highlighted by cell-cycle stage") + theme_void()
B_cells_cell_cylce <- B_cells_cell_cylce + theme(plot.title = element_text(color="black", size=25, face="bold"))

##Cell numbers and percent by genotype by clusters - harmony.wnn)
cell.numbers <- table(B_cells@meta.data$seurat_clusters, B_cells@meta.data$Genotype)
cell.numbers <- as.data.frame.matrix(cell.numbers)

cell_number_heatmap <- pheatmap::pheatmap(t(cell.numbers), cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                                          cellwidth = 30,cellheight = 30, angle_col = 45, col = col_con)

B_cells_meta <- B_cells@meta.data
genotype.numbers <- B_cells_meta %>% dplyr::count(Genotype)
genotype.numbers.vector <- genotype.numbers %>% pull(n)
str(genotype.numbers.vector)
cell.percent <- sweep(cell.numbers, 2, genotype.numbers.vector, "/")
cell.percent <- cell.percent*100

cell_percent_heatmap <- pheatmap::pheatmap(t(cell.percent), cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                                           cellwidth = 30,cellheight = 30, angle_col = 45, col = col_con)

cell.percent <- tibble::rownames_to_column(cell.numbers, "Cluster")
write_xlsx(cell.numbers, "cell.numbers.xlsx")

##Looking for specific ADT markers
DefaultAssay(B_cells) <- "ADT"

B220 <- FeaturePlot(B_cells, features = "B220", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.5) + theme_void() + ggtitle("B220") + theme(plot.title = element_text(color="black", size=15, face="bold"))
Cd19 <- FeaturePlot(B_cells, features = "Cd19", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD19")
IgM <- FeaturePlot(B_cells, features = "Igm", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("IgM")
IgD <- FeaturePlot(B_cells, features = "Igd", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("IgD")
CD23 <- FeaturePlot(B_cells, features = "Cd23", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD23")
CD21 <- FeaturePlot(B_cells, features = "Cd21", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD21")
CD93 <- FeaturePlot(B_cells, features = "Cd93", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD93")
CD40 <- FeaturePlot(B_cells, features = "Cd40", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD40")


PD_L1 <- FeaturePlot(B_cells, features = "Pd-L1", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("PD-L1")
PD_L2 <- FeaturePlot(B_cells, features = "Pd-L2", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("PD-L2")
PD_1 <- FeaturePlot(B_cells, features = "Pd-1", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("PD-1")
CTLA4 <- FeaturePlot(B_cells, features = "Ctla4", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CTLA-4") + theme(plot.title = element_text(color="black", size=15, face="bold"))
CD95 <- FeaturePlot(B_cells, features = "Cd95", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD95")
CXCR5 <- FeaturePlot(B_cells, features = "Cxcr5", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CXCR5")
CD38 <- FeaturePlot(B_cells, features = "Cd38", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD38")
CD80 <- FeaturePlot(B_cells, features = "Cd80", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD80")
CD83 <- FeaturePlot(B_cells, features = "Cd83", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD83")
CD86 <- FeaturePlot(B_cells, features = "Cd86", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD86")

CD4 <- FeaturePlot(B_cells, features = "Cd4", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD4")
CD8 <- FeaturePlot(B_cells, features = "Cd8a", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD8")
CD44 <- FeaturePlot(B_cells, features = "Cd44", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.5) + theme_void() + ggtitle("CD44")

##Cluster identification
B_cells_p6 + B_cells_Mouse + B_cells_cell_cylce
B_cells_p6 + B220
B_cells_p6 + Cd19 + CD4 + CD8
B_cells_p6 + IgM + IgD + CD23 + CD21 + CD93
B_cells_p6 + CD95 + CD38 + CD86 + CXCR4 + AID
B_cells_p6 + PD_L2 + CD80 + Nt5e
B_cells_p6 + CTLA4 + Ctla4
B_cells_p6 + CD80 + CD86 + PD_L1 + CXCR5 + CD44



Nt5e <- VlnPlot(B_cells, features = "Nt5e") + ggtitle("CD73") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Apoe <- VlnPlot(B_cells, features = "Apoe") + ggtitle("APOE") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Igkv5_39 <- VlnPlot(B_cells, features = "Igkv5-39") + ggtitle("Igkv5-39") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Fcer1g <- VlnPlot(B_cells, features = "Fcer1g") + ggtitle("Fcer1g") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Ighg3 <- VlnPlot(B_cells, features = "Ighg3") + ggtitle("Ighg3") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Cd300lf <- VlnPlot(B_cells, features = "Cd300lf") + ggtitle("") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Cd80 <- VlnPlot(B_cells, features = "Cd80") + ggtitle("") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Ctla4 <- VlnPlot(B_cells, features = "Ctla4") + ggtitle("CTLA-4") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Cxcr4 <- VlnPlot(B_cells, features = "Cxcr4") + ggtitle("Cxcr4") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Aid <- VlnPlot(B_cells, features = "Ctla4") + ggtitle("Aid") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Cyp11a1 <- VlnPlot(B_cells, features = "Cyp11a1") + ggtitle("Cyp11a1") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Ighm <- VlnPlot(B_cells, features = "Ighm") + ggtitle("Ighm") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
Ighd <- VlnPlot(B_cells, features = "Ighd") + ggtitle("Ighd") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))

##Visualise ADT markers with dotplot
DefaultAssay(B_cells) <- "ADT"
Abs_Bcells <- c("Cd19", "B220", "Igm", "Igd", "Cd21", "Cd23", "Cd93", "Cxcr5", "Cd38", "Cd86", "Cd80", "Cd83", "Pd-L1", "Pd-L2", "Pd-1", "Ctla4", "Cd95", "Cd40", "Cd44")

Abs_1 <- c("B220", "Igm","Igd", "Cd23", "Cd21", "Cd93", "Cd38", "Cd40")
Abs_c3 <- c("Cd19", "B220","Cd80", "Pd-L1", "Cd86", "Cd44", "Cxcr5")

DotPlot(B_cells, features = Abs_1, cluster.idents = TRUE) + RotatedAxis() + scale_colour_viridis()
DotPlot(B_cells, features = c("Fcer2", "Cd4"), cluster.idents = TRUE) + scale_colour_viridis()

DefaultAssay(B_cells) <- "RNA"
Abs_Bcells_genes <- c("Cd19", "Ighm", "Ighd", "Cr2", "Cd93", "Cxcr5", "Cd38", "Cd86", "Cd80", "Cd83", "Cd274", "Pdcd1lg2", "Pdcd1", "Ctla4", "Fas", "Cd40", "Cd44")
DotPlot(B_cells, features = Abs_Bcells_genes, cluster.idents = TRUE) + RotatedAxis() + scale_colour_viridis()

StackedVlnPlot(B_cells, features = c("Cd21","Cd23"), cols = col_con)

##Find markers in each cluster
head(B_cells[[]])
Idents(object = B_cells) <- B_cells$harmony.wsnn_res.0.6
B_cells$seurat_clusters <- B_cells$harmony.wsnn_res.0.6

Bcells_RNA_c0 <- FindMarkers(B_cells, ident.1 = 0, assay = "RNA")
Bcells_ADT_c0 <- FindMarkers(B_cells, ident.1 = 0, assay = "ADT")

Bcells_RNA_c1 <- FindMarkers(B_cells, ident.1 = 1, assay = "RNA")
Bcells_ADT_c1 <- FindMarkers(B_cells, ident.1 = 1, assay = "ADT")

Bcells_RNA_c2 <- FindMarkers(B_cells, ident.1 = 2, assay = "RNA")
Bcells_ADT_c2 <- FindMarkers(B_cells, ident.1 = 2, assay = "ADT")

Bcells_RNA_c3 <- FindMarkers(B_cells, ident.1 = 3, assay = "RNA")
Bcells_ADT_c3 <- FindMarkers(B_cells, ident.1 = 3, assay = "ADT")

Bcells_RNA_c4 <- FindMarkers(B_cells, ident.1 = 4, assay = "RNA")
Bcells_ADT_c4 <- FindMarkers(B_cells, ident.1 = 4, assay = "ADT")

Bcells_RNA_c5 <- FindMarkers(B_cells, ident.1 = 5, assay = "RNA")
Bcells_ADT_c5 <- FindMarkers(B_cells, ident.1 = 5, assay = "ADT")

Bcells_RNA_c6 <- FindMarkers(B_cells, ident.1 = 6, assay = "RNA")
Bcells_ADT_c6 <- FindMarkers(B_cells, ident.1 = 6, assay = "ADT")

Bcells_RNA_c7 <- FindMarkers(B_cells, ident.1 = 7, assay = "RNA")
Bcells_ADT_c7 <- FindMarkers(B_cells, ident.1 = 7, assay = "ADT")

Bcells_RNA_c8 <- FindMarkers(B_cells, ident.1 = 8, assay = "RNA")
Bcells_ADT_c8 <- FindMarkers(B_cells, ident.1 = 8, assay = "ADT")

Bcells_RNA_c9 <- FindMarkers(B_cells, ident.1 = 9, assay = "RNA")
Bcells_ADT_c9 <- FindMarkers(B_cells, ident.1 = 9, assay = "ADT")

Bcells_RNA_c10 <- FindMarkers(B_cells, ident.1 = 10, assay = "RNA")
Bcells_ADT_c10 <- FindMarkers(B_cells, ident.1 = 10, assay = "ADT")

Bcells_RNA_c11 <- FindMarkers(B_cells, ident.1 = 11, assay = "RNA")
Bcells_ADT_c11 <- FindMarkers(B_cells, ident.1 = 11, assay = "ADT")

Bcells_RNA_c12 <- FindMarkers(B_cells, ident.1 = 12, assay = "RNA")
Bcells_ADT_c12 <- FindMarkers(B_cells, ident.1 = 12, assay = "ADT")

Bcells_RNA_c13 <- FindMarkers(B_cells, ident.1 = 13, assay = "RNA")
Bcells_ADT_c13 <- FindMarkers(B_cells, ident.1 = 13, assay = "ADT")

Bcells_RNA_c14 <- FindMarkers(B_cells, ident.1 = 14, assay = "RNA")
Bcells_ADT_c14 <- FindMarkers(B_cells, ident.1 = 14, assay = "ADT")

Bcells_RNA_c15 <- FindMarkers(B_cells, ident.1 = 15, assay = "RNA")
Bcells_ADT_c15 <- FindMarkers(B_cells, ident.1 = 15, assay = "ADT")

Bcells_RNA_c16 <- FindMarkers(B_cells, ident.1 = 16, assay = "RNA")
Bcells_ADT_c16 <- FindMarkers(B_cells, ident.1 = 16, assay = "ADT")

Bcells_RNA_c0 <- tibble::rownames_to_column(Bcells_RNA_c0, "Genes")
write_xlsx(Bcells_RNA_c0, "Bcells_RNA_c0.xlsx")

##Most abundant genes
DefaultAssay(B_cells)<-"RNA"
Bcells_SCT_c13_TFIDF <- WhichCells(object = B_cells, ident = "13")
Bcells_SCT_c13_TFIDF_genes <- tfidf(GetAssayData(B_cells), Bcells_SCT_c13_TFIDF, colnames(B_cells))
Bcells_SCT_c3_TFIDF_genes <- tibble::rownames_to_column(Bcells_SCT_c3_TFIDF_genes, "Genes")
write_xlsx(Bcells_SCT_c3_TFIDF_genes, "Bcells_SCT_c3_TFIDF_genes.xlsx")

B_cells_p6 + B220
FeaturePlot(B_cells, features = "Sox4", reduction = "harmony.wnn.umap")
VlnPlot(B_cells, features = "Ybx3", cols = turbo(17)) + NoLegend()

##Cluster vs cluster
c6_vs_c3_abs <- FindMarkers(B_cells, ident.1 = 6, ident.2 = 3, assay = "ADT")
c3_vs_c6_abs <- FindMarkers(B_cells, ident.1 = 3, ident.2 = 6, assay = "RNA")

##Clonotype
DimPlot(B_cells, group.by = "cloneType", reduction = "harmony.wnn.umap")

##Looking for specific RNA markers
DefaultAssay(B_cells)<-"RNA"
AID <- FeaturePlot(B_cells, features = "Aicda", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
AHR <- FeaturePlot(B_cells, features = "Ahr", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
CD5 <- FeaturePlot(B_cells, features = "Cd5", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
NOTCH2 <- FeaturePlot(B_cells, features = "Notch2", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
IL10 <- FeaturePlot(B_cells, features = "Il10", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
TGFb <- FeaturePlot(B_cells, features = "Tgfb1", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
CD24 <- FeaturePlot(B_cells, features = "Cd24a", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
BCL6 <- FeaturePlot(B_cells, features = "Bcl6", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
TBX21 <- FeaturePlot(B_cells, features = "Tbx21", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("Tbx21") + theme(plot.title = element_text(color="black", size=15))
ITGAM <- FeaturePlot(B_cells, features = "Itgam", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
PRDM1 <- FeaturePlot(B_cells, features = "Prdm1", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() 
SPN <- FeaturePlot(B_cells, features = "Spn", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + ggtitle("CD43") + theme(plot.title = element_text(color="black", size=15))
TACI <- FeaturePlot(B_cells, features = "Tnfrsf13b", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
CCR6 <- FeaturePlot(B_cells, features = "Ccr6", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
CD9 <- FeaturePlot(B_cells, features = "Cd9", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
ZBTB32 <- FeaturePlot(B_cells, features = "Zbtb32", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
BACH2 <- FeaturePlot(B_cells, features = "Bach2", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
CXCR4 <- FeaturePlot(B_cells, features = "Cxcr4", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
NT5E <- FeaturePlot(B_cells, features = "Nt5e", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void()
APOE <- FeaturePlot(B_cells, features = "Apoe", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
Gene_1810046k07rik <-  FeaturePlot(B_cells, features = "1810046k07rik", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
SDC1 <- FeaturePlot(B_cells, features = "Sdc1", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
CYP11A1 <- FeaturePlot(B_cells, features = "Cyp11a1", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
IGHM <- FeaturePlot(B_cells, features = "Ighm", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
IGHD <- FeaturePlot(B_cells, features = "Ighd", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
SLPI <- FeaturePlot(B_cells, features = "Slpi", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
SOX4 <- FeaturePlot(B_cells, features = "Sox4", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15))
FeaturePlot(B_cells, features = "Dntt", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2) + theme_void() + theme(plot.title = element_text(color="black", size=15)) + B_cells_p6

VlnPlot(B_cells, features = "Sox4") + B_cells_p6


B_cells_p6

##Visualisation of genes with DotPlot
DefaultAssay(B_cells)<-"RNA"
interesting_Bcell_genes <- c("Aicda", "Ahr", "Cd5", "Notch2", "Il10", "Tgfb1","Cd24a", "Bcl6", "Tbx21", "Itgam", "Prdm1", "Spn", "Tnfrsf13b")
DotPlot(B_cells, features = interesting_Bcell_genes, cluster.idents = TRUE) + RotatedAxis() + scale_colour_viridis()
StackedVlnPlot(B_cells, features = c("Apoe", "Ahr", "Ccdc28b", "Fcer1g"))


##Cyp11a1 expression
FeaturePlot(B_cells, features = "Cyp11a1", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 3) +
  theme_void() + ggtitle("Cyp11a1") + theme(plot.title = element_text(color="black", size=30))

VlnPlot(B_cells, features = "Cyp11a1", pt.size  = 0.5, group.by = "Genotype") + NoLegend()
VlnPlot(B_cells, features = "Cyp11a1", pt.size  = 0.5, group.by = "seurat_clusters") +
  NoLegend() + xlab("") + ggtitle("")



##Number and percentage of cells expressing Cyp11a1 (across all cells)
sum(GetAssayData(object = B_cells, slot = "data")["Cyp11a1",]>0) #50
sum(GetAssayData(object = B_cells, slot = "data")["Cyp11a1",]>0)/nrow(B_cells@meta.data) #1.41%

sumCyp <- subset(B_cells, subset = Cyp11a1 > 0)

cell.numbers.cyp <- table(sumCyp@meta.data$seurat_clusters, sumCyp@meta.data$Genotype)
cell.numbers.cyp <- as.data.frame.matrix(cell.numbers.cyp)

cell_number_heatmap_cyp <- pheatmap::pheatmap(t(cell.numbers.cyp), cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                                              cellwidth = 30,cellheight = 30, angle_col = 45, col = col_con)

sumCyp_meta <- sumCyp@meta.data
genotype.numbers.cyp <- sumCyp_meta %>% dplyr::count(Genotype, .drop = FALSE)
genotype.numbers.vector.cyp <- genotype.numbers.cyp %>% pull(n)
str(genotype.numbers.vector.cyp)
cell.percent.cyp <- sweep(cell.numbers.cyp, 2, genotype.numbers.vector.cyp, "/")
cell.percent.cyp <- cell.percent.cyp*100

cell_percent_heatmap.cyp<- pheatmap::pheatmap(t(cell.percent.cyp), cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                                              cellwidth = 30,cellheight = 30, angle_col = 45, col = col_con)

##Identifying markers for cluster 3
FeaturePlot(B_cells, features = "S100a6", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 2) + theme_void()
VlnPlot(B_cells, features = "S100a6")
Top_10_c3 <- VlnPlot(B_cells, features = c("Zbtb32", "Ahr", "Cyp11a1", "Apoe"), group.by = "seurat_clusters", pt.size = 0.5)

##Differential gene expression between cluster 3 and cluster 8
Cluster_3_vs_8 <- FindMarkers(B_cells, ident.1 = 3, ident.2 = 8)
Cluster_8_vs_3 <- FindMarkers(B_cells, ident.1 = 8, ident.2 = 3)

###Pathway analysis
#Select gene set collection(s) (e.g. MySigDB "H" or "C7")
murine_gene_sets_c6 <- getGeneSets(species = "Mus musculus", library = "C6")


#Subset ADT cluster 4 for further downstream analysis
B_cells_c4 <- subset(B_cells, idents = "4")
DimPlot(B_cells_c4, reduction = "adt.umap")


en_B_cells_c4 <- enrichIt(obj = B_cells_c4, gene.sets = murine_gene_sets_c6, groups = 1000, cores = 4)
B_cells_c4 <- Seurat::AddMetaData(B_cells_c4, en_B_cells_c4)
head(B_cells_c4@meta.data)
dittoHeatmap(B_cells_c4, genes = NULL, metas = names(en_B_cells_c4), 
             fontsize = 7, 
             cluster_cols = TRUE,
             heatmap.colors = col_con)

###GSEA
Book1 <- read.csv("Book1.csv",header = T, sep = ',')
Idents(B_cells) <-  B_cells$seurat_clusters

##B220- bulk RNASeq data
E1_B220_neg_vs_pos <- list(Book1$E1_B220_neg_vs_pos)
E1_B220_neg_vs_pos <- lapply(E1_B220_neg_vs_pos, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = E1_B220_neg_vs_pos, name = "E1_B220_neg_vs_pos", n = 22, search = TRUE)
plot1 <- VlnPlot(B_cells, c("E1_B220_neg_vs_pos1"), pt.size  = 0.5) + NoLegend() +
  ggtitle("Enrichment for E1020K B220- vs E1020K B220+ differentially expressed genes") + 
  theme(plot.title = element_text(color="black", size=30, face="bold")) +
  xlab("") + ylab("Enrichment score")
plot2.1 <- FeaturePlot(B_cells, c("E1_B220_neg_vs_pos1"), cols=col_con, reduction = "harmony.wnn.umap",pt.size  = 2.8) +
  ggtitle("Enrichment for E1020K B220- vs E1020K B220+ differentially expressed genes") +
  theme_void() +
  theme(plot.title = element_text(color="black", size=30, face="bold"),
        legend.text = element_text(size = 15))

WT_B220_neg_vs_pos <- list(Book1$WT_B220_neg_vs_pos)
WT_B220_neg_vs_pos <- lapply(WT_B220_neg_vs_pos, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = WT_B220_neg_vs_pos, name = "WT_B220_neg_vs_pos", n = 22, search = TRUE)
plot2 <- VlnPlot(B_cells, c("WT_B220_neg_vs_pos1"), pt.size  = 0.5) + NoLegend() +
  ggtitle("Enrichment for WT B220- vs WT B220+ differentially expressed genes") +
  theme(plot.title = element_text(color="black", size=30, face="bold")) +
  xlab("") + ylab("Enrichment score")
plot2.2 <- FeaturePlot(B_cells, c("WT_B220_neg_vs_pos1"), cols=col_con, reduction = "harmony.wnn.umap",pt.size  = 2.8) +
  ggtitle("Enrichment for WT B220- vs WT B220+ differentially expressed genes") +
  theme_void() +
  theme(plot.title = element_text(color="black", size=30, face="bold"),
        legend.text = element_text(size = 15))


B220_neg_E1_vs_WT <- list(Book1$B220_neg_E1_vs_WT)
B220_neg_E1_vs_WT <- lapply(B220_neg_E1_vs_WT, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = B220_neg_E1_vs_WT, name = "B220_neg_E1_vs_WT", n = 22, search = TRUE)
plot3 <- VlnPlot(B_cells, c("B220_neg_E1_vs_WT1"), pt.size  = 0.5) + NoLegend()  +
  ggtitle("Enrichment for E1020K B220- vs WT B220- differentially expressed genes") +
  theme(plot.title = element_text(color="black", size=30, face="bold")) +
  xlab("") + ylab("Enrichment score")
plot2.3 <- FeaturePlot(B_cells, c("B220_neg_E1_vs_WT1"), cols=col_con, reduction = "harmony.wnn.umap",pt.size  = 2.8) +
  ggtitle("Enrichment for E1020K B220- vs WT B220- differentially expressed genes") +
  theme_void() +
  theme(plot.title = element_text(color="black", size=30, face="bold"),
        legend.text = element_text(size = 15))

plot2.1 + B_cells_p6 + plot2
plot1
grid.arrange(plot2.3, B_cells_p6, ncol = 2)

#Heatmap for bulk RNASeq signatures
Data <- B_cells@meta.data
Forest <- Data %>% 
  group_by(`seurat_clusters`) %>% 
  summarize_all(mean)
Forest <- as.data.frame(Forest)
rownames(Forest) <- Forest$seurat_clusters
colnames(Forest)
Forest <- Forest[ , c(28:30)]
colnames(Forest) <-  c("E1020K B220- vs E1020K B220+", "WT B220- vs WT B220+", "E1020K B220- vs WT B220-")
final <- t(Forest)

Data <- B_cells@meta.data
Forest <- Data %>% 
  group_by(`Genotype`) %>% 
  summarize_all(mean)
Forest <- as.data.frame(Forest)
rownames(Forest) <- Forest$Genotype
Forest <- Forest[ , c(28:30)]
colnames(Forest) <-  c("E1020K B220- vs E1020K B220+", "WT B220- vs WT B220+", "E1020K B220- vs WT B220-")
final <- t(Forest)


heatmap_1 <- pheatmap::pheatmap(final, cluster_rows = T, cluster_cols = T, show_rownames = T, show_colnames = T,
                                cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)

ggsave("Heatmap_cluster_bulkRNAsig.png", width = 20, height = 10)

##Memory vs naive
GSE11386_NAIVE_VS_MEMORY_BCELL_UP <- list(Book1$GSE11386_NAIVE_VS_MEMORY_BCELL_UP)
GSE11386_NAIVE_VS_MEMORY_BCELL_UP <- lapply(GSE11386_NAIVE_VS_MEMORY_BCELL_UP, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE11386_NAIVE_VS_MEMORY_BCELL_UP, name = "GSE11386_NAIVE_VS_MEMORY_BCELL_UP", n = 22, search = TRUE)
VlnPlot(B_cells, c("GSE11386_NAIVE_VS_MEMORY_BCELL_UP1"), pt.size  = 0.5) + NoLegend()
FeaturePlot(B_cells, c("GSE11386_NAIVE_VS_MEMORY_BCELL_UP1"), cols=col_con, reduction = "harmony.wnn.umap",pt.size  = 1.5)

GSE11386_NAIVE_VS_MEMORY_BCELL_DN <- list(Book1$GSE11386_NAIVE_VS_MEMORY_BCELL_DN)
GSE11386_NAIVE_VS_MEMORY_BCELL_DN <- lapply(GSE11386_NAIVE_VS_MEMORY_BCELL_DN, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE11386_NAIVE_VS_MEMORY_BCELL_DN, name = "GSE11386_NAIVE_VS_MEMORY_BCELL_DN", n = 22, search = TRUE)
VlnPlot(B_cells, c("GSE11386_NAIVE_VS_MEMORY_BCELL_DN1"), pt.size  = 0.5) + NoLegend()
FeaturePlot(B_cells, c("GSE11386_NAIVE_VS_MEMORY_BCELL_DN1"), cols=col_con, reduction = "harmony.wnn.umap",pt.size  = 1.5)

GSE11961_MEMORY_BCELL_DAY7_VS_GERMINAL_CENTER_BCELL_DAY7_UP <- list(Book1$GSE11961_MEMORY_BCELL_DAY7_VS_GERMINAL_CENTER_BCELL_DAY7_UP)
GSE11961_MEMORY_BCELL_DAY7_VS_GERMINAL_CENTER_BCELL_DAY7_UP <- lapply(GSE11961_MEMORY_BCELL_DAY7_VS_GERMINAL_CENTER_BCELL_DAY7_UP, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE11961_MEMORY_BCELL_DAY7_VS_GERMINAL_CENTER_BCELL_DAY7_UP, name = "GSE11961_MEMORY_BCELL_DAY7_VS_GERMINAL_CENTER_BCELL_DAY7_UP", n = 22, search = FALSE)
VlnPlot(B_cells, c("GSE11961_MEMORY_BCELL_DAY7_VS_GERMINAL_CENTER_BCELL_DAY7_UP1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("GSE11961_MEMORY_BCELL_DAY7_VS_GERMINAL_CENTER_BCELL_DAY7_UP1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

PreMemory_precursor <- list(Book1$PreMemory_precursor)
PreMemory_precursor <- lapply(PreMemory_precursor, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = PreMemory_precursor, name = "PreMemory_precursor", n = 22, search = FALSE)
VlnPlot(B_cells, c("PreMemory_precursor1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("PreMemory_precursor1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

Memory_B <- list(Book1$Memory_B)
Memory_B <- lapply(Memory_B, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Memory_B, name = "Memory_B", n = 22, search = FALSE)
VlnPlot(B_cells, c("Memory_B1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("Memory_B1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
FeaturePlot(B_cells, c("Fah"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

##Plasmablast signatures
Plasmablasts <- list(Book1$Plasmablasts)
Plasmablasts<- lapply(Plasmablasts, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Plasmablasts, name = "Plasmablasts", n = 22, search = FALSE)
VlnPlot(B_cells, c("Plasmablasts1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("Plasmablasts1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

##B1a B cell signatures
B1a <- list(Book1$B1a)
B1a<- lapply(B1a, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = B1a, name = "B1a", n = 22, search = FALSE)
B1a_plot1 <- VlnPlot(B_cells, c("B1a1"), pt.size  = 0.5, cols = turbo(17)) + NoLegend() + ggtitle("Enrichment for B1a signature") + theme(plot.title = element_text(color="black", size=15, face="bold"))
B1a_plot2 <- FeaturePlot(B_cells, c("B1a1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for B1a signature") + theme(plot.title = element_text(color="black", size=15, face="bold"))

B1_Bcells <- list(Book1$B1_Bcells)
B1_Bcells <- lapply(B1_Bcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = B1_Bcells, name = "B1_Bcells", n = 22, search = FALSE)
B1_Bcells_plot1 <- VlnPlot(B_cells, c("B1_Bcells1"), pt.size  = 0.5, cols = turbo(17)) + NoLegend() + ggtitle("Enrichment for B1 signature") + theme(plot.title = element_text(color="black", size=15, face="bold"))
B1_Bcells_plot2 <- FeaturePlot(B_cells, c("B1_Bcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for B1 signature") + theme(plot.title = element_text(color="black", size=15, face="bold"))
grid.arrange(B1_Bcells_plot2, B1_Bcells_plot1, ncol = 2)

B_cells_p6

##Breg signatures
Breg <- list(Book1$Breg)
Breg <- lapply(Breg, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Breg, name = "Breg", n = 22, search = TRUE)
Breg_plot1 <- VlnPlot(B_cells, c("Breg1"), pt.size  = 0.5, cols = turbo(17)) + NoLegend() + ggtitle("Enrichment for Breg signature") + theme(plot.title = element_text(color="black", size=15, face="bold"))
Breg_plot2 <-FeaturePlot(B_cells, c("Breg1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)  + theme_void() + ggtitle("Enrichment for Breg signature") + theme(plot.title = element_text(color="black", size=15, face="bold"))

grid.arrange(Breg_plot2, Breg_plot1, ncol = 2)

Breg_overlapping <- list(Book1$Breg_overlapping)
Breg_overlapping <- lapply(Breg_overlapping, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Breg_overlapping, name = "Breg_overlapping", n = 22, search = FALSE)
VlnPlot(B_cells, c("Breg_overlapping1"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(B_cells, c("Breg_overlapping1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

Breg_neg_imm_reg <- list(Book1$Breg_neg_imm_reg)
Breg_neg_imm_reg <- lapply(Breg_neg_imm_reg, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Breg_neg_imm_reg, name = "Breg_neg_imm_reg", n = 22, search = FALSE)
VlnPlot(B_cells, c("Breg_neg_imm_reg1"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(B_cells, c("Breg_neg_imm_reg1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

B10 <- list(Book1$B10)
B10 <- lapply(B10, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = B10, name = "B10", n = 22, search = FALSE)
VlnPlot(B_cells, c("B101"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(B_cells, c("B101"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

Data <- B_cells@meta.data
Forest <- Data %>% 
  group_by(`Genotype`) %>% 
  summarize_all(mean)
Forest <- as.data.frame(Forest)
rownames(Forest) <- Forest$Genotype
colnames(Forest)
Forest <- Forest[ , c(31:34)]
colnames(Forest) <-  c("Common Breg genes across tissues","Overlapping between negative regulator and B10 genes", "Negative regulation of the immune system", "B10 Breg signature")
final <- t(Forest)


heatmap_1 <- pheatmap::pheatmap(final, cluster_rows = T, cluster_cols = T, show_rownames = T, show_colnames = T,
                                cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)

##Follicular B cells
GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY40_UP <- list(Book1$GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY40_UP)
GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY40_UP <- lapply(GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY40_UP, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY40_UP, name = "GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY40_UP", n = 22, search = TRUE)
VlnPlot(B_cells, c("GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY40_UP1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY40_UP1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY7_UP <- list(Book1$GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY7_UP)
GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY7_UP <- lapply(GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY7_UP, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY7_UP, name = "GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY7_UP", n = 22, search = TRUE)
VlnPlot(B_cells, c("GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY7_UP1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("GSE11961_FOLLICULAR_BCELL_VS_GERMINAL_CENTER_BCELL_DAY7_UP1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN <- list(Book1$GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN)
GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN <- lapply(GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN, name = "GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN", n = 22, search = FALSE)
VlnPlot(B_cells, c("GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

GSE39916_B_CELL_SPLEEN_VS_PLASMA_CELL_BONE_MARROW_UP <- list(Book1$GSE39916_B_CELL_SPLEEN_VS_PLASMA_CELL_BONE_MARROW_UP)
GSE39916_B_CELL_SPLEEN_VS_PLASMA_CELL_BONE_MARROW_UP <- lapply(GSE39916_B_CELL_SPLEEN_VS_PLASMA_CELL_BONE_MARROW_UP, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE39916_B_CELL_SPLEEN_VS_PLASMA_CELL_BONE_MARROW_UP, name = "GSE39916_B_CELL_SPLEEN_VS_PLASMA_CELL_BONE_MARROW_UP", n = 22, search = FALSE)
VlnPlot(B_cells, c("GSE39916_B_CELL_SPLEEN_VS_PLASMA_CELL_BONE_MARROW_UP1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("GSE39916_B_CELL_SPLEEN_VS_PLASMA_CELL_BONE_MARROW_UP1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

GSE3039_B2_VS_B1_BCELL_DN <- list(Book1$GSE3039_B2_VS_B1_BCELL_DN)
GSE3039_B2_VS_B1_BCELL_DN <- lapply(GSE3039_B2_VS_B1_BCELL_DN, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE3039_B2_VS_B1_BCELL_DN, name = "GSE3039_B2_VS_B1_BCELL_DN", n = 22, search = FALSE)
VlnPlot(B_cells, c("GSE3039_B2_VS_B1_BCELL_DN1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("GSE3039_B2_VS_B1_BCELL_DN1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

##Marginal Zone B cells
GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_UP <- list(Book1$GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_UP)
GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_UP <- lapply(GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_UP, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_UP, name = "GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_UP", n = 22, search = FALSE)
VlnPlot(B_cells, c("GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_UP1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_UP1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

GSE24972_MARGINAL_ZONE_BCELL_VS_FOLLICULAR_BCELL_UP <- list(Book1$GSE24972_MARGINAL_ZONE_BCELL_VS_FOLLICULAR_BCELL_UP)
GSE24972_MARGINAL_ZONE_BCELL_VS_FOLLICULAR_BCELL_UP <- lapply(GSE24972_MARGINAL_ZONE_BCELL_VS_FOLLICULAR_BCELL_UP, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE24972_MARGINAL_ZONE_BCELL_VS_FOLLICULAR_BCELL_UP, name = "GSE24972_MARGINAL_ZONE_BCELL_VS_FOLLICULAR_BCELL_UP", n = 22, search = FALSE)
VlnPlot(B_cells, c("GSE24972_MARGINAL_ZONE_BCELL_VS_FOLLICULAR_BCELL_UP1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("GSE24972_MARGINAL_ZONE_BCELL_VS_FOLLICULAR_BCELL_UP1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

MZ_Bcells <- list(Book1$MZ_Bcells)
MZ_Bcells <- lapply(MZ_Bcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = MZ_Bcells, name = "MZ_Bcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("MZ_Bcells1"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(B_cells, c("MZ_Bcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

##Germinal Centre B cells
GSE23925_LIGHT_ZONE_VS_NAIVE_BCELL_UP <- list(Book1$GSE23925_LIGHT_ZONE_VS_NAIVE_BCELL_UP)
GSE23925_LIGHT_ZONE_VS_NAIVE_BCELL_UP <- lapply(GSE23925_LIGHT_ZONE_VS_NAIVE_BCELL_UP, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE23925_LIGHT_ZONE_VS_NAIVE_BCELL_UP, name = "GSE23925_LIGHT_ZONE_VS_NAIVE_BCELL_UP", n = 22, search = FALSE)
VlnPlot(B_cells, c("GSE23925_LIGHT_ZONE_VS_NAIVE_BCELL_UP1"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(B_cells, c("GSE23925_LIGHT_ZONE_VS_NAIVE_BCELL_UP1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

GSE23925_DARK_ZONE_VS_NAIVE_BCELL_UP <- list(Book1$GSE23925_DARK_ZONE_VS_NAIVE_BCELL_UP)
GSE23925_DARK_ZONE_VS_NAIVE_BCELL_UP <- lapply(GSE23925_DARK_ZONE_VS_NAIVE_BCELL_UP, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = GSE23925_DARK_ZONE_VS_NAIVE_BCELL_UP, name = "GSE23925_DARK_ZONE_VS_NAIVE_BCELL_UP", n = 22, search = FALSE)
VlnPlot(B_cells, c("GSE23925_DARK_ZONE_VS_NAIVE_BCELL_UP1"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(B_cells, c("GSE23925_DARK_ZONE_VS_NAIVE_BCELL_UP1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

DZ_a <- list(Book1$DZ_a)
DZ_a <- lapply(DZ_a, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = DZ_a, name = "DZ_a", n = 22, search = FALSE)
Dz_a1 <- VlnPlot(B_cells, c("DZ_a1"), pt.size  = 0.5, cols = turbo(17)) + ggtitle("DZ B cells") + NoLegend()
DZ_a1 <- FeaturePlot(B_cells, c("DZ_a1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + ggtitle("DZ B cells")

LZ_a <- list(Book1$LZ_a)
LZ_a <- lapply(LZ_a, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = LZ_a, name = "LZ_a", n = 22, search = FALSE)
Lz_a1 <- VlnPlot(B_cells, c("LZ_a1"), pt.size  = 0.5, cols = turbo(17)) + ggtitle("LZ B cells") + NoLegend()
LZ_a1 <- FeaturePlot(B_cells, c("LZ_a1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + ggtitle("LZ B cells")

LZ_b <- list(Book1$LZ_b)
LZ_b <- lapply(LZ_b, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = LZ_b, name = "LZ_b", n = 22, search = FALSE)
Lz_b <- VlnPlot(B_cells, c("LZ_b1"), pt.size  = 0.5, cols = turbo(17)) + ggtitle("LZ B cells") + NoLegend()
LZ_b <- FeaturePlot(B_cells, c("LZ_b1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + ggtitle("LZ B cells")

Int_a <- list(Book1$Int_a)
Int_a <- lapply(Int_a, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Int_a, name = "Int_a", n = 22, search = FALSE)
Int_a <- VlnPlot(B_cells, c("Int_a1"), pt.size  = 0.5, cols = turbo(17)) + ggtitle("Int_a") + NoLegend()
INT_a <- FeaturePlot(B_cells, c("Int_a1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + ggtitle("Int_a")

Int_c <- list(Book1$Int_c)
Int_c <- lapply(Int_c, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Int_c, name = "Int_c", n = 22, search = FALSE)
Int_c <- VlnPlot(B_cells, c("Int_c1"), pt.size  = 0.5, cols = turbo(17)) + ggtitle("Int_c") + NoLegend()
INT_c <- FeaturePlot(B_cells, c("Int_c1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + ggtitle("Int_c")

##Heatmap for B cell signatures
Data<-B_cells@meta.data
Forest <-Data %>% 
  group_by(`seurat_clusters`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$seurat_clusters
Forest<-Forest[ , c(28:35, 49:53)]
final<-t(Forest)

Data<-B_cells@meta.data
Forest <-Data %>% 
  group_by(`old.ident`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$old.ident
Forest<-Forest[ , c(28:35)]
final<-t(Forest)


pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)

##Cytokine responses from ImmGen - genes >= 2.0 log2fold change
#IFNalpha response - genes >= 2.0 log2fold change
IFNa_response <- list(Book1$IFNa_response)
IFNa_response <- lapply(IFNa_response, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IFNa_response, name = "IFNa_response", n = 22, search = FALSE)
VlnPlot(B_cells, c("IFNa_response1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IFNg_response1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

#IL-2 response - genes >= 1.0 log2fold change
IL2_response_FOBcells <- list(Book1$IL2_response_FOBcells)
IL2_response_FOBcells <- lapply(IL2_response_FOBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL2_response_FOBcells, name = "IL2_response_FOBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL2_response_FOBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL2_response_FOBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)


IL2_response_MZBcells <- list(Book1$IL2_response_MZBcells)
IL2_response_MZBcells <- lapply(IL2_response_MZBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL2_response_MZBcells, name = "IL2_response_MZBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL2_response_MZBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL2_response_MZBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

#IL-4 response - genes >= 2.0 log2fold change
IL4_response_FOBcells <- list(Book1$IL4_response_FOBcells)
IL4_response_FOBcells <- lapply(IL4_response_FOBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL4_response_FOBcells, name = "IL4_response_FOBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL4_response_FOBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL4_response_FOBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

IL4_response_MZBcells <- list(Book1$IL4_response_MZBcells)
IL4_response_MZBcells <- lapply(IL4_response_MZBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL4_response_MZBcells, name = "IL4_response_MZBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL4_response_MZBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL4_response_MZBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

#IL-7 response - genes >= 1.0 log2fold change
IL7_response_FOBcells <- list(Book1$IL7_response_FOBcells)
IL7_response_FOBcells <- lapply(IL7_response_FOBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL7_response_FOBcells, name = "IL7_response_FOBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL7_response_FOBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL7_response_FOBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

IL7_response_MZBcells <- list(Book1$IL7_response_MZBcells)
IL7_response_MZBcells <- lapply(IL7_response_MZBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL7_response_MZBcells, name = "IL7_response_MZBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL7_response_MZBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL7_response_MZBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

#IL-9 response - genes >= 1.0 log2fold change
IL9_response_FOBcells <- list(Book1$IL9_response_FOBcells)
IL9_response_FOBcells <- lapply(IL9_response_FOBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL9_response_FOBcells, name = "IL9_response_FOBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL9_response_FOBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL9_response_FOBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

IL9_response_MZBcells <- list(Book1$IL9_response_MZBcells)
IL9_response_MZBcells <- lapply(IL9_response_MZBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL9_response_MZBcells, name = "IL9_response_MZBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL9_response_MZBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL9_response_MZBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

#IL-15 response - genes >= 1.0 log2fold change
IL15_response_FOBcells <- list(Book1$IL15_response_FOBcells)
IL15_response_FOBcells <- lapply(IL15_response_FOBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL15_response_FOBcells, name = "IL15_response_FOBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL15_response_FOBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL15_response_FOBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

IL15_response_MZBcells <- list(Book1$IL15_response_MZBcells)
IL15_response_MZBcells <- lapply(IL15_response_MZBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL15_response_MZBcells, name = "IL15_response_MZBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL15_response_MZBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL15_response_MZBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

#IL-21 response - genes >= 1.0 log2fold change
IL21_response_FOBcells <- list(Book1$IL21_response_FOBcells)
IL21_response_FOBcells <- lapply(IL21_response_FOBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL21_response_FOBcells, name = "IL21_response_FOBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL21_response_FOBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL21_response_FOBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

IL21_response_MZBcells <- list(Book1$IL21_response_MZBcells)
IL21_response_MZBcells <- lapply(IL21_response_MZBcells, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = IL21_response_MZBcells, name = "IL21_response_MZBcells", n = 22, search = FALSE)
VlnPlot(B_cells, c("IL21_response_MZBcells1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("IL21_response_MZBcells1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

B_cells_p6

#IFN
HALLMARK_INTERFERON_ALPHA_RESPONSE <- list(Book1$HALLMARK_INTERFERON_ALPHA_RESPONSE)
HALLMARK_INTERFERON_ALPHA_RESPONSE <- lapply(HALLMARK_INTERFERON_ALPHA_RESPONSE, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_INTERFERON_ALPHA_RESPONSE, name = "HALLMARK_INTERFERON_ALPHA_RESPONSE", n = 22, search = FALSE)
VlnPlot(B_cells, c("HALLMARK_INTERFERON_ALPHA_RESPONSE1"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(B_cells, c("HALLMARK_INTERFERON_ALPHA_RESPONSE1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

FeaturePlot(B_cells, features = "Myc", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2)
VlnPlot(B_cells, features = "Myc", pt.size  = 0.5, cols = turbo(14)) + NoLegend()

FeaturePlot(B_cells, features = "Cxcr4", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2)
VlnPlot(B_cells, features = "Cxcr4", pt.size  = 0.5, cols = turbo(14)) + NoLegend()

FeaturePlot(B_cells, features = "Cd86", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2)
VlnPlot(B_cells, features = "Cd86", pt.size  = 0.5, cols = turbo(14)) + NoLegend()

FeaturePlot(B_cells, features = "Irf4", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2)
VlnPlot(B_cells, features = "Irf4", pt.size  = 0.5, cols = turbo(14)) + NoLegend()

FeaturePlot(B_cells, features = "Cd40", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2)
VlnPlot(B_cells, features = "Cd40", pt.size  = 0.5, cols = turbo(14)) + NoLegend()

FeaturePlot(B_cells, features = "Foxo1", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2)
VlnPlot(B_cells, features = "Foxo1", pt.size  = 0.5, cols = turbo(14)) + NoLegend()

FeaturePlot(B_cells, features = "Fcer1g", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.2)
VlnPlot(B_cells, features = "Foxo1", pt.size  = 0.5, cols = turbo(14)) + NoLegend()

#Heatmap for cytokine signatures
Data<-B_cells@meta.data
Forest <-Data %>% 
  group_by(`seurat_clusters`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$seurat_clusters
Forest<-Forest[ , c(36:48)]
final<-t(Forest)

Data<-B_cells@meta.data
Forest <-Data %>% 
  group_by(`old.ident`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$old.ident
Forest<-Forest[ , c(36:48)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)
##Inflammation
HALLMARK_INFLAMMATORY_RESPONSE <- list(Book1$HALLMARK_INFLAMMATORY_RESPONSE)
HALLMARK_INFLAMMATORY_RESPONSE <- lapply(HALLMARK_INFLAMMATORY_RESPONSE, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_INFLAMMATORY_RESPONSE, name = "HALLMARK_INFLAMMATORY_RESPONSE", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = HALLMARK_INFLAMMATORY_RESPONSE, name = "HALLMARK_INFLAMMATORY_RESPONSE", n = 22, search = TRUE)
VlnPlot(B_cells, c("HALLMARK_INFLAMMATORY_RESPONSE1"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(B_cells, c("HALLMARK_INFLAMMATORY_RESPONSE1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("HALLMARK_INFLAMMATORY_RESPONSE1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(Cluster_3, c("HALLMARK_INFLAMMATORY_RESPONSE1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

VlnPlot(B_cells, c("Cd19"), pt.size  = 0.5, cols = turbo(14))

##IL6-JAK-STAT3 signalling
HALLMARK_IL6_JAK_STAT3_SIGNALING <- list(Book1$HALLMARK_IL6_JAK_STAT3_SIGNALING)
HALLMARK_IL6_JAK_STAT3_SIGNALING <- lapply(HALLMARK_IL6_JAK_STAT3_SIGNALING, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_IL6_JAK_STAT3_SIGNALING, name = "HALLMARK_IL6_JAK_STAT3_SIGNALING", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = HALLMARK_IL6_JAK_STAT3_SIGNALING, name = "HALLMARK_IL6_JAK_STAT3_SIGNALING", n = 22, search = TRUE)
VlnPlot(B_cells, c("HALLMARK_IL6_JAK_STAT3_SIGNALING1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_IL6_JAK_STAT3_SIGNALING1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("HALLMARK_IL6_JAK_STAT3_SIGNALING1"), pt.size  = 0.5, cols = turbo(14))

##Myc targets
HALLMARK_MYC_TARGETS_V1 <- list(Book1$HALLMARK_MYC_TARGETS_V1)
HALLMARK_MYC_TARGETS_V1 <- lapply(HALLMARK_MYC_TARGETS_V1, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_MYC_TARGETS_V1, name = "HALLMARK_MYC_TARGETS_V1", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = HALLMARK_MYC_TARGETS_V1, name = "HALLMARK_MYC_TARGETS_V1", n = 22, search = FALSE)
VlnPlot(B_cells, c("HALLMARK_MYC_TARGETS_V11"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_MYC_TARGETS_V11"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("HALLMARK_MYC_TARGETS_V11"), pt.size  = 0.5, cols = turbo(14))

HALLMARK_MYC_TARGETS_V2 <- list(Book1$HALLMARK_MYC_TARGETS_V2)
HALLMARK_MYC_TARGETS_V2 <- lapply(HALLMARK_MYC_TARGETS_V2, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_MYC_TARGETS_V2, name = "HALLMARK_MYC_TARGETS_V2", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = HALLMARK_MYC_TARGETS_V2, name = "HALLMARK_MYC_TARGETS_V2", n = 22, search = FALSE)
VlnPlot(B_cells, c("HALLMARK_MYC_TARGETS_V21"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_MYC_TARGETS_V21"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("HALLMARK_MYC_TARGETS_V21"), pt.size  = 0.5, cols = turbo(14))

##TNFA signalling
HALLMARK_TNFA_SIGNALING_VIA_NFKB <- list(Book1$HALLMARK_TNFA_SIGNALING_VIA_NFKB)
HALLMARK_TNFA_SIGNALING_VIA_NFKB <- lapply(HALLMARK_TNFA_SIGNALING_VIA_NFKB, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_TNFA_SIGNALING_VIA_NFKB, name = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = HALLMARK_TNFA_SIGNALING_VIA_NFKB, name = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", n = 22, search = FALSE)
VlnPlot(B_cells, c("HALLMARK_TNFA_SIGNALING_VIA_NFKB1"), pt.size  = 0.5, cols = turbo(17))
FeaturePlot(B_cells, c("HALLMARK_TNFA_SIGNALING_VIA_NFKB1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("HALLMARK_TNFA_SIGNALING_VIA_NFKB1"), pt.size  = 0.5, cols = turbo(14))

##WNT signalling pathways
WNT_signalling <- list(Book1$WNT_signalling)
WNT_signalling <- lapply(WNT_signalling, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = WNT_signalling, name = "WNT_signalling", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = WNT_signalling, name = "WNT_signalling", n = 22, search = TRUE)
VlnPlot(B_cells, c("WNT_signalling1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("WNT_signalling1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("WNT_signalling1"), pt.size  = 0.5, cols = turbo(14))

HALLMARK_WNT_BETA_CATENIN_SIGNALING <- list(Book1$HALLMARK_WNT_BETA_CATENIN_SIGNALING)
HALLMARK_WNT_BETA_CATENIN_SIGNALING <- lapply(HALLMARK_WNT_BETA_CATENIN_SIGNALING, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_WNT_BETA_CATENIN_SIGNALING, name = "HALLMARK_WNT_BETA_CATENIN_SIGNALING", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = HALLMARK_WNT_BETA_CATENIN_SIGNALING, name = "HALLMARK_WNT_BETA_CATENIN_SIGNALING", n = 22, search = TRUE)
VlnPlot(B_cells, c("HALLMARK_WNT_BETA_CATENIN_SIGNALING1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_WNT_BETA_CATENIN_SIGNALING1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("HALLMARK_WNT_BETA_CATENIN_SIGNALING1"), pt.size  = 0.5, cols = turbo(14))

##Hedgehog signalling pathways
HALLMARK_HEDGEHOG_SIGNALING <- list(Book1$HALLMARK_HEDGEHOG_SIGNALING)
HALLMARK_HEDGEHOG_SIGNALING <- lapply(HALLMARK_HEDGEHOG_SIGNALING, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_HEDGEHOG_SIGNALING, name = "HALLMARK_HEDGEHOG_SIGNALING", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = HALLMARK_HEDGEHOG_SIGNALING, name = "HALLMARK_HEDGEHOG_SIGNALING", n = 22, search = TRUE)
VlnPlot(B_cells, c("HALLMARK_HEDGEHOG_SIGNALING1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_HEDGEHOG_SIGNALING1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("HALLMARK_HEDGEHOG_SIGNALING1"), pt.size  = 0.5, cols = turbo(14))

##Cholesterol Metabolism
HALLMARK_CHOLESTEROL_HOMEOSTASIS <- list(Book1$HALLMARK_CHOLESTEROL_HOMEOSTASIS)
HALLMARK_CHOLESTEROL_HOMEOSTASIS <- lapply(HALLMARK_CHOLESTEROL_HOMEOSTASIS, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_CHOLESTEROL_HOMEOSTASIS, name = "HALLMARK_CHOLESTEROL_HOMEOSTASIS", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = HALLMARK_CHOLESTEROL_HOMEOSTASIS, name = "HALLMARK_CHOLESTEROL_HOMEOSTASIS", n = 22, search = FALSE)
VlnPlot(B_cells, c("HALLMARK_CHOLESTEROL_HOMEOSTASIS1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_CHOLESTEROL_HOMEOSTASIS1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("HALLMARK_CHOLESTEROL_HOMEOSTASIS1"), pt.size  = 0.5, cols = turbo(14))

Cholsterol_metabolism <- list(Book1$Cholsterol_metabolism)
Cholsterol_metabolism <- lapply(Cholsterol_metabolism, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Cholsterol_metabolism, name = "Cholsterol_metabolism", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Cholsterol_metabolism, name = "Cholsterol_metabolism", n = 22, search = TRUE)
VlnPlot(B_cells, c("Cholsterol_metabolism1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("Cholsterol_metabolism1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("Cholsterol_metabolism1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(Cluster_3, c("Cholsterol_metabolism1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

##Glucocorticoid (Nr3c1)-signalling
Nr3c1_signalling <- list(Book1$Nr3c1_signalling)
Nr3c1_signalling <- lapply(Nr3c1_signalling, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Nr3c1_signalling, name = "Nr3c1_signalling", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Nr3c1_signalling, name = "Nr3c1_signalling", n = 22, search = TRUE)
VlnPlot(B_cells, c("Nr3c1_signalling1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("Nr3c1_signalling1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("Nr3c1_signalling1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(Cluster_3, c("Nr3c1_signalling1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, features = "Cyp19a1", pt.size  = 0.5, cols = turbo(14))

##Metabolism
HALLMARK_FATTY_ACID_METABOLISM <- list(Book1$HALLMARK_FATTY_ACID_METABOLISM)
HALLMARK_FATTY_ACID_METABOLISM <- lapply(HALLMARK_FATTY_ACID_METABOLISM, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_FATTY_ACID_METABOLISM, name = "HALLMARK_FATTY_ACID_METABOLISM", n = 22, search = TRUE)
VlnPlot(B_cells, c("HALLMARK_FATTY_ACID_METABOLISM1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_FATTY_ACID_METABOLISM1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

Lipid_metabolism <- list(Book1$Lipid_metabolism)
Lipid_metabolism <- lapply(Lipid_metabolism, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Lipid_metabolism, name = "Lipid_metabolism", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Lipid_metabolism, name = "Lipid_metabolism", n = 22, search = TRUE)
VlnPlot(B_cells, c("Lipid_metabolism1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("Lipid_metabolism1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("Lipid_metabolism1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(Cluster_3, c("Lipid_metabolism1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

HALLMARK_GLYCOLYSIS <- list(Book1$HALLMARK_GLYCOLYSIS)
HALLMARK_GLYCOLYSIS <- lapply(HALLMARK_GLYCOLYSIS, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_GLYCOLYSIS, name = "HALLMARK_GLYCOLYSIS", n = 22, search = TRUE)
VlnPlot(B_cells, c("HALLMARK_GLYCOLYSIS1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_GLYCOLYSIS1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

HALLMARK_OXIDATIVE_PHOSPHORYLATION <- list(Book1$HALLMARK_OXIDATIVE_PHOSPHORYLATION)
HALLMARK_OXIDATIVE_PHOSPHORYLATION <- lapply(HALLMARK_OXIDATIVE_PHOSPHORYLATION, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_OXIDATIVE_PHOSPHORYLATION, name = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", n = 22, search = FALSE)
VlnPlot(B_cells, c("HALLMARK_OXIDATIVE_PHOSPHORYLATION1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_OXIDATIVE_PHOSPHORYLATION1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

HALLMARK_PI3K_AKT_MTOR_SIGNALING <- list(Book1$HALLMARK_PI3K_AKT_MTOR_SIGNALING)
HALLMARK_PI3K_AKT_MTOR_SIGNALING <- lapply(HALLMARK_PI3K_AKT_MTOR_SIGNALING, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = HALLMARK_PI3K_AKT_MTOR_SIGNALING, name = "HALLMARK_PI3K_AKT_MTOR_SIGNALING", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = HALLMARK_PI3K_AKT_MTOR_SIGNALING, name = "HALLMARK_PI3K_AKT_MTOR_SIGNALING", n = 22, search = FALSE)
VlnPlot(B_cells, c("HALLMARK_PI3K_AKT_MTOR_SIGNALING1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("HALLMARK_PI3K_AKT_MTOR_SIGNALING1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)
VlnPlot(Cluster_3, c("HALLMARK_PI3K_AKT_MTOR_SIGNALING1"), pt.size  = 0.5, cols = turbo(14))

##Steroidogenesis
Steroidogenesis <- list(Book1$Steroidogenesis)
Steroidogenesis <- lapply(Steroidogenesis, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Steroidogenesis, name = "Steroidogenesis", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Steroidogenesis, name = "Steroidogenesis", n = 22, search = FALSE)
VlnPlot(B_cells, c("Steroidogenesis1"), pt.size  = 0.5, cols = turbo(14)) + NoLegend() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(B_cells, c("Steroidogenesis1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
VlnPlot(Cluster_3, c("Steroidogenesis1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(Cluster_3, c("Steroidogenesis1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))

Steroidogenesis_mouse <- list(Book1$Steroidogenesis_mouse)
Steroidogenesis_mouse <- lapply(Steroidogenesis_mouse, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Steroidogenesis_mouse, name = "Steroidogenesis_mouse", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Steroidogenesis_mouse, name = "Steroidogenesis_mouse", n = 22, search = FALSE)
VlnPlot(B_cells, c("Steroidogenesis_mouse1"), pt.size  = 0.5, cols = turbo(17)) + NoLegend() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(B_cells, c("Steroidogenesis_mouse1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
VlnPlot(Cluster_3, c("Steroidogenesis_mouse1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(Cluster_3, c("Steroidogenesis_mouse1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))

Nr3c1_supressed <- list(Book1$Nr3c1_supressed)
Nr3c1_supressed <- lapply(Nr3c1_supressed, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Nr3c1_supressed, name = "Nr3c1_supressed", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Nr3c1_supressed, name = "Nr3c1_supressed", n = 22, search = FALSE)
VlnPlot(B_cells, c("Nr3c1_supressed1"), pt.size  = 0.5, cols = turbo(14)) + NoLegend() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(B_cells, c("Nr3c1_supressed1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
VlnPlot(Cluster_3, c("Nr3c1_supressed1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(Cluster_3, c("Nr3c1_supressed1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))

StAR_proteins <- list(Book1$StAR_proteins)
StAR_proteins <- lapply(StAR_proteins, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = StAR_proteins, name = "StAR_proteins", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = StAR_proteins, name = "StAR_proteins", n = 22, search = FALSE)
VlnPlot(B_cells, c("StAR_proteins1"), pt.size  = 0.5, cols = turbo(14)) + NoLegend() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(B_cells, c("StAR_proteins1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
VlnPlot(Cluster_3, c("StAR_proteins1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(Cluster_3, c("StAR_proteins1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))

Nr3c2_associated_genes <- list(Book1$Nr3c2_associated_genes)
Nr3c2_associated_genes <- lapply(Nr3c2_associated_genes, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Nr3c2_associated_genes, name = "Nr3c2_associated_genes", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Nr3c2_associated_genes, name = "Nr3c2_associated_genes", n = 22, search = FALSE)
VlnPlot(B_cells, c("Nr3c2_associated_genes1"), pt.size  = 0.5, cols = turbo(14)) + NoLegend() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(B_cells, c("Nr3c2_associated_genes1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
VlnPlot(Cluster_3, c("Nr3c2_associated_genes1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(Cluster_3, c("Nr3c2_associated_genes1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))

Nr3c2_associated_genes2 <- list(Book1$Nr3c2_associated_genes2)
Nr3c2_associated_genes2 <- lapply(Nr3c2_associated_genes2, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Nr3c2_associated_genes2, name = "Nr3c2_associated_genes2", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Nr3c2_associated_genes2, name = "Nr3c2_associated_genes2", n = 22, search = FALSE)
Banana <- VlnPlot(B_cells, c("Nr3c2_associated_genes21"), pt.size  = 0.5, cols = turbo(14)) + NoLegend() + ggtitle("Enrichment for Nr3c2-associated genes") + NoLegend()+ theme(plot.title = element_text(color="black", size=15, face="bold")) + xlab("Cluster")
Apple <- FeaturePlot(B_cells, c("Nr3c2_associated_genes21"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for Nr3c2-associated genes") + theme(plot.title = element_text(color="black", size=15, face="bold"))
Grape <- VlnPlot(Cluster_3, c("Nr3c2_associated_genes21"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Enrichment for Nr3c2-associated genes") + NoLegend()+ theme(plot.title = element_text(color="black", size=15, face="bold")) + xlab("Cyp11a1")
grid.arrange(Apple, Banana, Grape, ncol = 3)

Cyp11a1_associated_genes <- list(Book1$Cyp11a1_associated_genes)
Cyp11a1_associated_genes <- lapply(Cyp11a1_associated_genes, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Cyp11a1_associated_genes, name = "Cyp11a1_associated_genes", n = 22, search = FALSE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Cyp11a1_associated_genes, name = "Cyp11a1_associated_genes", n = 22, search = FALSE)
VlnPlot(B_cells, c("Cyp11a1_associated_genes1"), pt.size  = 0.5, cols = turbo(14)) + NoLegend() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(B_cells, c("Cyp11a1_associated_genes1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))
VlnPlot(Cluster_3, c("Cyp11a1_associated_genes1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Enrichment for genes in the steroidogenesis pathway") + theme(plot.title = element_text(color="black", size=15, face="bold"))

Nr3c1_plot1 <- FeaturePlot(B_cells, c("Nr3c1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Nr3c1") + theme(plot.title = element_text(color="black", size=15, face="bold"))
Nr3c1_plot2 <- VlnPlot(B_cells, c("Nr3c1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Nr3c1") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold")) + xlab("Cluster")
Nr3c1_plot3 <- VlnPlot(Cluster_3, c("Nr3c1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Nr3c1 expression in Cluster 3") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold")) + xlab("Cyp11a1")
grid.arrange(Nr3c1_plot1, Nr3c1_plot2, Nr3c1_plot3,  ncol = 3)

Nr3c2_plot1 <- FeaturePlot(B_cells, c("Nr3c2"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("Nr3c2") + theme(plot.title = element_text(color="black", size=15, face="bold"))
Nr3c2_plot2 <- VlnPlot(B_cells, c("Nr3c2"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Nr3c2") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold")) + xlab("Cluster")
Nr3c2_plot3 <- VlnPlot(Cluster_3, c("Nr3c2"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("Nr3c2 expression in Cluster 3") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold")) + xlab("Cyp11a1")
grid.arrange(Nr3c2_plot1, Nr3c2_plot2, Nr3c2_plot3,  ncol = 3)

##Cyp11a1 co-expressed and interactors
Cyp11a1_co_expressed <- list(Book1$Cyp11a1_co_expressed)
Cyp11a1_co_expressed <- lapply(Cyp11a1_co_expressed, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Cyp11a1_co_expressed, name = "Cyp11a1_co_expressed", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Cyp11a1_co_expressed, name = "Cyp11a1_co_expressed", n = 22, search = FALSE)
VlnPlot(B_cells, c("Cyp11a1_co_expressed1"), pt.size  = 0.5, cols = turbo(14)) + NoLegend() + ggtitle("") + NoLegend() + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(B_cells, c("Cyp11a1_co_expressed1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("") + theme(plot.title = element_text(color="black", size=15, face="bold"))
VlnPlot(Cluster_3, c("Cyp11a1_co_expressed1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("") + NoLegend()+ theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(Cluster_3, c("Cyp11a1_co_expressed1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("") + theme(plot.title = element_text(color="black", size=15, face="bold"))

Cyp11a1_interactors <- list(Book1$Cyp11a1_interactors)
Cyp11a1_interactors <- lapply(Cyp11a1_interactors, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Cyp11a1_interactors, name = "Cyp11a1_interactors", n = 22, search = TRUE)
Cluster_3 <-AddModuleScore(Cluster_3, features = Cyp11a1_interactors, name = "Cyp11a1_interactors", n = 22, search = FALSE)
VlnPlot(B_cells, c("Cyp11a1_interactors1"), pt.size  = 0.5, cols = turbo(14)) + NoLegend() + ggtitle("") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(B_cells, c("Cyp11a1_interactors1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("") + theme(plot.title = element_text(color="black", size=15, face="bold"))
VlnPlot(Cluster_3, c("Cyp11a1_interactors1"), pt.size  = 0.5, cols = turbo(14)) + ggtitle("") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(Cluster_3, c("Cyp11a1_interactors1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("") + theme(plot.title = element_text(color="black", size=15, face="bold"))

##Unknown cluster
Unknown_cluster <- list(Book1$Unknown_cluster)
Unknown_cluster <- lapply(Unknown_cluster, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Unknown_cluster, name = "Unknown_cluster", n = 22, search = TRUE)
VlnPlot(B_cells, c("Unknown_cluster1"), pt.size  = 0.5, cols = turbo(17)) + NoLegend() + ggtitle("") + theme(plot.title = element_text(color="black", size=15, face="bold"))
FeaturePlot(B_cells, c("Unknown_cluster1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5) + theme_void() + ggtitle("") + theme(plot.title = element_text(color="black", size=15, face="bold"))


####Annotation by reference datasets####
#Load reference
load(file.choose("NicheData10x.rda"), verbose = TRUE)
head(NicheData10x[[]])
colnames(NicheData10x[[]])
DimPlot(NicheData10x)
FeaturePlot(NicheData10x, features = "Cd19")
reference <- NicheData10x
reference[["old.ident"]] <- Idents(reference)
Idents(reference) <- levels(reference)

#DEG list from clusters of interest: B cell, large-pre B, pro-B, small pre-B
Bcells_ref <- FindMarkers(reference, ident.1 = "B cell", assay = "RNA")
largepreB_ref <- FindMarkers(reference, ident.1 = "large pre-B.", assay = "RNA")
proB_ref <- FindMarkers(reference, ident.1 = "pro-B", assay = "RNA")
smallpreB_ref <- FindMarkers(reference, ident.1 = "small pre-B.", assay = "RNA")

#Export as Excel file and add to book.1 for GSEA
Bcells_ref <- tibble::rownames_to_column(Bcells_ref, "Genes")
write_xlsx(Bcells_ref, "Bcells_ref.xlsx")

largepreB_ref <- tibble::rownames_to_column(largepreB_ref, "Genes")
write_xlsx(largepreB_ref, "largepreB_ref.xlsx")

proB_ref <- tibble::rownames_to_column(proB_ref, "Genes")
write_xlsx(proB_ref, "proB_ref.xlsx")

smallpreB_ref <- tibble::rownames_to_column(smallpreB_ref, "Genes")
write_xlsx(smallpreB_ref, "smallpreB_ref.xlsx")

#GSEA
Book1 <- read.csv("Book1.csv",header = T, sep = ',')

#Test on reference
Bcell_ref <- list(Book1$Bcell_ref)
Bcell_ref <- lapply(Bcell_ref, str_to_title)
DefaultAssay(reference) <- "RNA"
reference <-AddModuleScore(reference, features = Bcell_ref, name = "Bcell_ref", n = 22, search = FALSE)
VlnPlot(reference, c("Bcell_ref1"), pt.size  = 0.5, cols = turbo(32))
FeaturePlot(reference, c("Bcell_ref1"), cols=col_con, pt.size  = 1.5)
DimPlot(reference, label = TRUE)

#Test on B_cell dataset
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = Bcell_ref, name = "Bcell_ref", n = 22, search = FALSE)
VlnPlot(B_cells, c("Bcell_ref1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("Bcell_ref1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

largepreB_ref<- list(Book1$largepreB_ref)
largepreB_ref <- lapply(largepreB_ref, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = largepreB_ref, name = "largepreB_ref", n = 22, search = FALSE)
VlnPlot(B_cells, c("largepreB_ref1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("largepreB_ref1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

proB_ref <- list(Book1$proB_ref)
proB_ref <- lapply(proB_ref, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = proB_ref, name = "proB_ref", n = 22, search = FALSE)
VlnPlot(B_cells, c("proB_ref1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("proB_ref1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

smallpreB_ref <- list(Book1$smallpreB_ref)
smallpreB_ref <- lapply(smallpreB_ref, str_to_title)
DefaultAssay(B_cells) <- "RNA"
B_cells <-AddModuleScore(B_cells, features = smallpreB_ref, name = "smallpreB_ref", n = 22, search = FALSE)
VlnPlot(B_cells, c("smallpreB_ref1"), pt.size  = 0.5, cols = turbo(14))
FeaturePlot(B_cells, c("smallpreB_ref1"), cols=col_con, reduction = "harmony.wnn.umap", pt.size  = 1.5)

#Heatmap for early prgenitor B cell GSEA
Data<-B_cells@meta.data
colnames(Data)
Forest <-Data %>% 
  group_by(`seurat_clusters`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$seurat_clusters
Forest<-Forest[ , c(54:57)]
final<-t(Forest)

colnames(Data)

Data<-B_cells@meta.data
Forest <-Data %>% 
  group_by(`old.ident`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$old.ident
Forest<-Forest[ , c(54:57)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T,
                   cellwidth = 30,cellheight = 30, angle_col = 45, color = col_con)


Book1 <- read.csv("Book1.csv",header = T, sep = ',')
Idents(B_cells) <-  B_cells$seurat_clusters

#DEGs Cyp11a1 expressing cells vs all B cells
?FindMarkers
DEG_Cyp11a1_cells <- FindMarkers(B_cells, )
B_cells_1 <- B_cells
Cyp11a1_cells <- subset(B_cells_1, subset = Cyp11a1 >0)
cellNames <- rownames(Cyp11a1_cells@meta.data)
B_cells_1@meta.data <- B_cells_1@meta.data %>% mutate(Cyp11a1 = ifelse((B_cells_1$barcode %in% cellNames), "pos", "neg"))
DimPlot(B_cells_1, reduction = "harmony.wnn.umap", group.by = "Cyp11a1")
Idents(B_cells_1) <- B_cells_1@meta.data$Cyp11a1
Cyp11a1_DEGs <- FindMarkers(B_cells_1, ident.1 = "pos", ident.2 = "neg", assay = "RNA")
Cyp11a1_DEGs_Abs <- FindMarkers(B_cells_1, ident.1 = "pos", ident.2 = "neg", assay = "ADT")
Cyp11a1_DEGs_Abs <- tibble::rownames_to_column(Cyp11a1_DEGs_Abs, "Genes")

DefaultAssay(B_cells_1)<-"RNA"
Idents(B_cells_1) <-  B_cells_1$Cyp11a1
Unique_Cyp11a1 <- WhichCells(object = B_cells_1, ident = "pos")
Unique_Cyp11a1 <- tfidf(GetAssayData(B_cells_1), Unique_Cyp11a1, colnames(B_cells_1))
Unique_Cyp11a1 <- tibble::rownames_to_column(Unique_Cyp11a1, "Genes")
write_xlsx(Unique_Cyp11a1, "Unique_Cyp11a1.xlsx")

Cluster_3 <- subset(B_cells_1, subset = seurat_clusters == "Cluster 3")
Idents(Cluster_3) <- Cluster_3$Cyp11a1
DimPlot(Cluster_3, reduction = "harmony.wnn.umap")
Cyp11a1_DEGs_c3 <- FindMarkers(Cluster_3, ident.1 = "pos", ident.2 = "neg", assay = "RNA")
Cyp11a1_DEGs_c3 <- tibble::rownames_to_column(Cyp11a1_DEGs_c3, "Genes1")
write_xlsx(Cyp11a1_DEGs_c3, "Cyp11a1_DEGs_c3.xlsx")
Cyp11a1_DEGs_Abs_c3 <- FindMarkers(B_cells_1, ident.1 = "pos", ident.2 = "neg", assay = "ADT")
Cyp11a1_DEGs_Abs_c3 <- tibble::rownames_to_column(Cyp11a1_DEGs_Abs_c3, "Genes")

Cyp11a1_DEGs_c3_sig <- Cyp11a1_DEGs_c3 %>%
  filter(p_val_adj < 0.05)
Cyp11a1_DEGs_c3_sig <- tibble::rownames_to_column(Cyp11a1_DEGs_c3_sig, "Genes")
write_xlsx(Cyp11a1_DEGs_c3_sig, "Cyp11a1_DEGs_c3_sig.xlsx")

FeaturePlot(B_cells, features = "Cd200r3", reduction = "harmony.wnn.umap")
FeaturePlot(Cluster_3, features = "Cd200r3", reduction = "harmony.wnn.umap")

VlnPlot(B_cells, features = "Osm")
VlnPlot(Cluster_3, features = "Nr3c2")
FeaturePlot(Cyp11a1_cells, features = "Lrp6", reduction = "harmony.wnn.umap")
DimPlot(Cyp11a1_cells, group.by = "Genotype")
DimPlot(Cluster_3, group.by = "Genotype", reduction = "harmony.wnn.umap")

DefaultAssay(Cluster_3)<-"RNA"
Unique_Cyp11a1_c3 <- WhichCells(object = Cluster_3, ident = "pos")
Unique_Cyp11a1_c3 <- tfidf(GetAssayData(Cluster_3), Unique_Cyp11a1_c3, colnames(Cluster_3))
Unique_Cyp11a1_c3 <- tibble::rownames_to_column(Unique_Cyp11a1_c3, "Genes")
write_xlsx(Unique_Cyp11a1_c3, "Unique_Cyp11a1_c3.xlsx")

##ClusterProfiler analysis
#Cyp11a1 DEGs - Cyp11a1 expressing vs all cells
top100_Cyp11a1_DEGs <- Cyp11a1_DEGs %>% arrange(desc(avg_log2FC)) %>% top_n(n = 100, wt = avg_log2FC)
top100Cyp11a1_DEGs_pval <- subset(top100_Cyp11a1_DEGs, rowSums(top100_Cyp11a1_DEGs[6] < 0.05) > 0)
dim(top100Cyp11a1_DEGs_pval)
df <- top100Cyp11a1_DEGs_pval[,1:3]
df$Genes = bitr(df$Genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)

genelist2 <- list("0" = df$Genes$ENTREZID)
GOclusterplot <- compareCluster(geneClusters = genelist2, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
dotplot(GOclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist2, fun = "enrichKEGG", organism = "mmu")
dotplot(GOclusterplot) + dotplot(KEGGclusterplot)

#Cyp11a1 DEGs - Cyp11a1 expressing vs all cells in cluster 3
top100_Cyp11a1_DEGS_c3 <- Cyp11a1_DEGs_c3 %>% arrange(desc(avg_log2FC)) %>% top_n(n = 100, wt = avg_log2FC)
top100Cyp11a1_DEGs__c3_pval <- subset(top100_Cyp11a1_DEGS_c3, rowSums(top100_Cyp11a1_DEGS_c3[6] < 0.05) > 0)
dim(top100Cyp11a1_DEGs_pval)
df <- top100Cyp11a1_DEGs_pval[,1:3]
df$Genes = bitr(df$Genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)

genelist2 <- list("0" = df$Genes$ENTREZID)
GOclusterplot <- compareCluster(geneClusters = genelist2, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
dotplot(GOclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist2, fun = "enrichKEGG", organism = "mmu")
dotplot(GOclusterplot) + dotplot(KEGGclusterplot)

FeaturePlot(Cluster_3, features = "Il10", reduction = "harmony.wnn.umap")
FeaturePlot(Cluster_3, features = c("Spn", "Cyp11a1"), reduction = "harmony.wnn.umap")


Idents(B_cells_1) <- B_cells$seurat_clusters
Bcells.markers <- FindAllMarkers(B_cells_1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
Bcells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top100_all <- Bcells.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% arrange(desc(avg_log2FC))
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
dfsample$`Cluster 10` = bitr(dfsample$`Cluster 10`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 11` = bitr(dfsample$`Cluster 11`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 12` = bitr(dfsample$`Cluster 12`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)
dfsample$`Cluster 13` = bitr(dfsample$`Cluster 13`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", drop = FALSE)

genelist <- list("0" = dfsample$`Cluster 0`$ENTREZID,"1" = dfsample$`Cluster 1`$ENTREZID, "3" = dfsample$`Cluster 3`$ENTREZID, "4" = dfsample$`Cluster 4`$ENTREZID,
                 "5" = dfsample$`Cluster 5`$ENTREZID, "6" = dfsample$`Cluster 6`$ENTREZID, "7" = dfsample$`Cluster 7`$ENTREZID, "8" = dfsample$`Cluster 8`$ENTREZID, "9" = dfsample$`Cluster 9`$ENTREZID,
                 "10" = dfsample$`Cluster 10`$ENTREZID, "11" = dfsample$`Cluster 11`$ENTREZID, "12" = dfsample$`Cluster 12`$ENTREZID, "13" = dfsample$`Cluster 13`$ENTREZID)

GOclusterplot <- compareCluster(geneClusters = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
dotplot(GOclusterplot, font.size = 8, label_format = 100)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", organism = "mmu")
dotplot(KEGGclusterplot,  font.size = 8, label_format = 100)

####Clonotype analysis####
##Load contig files
c1_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/C1_VDJ/outs/filtered_contig_annotations.csv")
c2_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/C2_VDJ/outs/filtered_contig_annotations.csv")
d1_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/D1_VDJ/outs/filtered_contig_annotations.csv")
d2_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_1_files/VDJ/D2_VDJ/outs/filtered_contig_annotations.csv")
a_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/A_WT_VDJ/outs/filtered_contig_annotations.csv")
b_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/B_WT_VDJ/outs/filtered_contig_annotations.csv")
f_cl.data <- read.csv("~/Desktop/CITE-Sequencing_Data/CITE_Seq_2_files/VDJ_batch2/F_E1020K_VDJ/outs/filtered_contig_annotations.csv")

##Match barcode names with GE and ADT data
c1_cl.data$barcode=gsub("-1","_c1",c1_cl.data$barcode)
c2_cl.data$barcode=gsub("-1","_c2",c2_cl.data$barcode)
d1_cl.data$barcode=gsub("-1","_d1",d1_cl.data$barcode)
d2_cl.data$barcode=gsub("-1","_d2",d2_cl.data$barcode)
a_cl.data$barcode=gsub("-1","_a",a_cl.data$barcode)
b_cl.data$barcode=gsub("-1","_b",b_cl.data$barcode)
f_cl.data$barcode=gsub("-1","_f",f_cl.data$barcode)

contig_list <- list(c1_cl.data, c2_cl.data, d1_cl.data, d2_cl.data, a_cl.data, b_cl.data, f_cl.data)
head(contig_list[[1]])

######### Exercise to train your dplyr/ggplot skills 
# check how many TRA/TRB chains are detected in cells
#For each cell barcode, check how many TRA and TRB chains are available in the first sample
contig_list[[2]] %>%
  filter(cdr3!="None",raw_consensus_id!="None")%>%
  group_by(barcode)%>%
  summarise(number_of_IGHs = sum(chain=="IGH"),
            number_of_IGLs = sum(chain=="IGL"),
            number_of_IGKs = sum(chain=="IGK"))

#Check how many cells have each combination of chain number (0 heavy - 1 beta, 1 alpha - 1 beta etc.).
contig_list[[2]] %>%
  filter(cdr3!="None",raw_consensus_id!="None")%>%
  group_by(barcode)%>%
  summarise(number_of_IGHs = sum(chain=="IGH"), number_of_IGLs = sum(chain=="IGL"),number_of_IGKs = sum(chain=="IGK"))%>%
  group_by(number_of_IGHs, number_of_IGLs, number_of_IGKs)%>%
  summarise(number_of_cells = n())


#Repeat the same for all samples
lapply(contig_list, function(one_rep){
  one_rep%>%
    filter(cdr3!="None",raw_consensus_id!="None")%>%
    group_by(barcode)%>%
    summarise(number_of_IGHs = sum(chain=="IGH"), number_of_IGLs = sum(chain=="IGL"),number_of_IGKs = sum(chain=="IGK"))%>%
    group_by(number_of_IGHs, number_of_IGLs, number_of_IGKs)%>%
    summarise(number_of_cells = n())
})

#Iterating over sample numbers, instead of samples - to add this ID to the final data frame
lapply(1:length(contig_list), function(one_rep){
  contig_list[[one_rep]]%>%
    filter(cdr3!="None",raw_consensus_id!="None")%>%
    group_by(barcode)%>%
    summarise(number_of_IGHs = sum(chain=="IGH"), number_of_IGLs = sum(chain=="IGL"),number_of_IGKs = sum(chain=="IGK"))%>%
    group_by(number_of_IGHs, number_of_IGLs, number_of_IGKs)%>%
    summarise(number_of_cells = n())%>%
    mutate(ID=paste0("sample", one_rep))
})

#Joining it all in one data frame
chains_per_samples <- lapply(1:length(contig_list), function(one_rep){
  contig_list[[one_rep]]%>%
    filter(cdr3!="None",raw_consensus_id!="None")%>%
    group_by(barcode)%>%
    summarise(number_of_IGHs = sum(chain=="IGH"), number_of_IGLs = sum(chain=="IGL"),number_of_IGKs = sum(chain=="IGK"))%>%
    group_by(number_of_IGHs, number_of_IGLs, number_of_IGKs)%>%
    summarise(number_of_cells = n())%>%
    mutate(ID = paste0("sample", one_rep))}) %>%
  bind_rows()%>%
  ungroup()

#Plot for each sample number of cells with each combination of chain numbers
chains_per_samples %>%
  ggplot(aes(x=interaction(number_of_IGHs,number_of_IGLs, number_of_IGKs), y= number_of_cells, col=ID)) +
  geom_jitter(height=0, width=0.2, size=3, alpha=0.5)  +   #from here just plot polishing
  theme_bw()+
  xlab("IGH.IGL.IGK chains")+
  scale_color_discrete(name="Sample")

#Plot for each sample fraction of cells with each combination of chain numbers
chains_per_samples %>%
  group_by(ID)%>%
  mutate(fraction_of_cells = number_of_cells/sum(number_of_cells))%>%
  ggplot(aes(x=interaction(number_of_IGHs,number_of_IGLs, number_of_IGKs), y= fraction_of_cells, col=as.factor(ID))) +
  geom_jitter(height=0, width=0.2, size=3, alpha=0.5)+
  theme_bw()+
  xlab("IGH.IGL.IGK chains")+
  scale_color_discrete(name="Sample")


##Filter NAs
contig_list_filtered <-lapply(contig_list, function(one_rep){
  one_rep%>%
    filter(cdr3!="None",raw_consensus_id!="None")
  })

#Clonotypes with 1x IGH and 1x IGL or IGK
contig_list_filtered_2 <-lapply(contig_list, function(one_rep){
  one_rep%>%
    filter(cdr3!="None",raw_consensus_id!="None")%>%
    group_by(barcode)%>%
    mutate(number_of_IGHs = sum(chain=="IGH"), number_of_IGLs = sum(chain=="IGL"),number_of_IGKs = sum(chain=="IGK"))%>%
    filter(number_of_IGHs == 1 & number_of_IGLs == 1 & number_of_IGKs == 0| number_of_IGHs == 1 & number_of_IGLs == 0 & number_of_IGKs == 1)%>%
    ungroup()
})


contig_list[[2]]%>%
  filter(cdr3!="None",raw_consensus_id!="None")%>%
  group_by(barcode)%>%
  mutate(number_of_IGHs = sum(chain=="IGH"), number_of_IGLs = sum(chain=="IGL"),number_of_IGKs = sum(chain=="IGK"))%>%
  filter(number_of_IGHs == 1 & number_of_IGLs == 1 & number_of_IGKs == 0| number_of_IGHs == 1 & number_of_IGLs == 0 & number_of_IGKs == 1)%>%
  ungroup()

contig_list[[2]]%>%
  filter(cdr3!="None",raw_consensus_id!="None")%>%
  group_by(barcode)%>%
  summarise(number_of_IGHs = sum(chain=="IGH"), number_of_IGLs = sum(chain=="IGL"),number_of_IGKs = sum(chain=="IGK"))%>%
  filter(number_of_IGHs == 1 & number_of_IGLs == 1 & number_of_IGKs == 0| number_of_IGHs == 1 & number_of_IGLs == 0 & number_of_IGKs == 1)%>%
  ungroup()%>%
  summarise(number_of_cells = n())


##Generate combined object
combined <- combineBCR(contig_list, samples = c("c1", "c2", "d1", "d2", "a", "b", "f"),  removeNA = TRUE, removeMulti = TRUE)
combined[[7]]
str(combined)
head(combined)

##Percent/total number of unique clonotypes 
quantContig(combined, cloneCall = "gene+nt", scale = T) #percent of unique clonotypes of total size of the size of clonotyeps
quantContig(combined, cloneCall = "gene+nt", scale = F) #number of uniqe clonotypes

quantContig(combined, cloneCall = "gene+nt", scale = T, chain = "IGH") + ggtitle("IGH")#by IGH - percent
quantContig(combined, cloneCall = "gene+nt", scale = F, chain = "IGH") + ggtitle("IGH")#by IGH - number

quantContig(combined, cloneCall = "gene+nt", scale = T, chain = "IGL") + ggtitle("IGL")#by IGL - percent
quantContig(combined, cloneCall = "gene+nt", scale = F, chain = "IGL") + ggtitle("IGL")#by IGL - number


##Abundance of clonotypes
abundanceContig(combined, cloneCall = "gene+nt", scale = F) #How many clonotypes with 1,2,....clones
Abundance_clonotypes <- abundanceContig(combined, cloneCall = "aa", scale = T, exportTable = T)


##Length of clonotypes
abundance <- lengthContig(combined, cloneCall = "aa", chain = "IGH", scale = T)

##Compare clonotypes
compareClonotypes(combined, numbers = 10, samples = c("a", "c1"), cloneCall = "aa", graph = "alluvial", chain = "both")
compareClonotypes(combined, numbers = 10, samples = c("a", "c2"), cloneCall = "aa", graph = "alluvial", chain = "both")
compareClonotypes(combined, numbers = 10, samples = c("a", "f"), cloneCall = "aa", graph = "alluvial", chain = "both")
compareClonotypes(combined, numbers = 10, samples = c("a", "b"), cloneCall = "aa", graph = "alluvial", chain = "both")
compareClonotypes(combined, numbers = 10, samples = c("c1", "a"), cloneCall = "aa", graph = "alluvial", chain = "both")


##Visualise Gene Usage
vizGenes(combined, gene = "C", chain = "IGH", plot = "bar", order = "variance", scale = TRUE) + ggtitle("Antibdoy ispotype usage")
vizGenes(combined, gene = "C", chain = "IGH", plot = "heatmap", order = "variance", scale = TRUE)

vizGenes(combined, gene = "C", chain = "IGL", plot = "bar", order = "variance", scale = TRUE) + ggtitle("Antibdoy ispotype usage")
vizGenes(combined, gene = "C", chain = "IGL", plot = "heatmap", order = "variance", scale = TRUE)

vizGenes(combined, gene = "V", chain = "IGH", plot = "bar", order = "variance", scale = TRUE) + ggtitle("V gene usage")
vizGenes(combined, gene = "V", chain = "IGH", plot = "heatmap", order = "variance", scale = TRUE)

vizGenes(combined, gene = "V", chain = "IGL", plot = "bar", order = "variance", scale = TRUE) + ggtitle("V gene usage")
vizGenes(combined, gene = "V", chain = "IGL", plot = "heatmap", order = "variance", scale = TRUE)

vizGenes(combined, gene = "D", chain = "IGH", plot = "bar", order = "variance", scale = TRUE) + ggtitle("D gene usage")
vizGenes(combined, gene = "D", chain = "IGH", plot = "heatmap", order = "variance", scale = TRUE)

vizGenes(combined, gene = "J", chain = "IGH", plot = "bar", order = "variance", scale = TRUE) + ggtitle("J gene usage")
vizGenes(combined, gene = "J", chain = "IGH", plot = "heatmap", order = "variance", scale = TRUE)

vizGenes(combined, gene = "J", chain = "IGL", plot = "bar", order = "variance", scale = TRUE) + ggtitle("J gene usage")
vizGenes(combined, gene = "J", chain = "IGL", plot = "heatmap", order = "variance", scale = TRUE)

#Antibody isotypes
grid.arrange(IgM, IgD, IGHM, IGHD, ncol = 2)
FeaturePlot(B_cells, features = "Ighg1", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Ighg2b", reduction = "harmony.wnn.umap", cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Ighg2c", reduction = "harmony.wnn.umap",  cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Ighg3", reduction = "harmony.wnn.umap",  cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Igha", reduction = "harmony.wnn.umap",  cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Ighe", reduction = "harmony.wnn.umap",  cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Jchain", reduction = "harmony.wnn.umap",  cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Igkc", reduction = "harmony.wnn.umap",  cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Iglc1", reduction = "harmony.wnn.umap",  cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Iglc2", reduction = "harmony.wnn.umap",  cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = "Iglc3", reduction = "harmony.wnn.umap",  cols = col_con, pt.size = 1.5)
FeaturePlot(B_cells, features = c("Ighm", "Ighd"), blend = TRUE, reduction = "harmony.wnn.umap")
FeaturePlot(B_cells, features = c("Cd79a", "Cd79b"), blend = TRUE, reduction = "harmony.wnn.umap")

#Clonal overlap
clonalOverlap(combined, cloneCall = "gene+nt", 
              method = "morisita")

#Clonotype proportion
clonalProportion(combined, cloneCall = "gene")
clonalProportion(combined, cloneCall = "nt")

###Clonal Homeostasis
clonalHomeostasis(combined, cloneCall = "gene")
clonalHomeostasis(combined, cloneCall = "nt")

##Clonotype and Seurat
head(B_cells[[]])
DimPlot(B_cells, group.by = "cloneType", reduction = "harmony.wnn.umap")
clonalOverlay(B_cells, reduction = "harmony.wnn.umap",freq.cutpoint = 10, bins = 10, facet = "orig.ident") + guides(color = "none")

clonalNetwork(B_cells, 
              reduction = "harmony.wnn.umap", 
              identity = "cluster",
              filter.clones = 10,
              filter.identity = "Cluster 4",
              cloneCall = "aa") + B_cells_p6

###Highlight certain clonotypes
experiment1 <- highlightClonotypes(experiment, cloneCall= "aa", sequence = c("IGHM"))
DimPlot(experiment1, group.by = "highlight")

?highlightClonotypes()
experiment1[[]]

occupiedscRepertoire(experiment, x.axis = "cluster")

####Clonal Network analysis####

clonalNetwork(experiment1, 
              reduction = "wnn.umap", 
              identity = "cluster",
              filter.identity = 15,
              cloneCall = "gene")

clonalNetwork(experiment1, 
              reduction = "wnn.umap", 
              identity = "cluster",
              filter.clones = 10,
              filter.identity = NULL,
              cloneCall = "aa")


