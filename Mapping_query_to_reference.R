## ##Mapping query dataset to the reference## ##
library(SeuratDisk)
library(anndata)
library(SummarizedExperiment)
library(TabulaMurisData)
library(patchwork)
library(scater)
library(ExperimentHub)
library(SingleCellExperiment)
library(iSEE)
library(scran)
library(GEOquery)
library(Seurat)
library(RcmdrMisc)

####Azimuth####
###Save Seurat objects in smaller size format
experiment.diet <- DietSeurat(experiment)
head(experiment.diet[[]])
SaveH5Seurat(experiment.diet, "experiment.diet", overwrite = TURE)

####Load Seurat objects####
#example from https://figshare.com/articles/dataset/Expression_of_97_surface_markers_and_462_mRNAs_in_31586_cells_from_15_leukemic_human_bone_marrow/14780127
bone.marrow.data <- readRDS("AMLs_Scano_projected.rds")
head(bone.marrow.data[[]])
Idents(bone.marrow.data) <- bone.marrow.data$ct

reference <- bone.marrow.data

####Transform AnnData files into Seurat objects - inter-compatibiltiy between scnapy (Python) and Seurat (R)
Convert("tabula-muris-senis-droplet-official-raw-obj.h5ad", dest = "h5Seurat", overwrite = TRUE)
reference1 <- LoadH5Seurat("tabula-muris-senis-droplet-official-raw-obj.h5seurat")
head(reference1[[]])

####TabulaMuris dataset####
eh <- ExperimentHub()
query(eh, "TabulaMurisData")

##Droplet data
droplet <- eh[["EH1617"]]
droplet <- TabulaMurisDroplet()

colnames(colData(droplet))

##Droplet data all
droplet <- scran::computeSumFactors(droplet)
droplet <- scater::logNormCounts(droplet)
droplet <- scater::runPCA(droplet)
droplet <- scater::runUMAP(droplet)
droplet <- scater::runTSNE(droplet)

if (require(iSEE)) {
  iSEE(droplet)
}

droplet.seurat <- as.Seurat(droplet, count = "counts", data = "logcounts")
head(droplet.seurat[[]])

SaveH5Seurat(droplet.seurat, "droplet.seurat", overwrite = TURE)
head(droplet.seurat[[]])

reference <- droplet.seurat

##Droplet data spleen
spleen <- subset(x = droplet.seurat, subset = tissue == "Spleen")
head(spleen[[]])
DimPlot(spleen, reduction = "PCA", group.by = "cell_ontology_class") + NoLegend()
Idents(object = spleen) <- "mouse_id"

####Datasets via Gene Expression Omnibus
meta.data <- read.delim(file.choose("GSE139833_series_matrix.txt"))
reference.data <- read.delim(file.choose("GSE139833_n_gc_m_lz_dz_tpm_grch38.txt"))
reference.data.seurat <- CreateSeuratObject(counts = reference.data, meta.data = meta.data)
head(reference.data.seurat[[]])


###Load suerat objects as .rda files
load(file.choose("NicheData10x.rda"), verbose = TRUE)
head(NicheData10x[[]])
colnames(NicheData10x[[]])
DimPlot(NicheData10x)
reference <- NicheData10x
Idents(reference) <- levels(NicheData10x)
head(NicheData10x[[]])

####Transfer anchors function####
#Load query datasets
experiment <- LoadH5Seurat("SeuratProject.h5Seurat")
query <- experiment

#Select reference dataset
head(reference[[]])
head(query[[]])

DefaultAssay(reference) <- "RNA"
DefaultAssay(query) <- "RNA"

#Perform same normalisation step as for query
reference = SCTransform(reference, verbose = TRUE)
reference[["SCT"]]

#Or

reference <- NormalizeData(reference)
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)

query <- NormalizeData(query)
query <- FindVariableFeatures(query)
query <- ScaleData(query)

DefaultAssay(reference) <- "SCT"
DefaultAssay(query) <- "SCT"

#Find anchors
anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:20, approx = FALSE)

anchors <- FindTransferAnchors(reference = reference, query = query, normalization.method = "SCT", dims = 1:405, approx = FALSE)

#Transfer labels
predictions <- TransferData(anchorset = anchors,refdata = reference$ct)

query <- AddMetaData(object = query, metadata = predictions)
DimPlot(query, label = TRUE, split.by = "orig.ident", repel = TRUE, reduction = "wnn.umap", raster = T) +  ggtitle("Query")
predictions$prediction.score.Classical.Monocytes
FeaturePlot(query,"prediction.score.Classical.Monocytes", reduction = "wnn.umap")
VlnPlot(query,"percent.mt", pt.size = 0 )+NoLegend()
head(query[[]])

query@meta.data[,'Annot'] <-Idents(query)
predictions$prediction.score.max
Data<-query@meta.data
Data$orig.ident
library(dplyr)
Forest <-Data %>% 
  group_by(`Annot`) %>% 
  summarize_all(mean)
Forest<-as.data.frame(Forest)
rownames(Forest)<-Forest$Annot
Forest<-Forest[,c(28:50)]
final<-t(Forest)

pheatmap::pheatmap(final, cluster_rows = T, cluster_cols = T,show_rownames = T, show_colnames = T,
                   cellwidth = 10,cellheight = 10, angle_col = 45)

number_perCluster<- table(query$Annot,query$predicted.id)

percents<-rowPercents(number_perCluster)
dim(percents)
percents<-percents[,c(-25,-26)]
my.colors <- viridisLite::viridis(n = 100, option = "D")
colfunc <- colorRampPalette(c("white" ,"navy"))
colors3<-colfunc(100)

pheatmap::pheatmap(t(percents), cluster_rows = T, cluster_cols = T,show_rownames = T, show_colnames = T,
                   cellwidth = 10,cellheight = 10, angle_col = 45, color = colors3)

number_perCluster<- table(query@meta.data$predicted.id,query@meta.data$orig.ident)

library(reshape2)
DF1 <- melt(number_perCluster, id.var="Cluster")
ggplot(DF1, aes(Var2,value , fill = Var1)) +
  geom_bar(stat= "identity", position = "fill") +
  xlab("Cluster") + ylab("Proportion of cluster")


####Assess query####
cluster.vs.id <- query[[c("predicted.id", "seurat_clusters")]]
DotPlot(query, features = ) + RotatedAxis()
?DotPlot

####Subset individual clusters####
Cluster.3 <- subset(x = experiment, subset = seurat_clusters == "3")

Cluster.4 <- subset(x = query, subset = seurat_clusters == "4")

query.c <- Cluster.4

##Cluster data processing
query.c = SCTransform(query.c, verbose = TRUE)
query.c[["SCT"]]
query.c <- RunPCA(query.c, verbose = FALSE, features = VariableFeatures(object = query.c))

pca_variance <- query.c@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #50


query.c <- FindNeighbors(query.c, dims = 1:50)
query.c <- FindClusters(query.c, resolution = 1.0, verbose = FALSE) #1.0 for the resolution
clustree(query.c, prefix = "SCT_snn_res.") + theme(legend.position="bottom")
query.c <- RunUMAP(query.c, dims = 1:50)
DimPlot(query.c, label = TRUE, reduction = "umap") +  ggtitle("RNA Clustering")

DefaultAssay(query.c) <- "ADT"
VariableFeatures(query.c) <- rownames(query.c[["ADT"]])
query.c <- NormalizeData(query.c, normalization.method = "CLR", margin = 2)
query.c <- ScaleData(query.c)
query.c <- RunPCA(query.c,reduction.name = 'apca')

apca_variance <- query.c@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #26

query.c <- FindMultiModalNeighbors(
  query.c, reduction.list = list("pca", "apca"), 
  dims.list = list(1:50, 1:26), modality.weight.name = "RNA.weight")


query.c <- RunUMAP(query.c, reduction = 'pca', dims = 1:50, assay = 'RNA', 
                      reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
query.c <- RunUMAP(query.c, reduction = 'apca', dims = 1:26, assay = 'ADT', 
                     reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
query.c <- RunUMAP(query.c, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
query.c <- FindClusters(query.c, graph.name = "wsnn", algorithm = 3, resolution = 1.0, verbose = TRUE)

clustree(query.c, prefix = "wsnn_res.") + theme(legend.position="bottom")

head(query.c[[]])


##Visualisation of clusters
DefaultAssay(query.c) <- "RNA"

DimPlot(query.c, label = TRUE, reduction = "wnn.umap", group.by = "predicted.id") +  ggtitle("wnn.umap")
FeaturePlot(query.c, features = "IL10", reduction = "wnn.umap")

Clus4.0 <- FindMarkers(query.c, ident.1 = 0, assay = "ADT")
Clus4.0 <- FindMarkers(query.c, ident.1 = 0, assay = "RNA")

Clus4.1 <- FindMarkers(query.c, ident.1 = 1, assay = "ADT")
Clus4.1 <- FindMarkers(query.c, ident.1 = 1, assay = "RNA")

Clus4.5 <- FindMarkers(query.c, ident.1 = 5, assay = "ADT")
Clus4.5 <- FindMarkers(query.c, ident.1 = 5, assay = "RNA")

head(Cluster.4[[]])
