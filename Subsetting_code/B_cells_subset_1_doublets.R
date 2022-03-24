####Doublet finder for B_cells####
library(DoubletFinder)

newBcells <- B_cells
head(B_cells[[]])


mark_doublets = function(seurat_object){
  
  message("Performing filter by number of genes and mitochondrial percentage.")
  seurat_object = subset(seurat_object, subset = nFeature_RNA > 200  & percent.Mt < 5)
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  
  # SCTransform
  # change the current plan to access parallelization
  plan("multiprocess", workers = 4)
  plan()
  # Give more memory to globals to prevent crashing by large Seurat objects
  options(future.globals.maxSize= 2097152000) # 2Gb
  
  seurat_object = SCTransform(seurat_object, verbose = FALSE)
  
  seurat_object = RunPCA(seurat_object, verbose = FALSE)
  
  seurat_object = FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object = FindClusters(seurat_object, verbose = FALSE)
  
  seurat_object = RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  
  sweep.res.list = paramSweep_v3(seurat_object, PCs = 1:10, sct = TRUE)
  sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn = find.pK(sweep.stats)
  
  # This plot will help you find the correct value for pK
  ggplot(bcmvn,aes(x=pK,y=BCmetric, group = 1)) + geom_line() + geom_point()
  # select pK that maximizes the BCmetric - if there are two peaks with the same height, choose the one with lowest pK
  bcmvn$pK = as.numeric(as.character(bcmvn$pK))
  
  pk = bcmvn[which(bcmvn$BCmetric == max(bcmvn$BCmetric)),"pK"]
  
  message("Your best pK value is ", pk)
  
  # Homotypic doublet proportion estimation - from number of clusters
  annotations =  as.character(seurat_object$seurat_clusters)
  homotypic.prop <- modelHomotypic(annotations)   
  
  # Estimate number of total expected doublets
  # Assuming 2.07% doublet formation rate (10X rate for 2700 recovered cells) - tailor for your dataset
  expected_doublet_rate = 0.027
  nExp_poi = round(expected_doublet_rate*nrow(seurat_object@meta.data)) 
  # Adjusting number based on homotypic doublet proportion expected
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # All confidence doublets - not adjusted for homotypic doublets
  seurat_object.doublet <-
    doubletFinder_v3(
      seurat_object,
      PCs = 1:10,
      pN = 0.25,
      pK = pk,
      nExp = nExp_poi, 
      reuse.pANN = FALSE,
      sct = TRUE
    )
  
  # High confidence doublets - adjusted for homotypic doublets
  
  seurat_object.doublet <-
    doubletFinder_v3(
      seurat_object.doublet,
      PCs = 1:10,
      pN = 0.25,
      pK = pk,
      nExp = nExp_poi.adj, 
      reuse.pANN = colnames(seurat_object.doublet@meta.data)[ncol(seurat_object.doublet@meta.data)-1],
      sct = TRUE
    )
  
  # Integrate info all and high confidence
  seurat_object.doublet@meta.data$DF_confidence = "Singlet"
  seurat_object.doublet@meta.data[seurat_object.doublet[[colnames(seurat_object.doublet@meta.data)[ncol(seurat_object.doublet@meta.data) -
                                                                                                     2]]] == "Doublet" &
                                    seurat_object.doublet[[colnames(seurat_object.doublet@meta.data)[ncol(seurat_object.doublet@meta.data)-1]]] == "Doublet",
                                  "DF_confidence"] = "High_confidence_doublet"
  
  seurat_object.doublet@meta.data[seurat_object.doublet[[colnames(seurat_object.doublet@meta.data)[ncol(seurat_object.doublet@meta.data) -
                                                                                                     2]]] == "Doublet" &
                                    seurat_object.doublet[[colnames(seurat_object.doublet@meta.data)[ncol(seurat_object.doublet@meta.data)-1]]] == "Singlet",
                                  "DF_confidence"] = "Low_confidence_doublet"
  
  return(seurat_object.doublet)
  
  
}

Bcells_doublet = mark_doublets(newBcells)
head(Bcells_doublet[[]])

# Saving only doublet information columns
write.csv(Bcells_doublet@meta.data[(ncol(Bcells_doublet@meta.data)-2):ncol(Bcells_doublet@meta.data)], 
          file = "~/Desktop/Bcells_double.csv")

# Read in results from doublet calling
doublets = read.csv("~/Desktop/Bcells_double.csv", row.names = 1)
# Add predicted doublets to metadata
newBcells@meta.data$DF_confidence = "Singlet"
newBcells@meta.data[rownames(doublets[doublets[ncol(doublets)] == "High_confidence_doublet",]),"DF_confidence"] = "High_confidence_doublet"
newBcells@meta.data[rownames(doublets[doublets[ncol(doublets)] == "Low_confidence_doublet",]),"DF_confidence"] = "Low_confidence_doublet"

message("There are ", sum(newBcells@meta.data$DF_confidence!="Singlet"), " predicted doublets in our data")
message("That's ", round((sum(newBcells@meta.data$DF_confidence!="Singlet")/nrow(newBcells@meta.data))*100), "% of our our data")


calculate_SCT_PCA_UMAP_neighbors_clusters = function(seurat_object){
  
  seurat_object = SCTransform(seurat_object, verbose = FALSE)
  seurat_object <- RunPCA(seurat_object, verbose = FALSE)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object <- FindClusters(seurat_object, verbose = FALSE)
  
  seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  
  
  
  return(seurat_object)
}

# Saving in a new object so we don't disturb the old un-transformed, unfiltered one
newBcells.doublet = calculate_SCT_PCA_UMAP_neighbors_clusters(newBcells)

p = DimPlot(
  newBcells.doublet,
  label = F,
  group.by = "DF_confidence",
  cols = c(
    'Singlet' = '#28C2FF',
    'High_confidence_doublet' = '#BD6B73',
    'Low_confidence_doublet' = '#FFA630'
  ),
  pt.size = 1.2,
  reduction = "harmony.wnn.umap",
  order = c('High_confidence_doublet','Low_confidence_doublet' ) # Plot Doublets on top
) +
  ggtitle("DoubletFinder doublets") +
  theme_void() +
  theme(plot.title = element_text(color="black", size=20, face="bold")) +
  scale_colour_manual(labels = c("Singlets", "High confidence doublets (HCDs)", "Low confidence doublets (LCDs)"), values = c("#28C2FF", "#BD6B73", "#FFA630"))


newBcells.doublet.high <- subset(newBcells.doublet, subset = DF_confidence == "High_confidence_doublet")
p_high <- FeaturePlot(newBcells.doublet.high, features = "Cyp11a1", reduction = "harmony.wnn.umap")
p_high <- p_high + theme_void() + ggtitle("Cyp11a1 expression among HDCs") + theme(plot.title = element_text(color="black", size=20, face="bold"))

newBcells.doublet.low <- subset(newBcells.doublet, subset = DF_confidence == "Low_confidence_doublet")
FeaturePlot(newBcells.doublet.low, features = "Cyp11a1", reduction = "harmony.wnn.umap")
p_low <- FeaturePlot(newBcells.doublet.low, features = "Cyp11a1", reduction = "harmony.wnn.umap")
p_low <- p_low + theme_void() + ggtitle("Cyp11a1 expression among LDCs") + theme(plot.title = element_text(color="black", size=20, face="bold"))

B_cells_p6 + p + p_high + p_low
