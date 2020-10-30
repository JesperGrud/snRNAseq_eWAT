### DISCLAIMER
Some of the algorithms are non-deterministic making the results slightly different from run to run.
Many of the algorithms are frequently updated with new releases making the results slightly different from version to version.
Depending on your system, this code below may not produce an exact 1:1 copy of the results. 
This adversely affects clustering, filtering etc. if the code is naively copy-pasted.
For better reproducibility, please download quality-filtered Rds objects and work on downstream analyses. 
<br>

```R
## Load libraries
library(SingleCellExperiment) 
library(simpleSingleCell) 
library(BiocSingular) 
library(umap)
library(Seurat) 
library(Matrix) 
library(scran) 
library(scater) 
library(DropletUtils) 
library(batchelor)
library(harmony)
library(ComplexHeatmap)
library(circlize)
library(MAST)
library(limma)
library(RANN)
library(biomaRt)
library(kBET)
library(LISI)
library(org.Mm.eg.db)
library(dplyr)
library(clusterProfiler)
library(rWikiPathways)
library(GEOquery)
library(edgeR)
library(glmnet)
library(velocyto.R)
library(phateR)
library(ElPiGraph.R)
library(slingshot)
library(TSCAN)
library(monocle3)
library(tradeSeq)
library(seriation)

## Import data
eWAT <- readRDS("eWAT_Annotated.Rds")
eWAT_R3 <- readRDS("eWAT_R3.Rds")
Amb <- readRDS("Ambient.Rds")
Markers <- readRDS("eWAT_Overall_Markers.Rds") # Generated during analysis of the main datasets

## Set parameters for the umap package to be equivalent to those used in Seurat
custom_settings <- umap.defaults
custom_settings$n_epochs <- 200
custom_settings$n_neighbors <- 30
custom_settings$metric <- "correlation"
custom_settings$min_dist <- 0.3
custom_settings$transform_state <- 42
custom_settings$random_state <- 42

## Transfer labels from the main dataset to the new dataset by using the PCA and UMAP from the main datasets
# Transfer variable features and scale the data
VariableFeatures(eWAT_R3) <- VariableFeatures(eWAT) 
eWAT_R3 <- ScaleData(eWAT_R3)

# Create a UMAP embedding of the main dataset. Identical to the one generated using Seurat, but the umap package saves the model to allow for future predictions. 
set.seed(42)
UModel <- umap::umap(d = eWAT@reductions$harmony@cell.embeddings[, 1:15], config = custom_settings, method = "umap-learn") 

# Re-use the loadings of features in the main dataset to perform an identical PC rotation in R3
Loadings <- eWAT@reductions$pca@feature.loadings
Loadings <- Loadings[order(rownames(Loadings)), ]
Counts <- eWAT_R3@assays$RNA@scale.data
Counts <- Counts[order(rownames(Counts)), ]
PCA <- t(Counts) %*% Loadings
colnames(PCA) <- paste0("PCA_", 1:50)
eWAT_R3[["pca"]] <- CreateDimReducObject(embeddings = PCA, key = "PCA_", assay = DefaultAssay(eWAT_R3))  

# Harmonize the R3 PCs across diets. NOTE: This step is non-deterministic! Results vary from run-to-run!
eWAT_R3 <- RunHarmony(object = eWAT_R3, group.by.vars = "Diet", dims.use = 1:20)

# Predict UMAP coordinates for R3 using the model from the main dataset
tmp2 <- predict(UModel, Embeddings(eWAT_R3, "harmony")[, c(1:15)])
colnames(tmp2) <- paste0("UMAP_", 1:2)
eWAT_R3[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = DefaultAssay(eWAT_R3))

# For each nuclei in R3, use the predicted UMAP coordinates to find the closest nuclei in the main data in UMAP space
IDs <- data.frame(ID = nn2(data = UModel$layout, query = tmp2, k = 1)$nn.idx)
IDs$cluster <- 10
for (i in 1:nrow(IDs)) {
  IDs[i, 2] <- as.character(eWAT@meta.data[IDs[i, 1], "seurat_clusters"])
}

# Transfer the labels to the R3 object
eWAT_R3$seurat_clusters <- IDs$cluster

# Set additional annotation
Labels <- eWAT@meta.data[ ,c("Annotation","seurat_clusters","Label")]
Labels <- Labels[ duplicated(Labels$seurat_clusters)==F,]
for (i in 1:nrow(Labels)) {
  eWAT_R3@meta.data[ eWAT_R3@meta.data$seurat_clusters %in% Labels[i,"seurat_clusters"],"Annotation"] <- Labels[i,"Annotation"]
  eWAT_R3@meta.data[ eWAT_R3@meta.data$seurat_clusters %in% Labels[i,"seurat_clusters"],"Label"] <- Labels[i,"Label"]
}
eWAT_R3$Subtype <- eWAT_R3$Annotation
eWAT_R3 <- SetIdent(eWAT_R3, value = eWAT_R3$Annotation)

# Save the object
saveRDS(eWAT_R3, "eWAT_R3.Rds") ## This object is downloadable from Open Science Framework, with futher annotations as subsequent scripts add to the file

## DimPlots
DimPlot(eWAT_R3, label=T)
DimPlot(eWAT_R3, label=T, split.by="Diet")

## Cell type proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(eWAT_R3$Annotation)), nrow=2))
colnames(Composition) <- unique(eWAT_R3$Annotation)
for (i in 1:length(unique(eWAT_R3$Annotation))) {
  Composition[1,i] <- nrow(eWAT_R3@meta.data[ eWAT_R3@meta.data$Dataset == "LFD_R3" & eWAT_R3@meta.data$Annotation == colnames(Composition)[i],])
  Composition[2,i] <- nrow(eWAT_R3@meta.data[ eWAT_R3@meta.data$Dataset == "HFD_R3" & eWAT_R3@meta.data$Annotation == colnames(Composition)[i],])
}
Composition <- Composition/rowSums(Composition)
rownames(Composition) <- c("LFD","HFD")

# Plot it
barplot(as.matrix(Composition), beside=T, las=1, ylim = c(0, 0.5), xlab = "Cell types", ylab = "Fraction of nuclei")

## Gene modules
# Merge the datasets
eWAT <- merge(eWAT, eWAT_R3)

# Extract marker gene results
Clust1 <- Markers[[1]] 
Clust2 <- Markers[[2]] 
Clust3 <- Markers[[3]] 
Clust4 <- Markers[[4]] 
Clust5 <- Markers[[5]] 
Clust6 <- Markers[[6]] 
Clust7 <- Markers[[7]]

# Subset to genes that are enriched or exclusive (!= NS)
Clust1 <- Clust1[Clust1$Marker != "NS", ]
Clust2 <- Clust2[Clust2$Marker != "NS", ]
Clust3 <- Clust3[Clust3$Marker != "NS", ]
Clust4 <- Clust4[Clust4$Marker != "NS", ]
Clust5 <- Clust5[Clust5$Marker != "NS", ]
Clust6 <- Clust6[Clust6$Marker != "NS", ]
Clust7 <- Clust7[Clust7$Marker != "NS", ]

# Extracting the minimum logFC across replicates in pairwise cluster comparisons
# RATIONALE: The genes with the highest minimum fold change are the most specific ones for any given cluster
Clust1$logFC_OvO <- apply(Clust1[, grep("logFC_Cluster", colnames(Clust1))], 1, FUN = "min")
Clust2$logFC_OvO <- apply(Clust2[, grep("logFC_Cluster", colnames(Clust2))], 1, FUN = "min")
Clust3$logFC_OvO <- apply(Clust3[, grep("logFC_Cluster", colnames(Clust3))], 1, FUN = "min")
Clust4$logFC_OvO <- apply(Clust4[, grep("logFC_Cluster", colnames(Clust4))], 1, FUN = "min")
Clust5$logFC_OvO <- apply(Clust5[, grep("logFC_Cluster", colnames(Clust5))], 1, FUN = "min")
Clust6$logFC_OvO <- apply(Clust6[, grep("logFC_Cluster", colnames(Clust6))], 1, FUN = "min")
Clust7$logFC_OvO <- apply(Clust7[, grep("logFC_Cluster", colnames(Clust7))], 1, FUN = "min")

# Computing gene module scores using the top 50 most specific marker genes. NOTE: THIS IS NON-DETERMINISTIC. THE RESULTS VARY SLIGHTLY FROM RUN-TO-RUN.
eWAT <- AddModuleScore(eWAT, features = list(Cluster1 = Clust1[order(-Clust1$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster2 = Clust2[order(-Clust2$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster3 = Clust3[order(-Clust3$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster4 = Clust4[order(-Clust4$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster5 = Clust5[order(-Clust5$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster6 = Clust6[order(-Clust6$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster7 = Clust7[order(-Clust7$logFC_OvO), ][1:50, ]$Symbol))

# Scale the module scores (across all nuclei)
eWAT$Cluster1 <- scale(eWAT$Cluster1)
eWAT$Cluster2 <- scale(eWAT$Cluster2)
eWAT$Cluster3 <- scale(eWAT$Cluster3)
eWAT$Cluster4 <- scale(eWAT$Cluster4)
eWAT$Cluster5 <- scale(eWAT$Cluster5)
eWAT$Cluster6 <- scale(eWAT$Cluster6)
eWAT$Cluster7 <- scale(eWAT$Cluster7)

# Computing mean scaled module scores in each cluster in each replicate
Averages <- as.data.frame(matrix(ncol=length(Markers), nrow=length(Markers)*3))

# Set column and row names
Labels <- eWAT@meta.data[, c("Annotation","Label")]
Labels <- Labels[ duplicated(Labels$Label)==F,]
Labels <- Labels[ order(Labels$Label),]
colnames(Averages) <- paste(Labels$Annotation,"module",sep="_")
rownames(Averages)[seq(1,(nrow(Averages)-2),by=3)] <- paste(Labels$Annotation,"R1",sep="_")
rownames(Averages)[seq(2,nrow(Averages)-1,by=3)] <- paste(Labels$Annotation,"R2",sep="_")
rownames(Averages)[seq(3,nrow(Averages),by=3)] <- paste(Labels$Annotation,"R3",sep="_")

# Calculate averages
for (module in 1:length(Markers)) {
  Counter <- 1
  for (label in 1:length(Markers)) {
    for (rep in c("R1","R2","R3")) {
      Averages[Counter,module] <- mean(eWAT@meta.data[ eWAT@meta.data$Label == label & eWAT$Replicate == rep,paste("Cluster",module,sep="")])
      Counter <- Counter + 1
    }
  }
}

## Plot mean scaled module scores
col <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(as.matrix(Averages), cluster_columns = F, cluster_rows = F, col = col)

## FeaturePlots
FeaturePlot(eWAT_R3, "Adgre1", max.cutoff = 3, pt.size = 0.001)
FeaturePlot(eWAT_R3, "Pdgfra", max.cutoff = 3, pt.size = 0.001)
FeaturePlot(eWAT_R3, "Lipe", max.cutoff = 3, pt.size = 0.001)
FeaturePlot(eWAT_R3, "Gpm6a", max.cutoff = 3, pt.size = 0.001)
```
[Back to start](../README.md)<br>



