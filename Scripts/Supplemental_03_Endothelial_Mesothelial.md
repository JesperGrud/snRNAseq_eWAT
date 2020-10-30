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

## Import data and subset
Endomeso <- readRDS("eWAT_Endothelial_Mesothelial.Rds")
eWAT_R3 <- readRDS("eWAT_R3.Rds")
Endomeso_R3 <- subset(eWAT_R3, subset = Annotation %in% c("Mesothelial", "Endothelial"))
Amb <- readRDS("Ambient.Rds")
Markers <- readRDS("eWAT_Endothelial_Mesothelial_Markers.Rds") # Generated during analysis of the main datasets
# Endomeso_R3 <- readRDS("eWAT_R3_Endothelial_Mesothelial.Rds")

## Set parameters for the umap package to be equivalent to those used in Seurat
custom_settings <- umap.defaults
custom_settings$n_epochs <- 500
custom_settings$n_neighbors <- 30
custom_settings$metric <- "correlation"
custom_settings$min_dist <- 0.3
custom_settings$transform_state <- 42
custom_settings$random_state <- 42

## Transfer labels from the main dataset to the new dataset by using the PCA and UMAP from the main datasets
# Transfer variable features and scale the data
VariableFeatures(Endomeso_R3) <- VariableFeatures(Endomeso) 
Endomeso_R3 <- ScaleData(Endomeso_R3)

# Create a UMAP embedding of the main dataset. Identical to the one generated using Seurat, but the umap package saves the model to allow for future predictions. 
set.seed(42)
UModel <- umap::umap(d = Endomeso@reductions$harmony@cell.embeddings[, 1:12], config = custom_settings, method = "umap-learn") 

# Re-use the loadings of features in the main dataset to perform an identical PC rotation in R3
Loadings <- Endomeso@reductions$pca@feature.loadings
Loadings <- Loadings[order(rownames(Loadings)), ]
Counts <- Endomeso_R3@assays$RNA@scale.data
Counts <- Counts[order(rownames(Counts)), ]
PCA <- t(Counts) %*% Loadings
colnames(PCA) <- paste0("PCA_", 1:50)
Endomeso_R3[["pca"]] <- CreateDimReducObject(embeddings = PCA, key = "PCA_", assay = DefaultAssay(Endomeso_R3))  

# Harmonize the R3 PCs across diets. NOTE: This step is non-deterministic! Results vary from run-to-run!
Endomeso_R3 <- RunHarmony(object = Endomeso_R3, group.by.vars = "Diet")

# Predict UMAP coordinates for R3 using the model from the main dataset
tmp2 <- predict(UModel, Embeddings(Endomeso_R3, "harmony")[, c(1:12)])
colnames(tmp2) <- paste0("UMAP_", 1:2)
Endomeso_R3[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = DefaultAssay(Endomeso_R3))

# For each nuclei in R3, use the predicted UMAP coordinates to find the closest nuclei in the main data in UMAP space
IDs <- data.frame(ID = nn2(data = UModel$layout, query = tmp2, k = 1)$nn.idx)
IDs$cluster <- 10
for (i in 1:nrow(IDs)) {
  IDs[i, 2] <- as.character(Endomeso@meta.data[IDs[i, 1], "seurat_clusters"])
}

# Transfer the labels to the R3 object
Endomeso_R3$seurat_clusters <- IDs$cluster

# Set additional annotation
Labels <- Endomeso@meta.data[ ,c("Subtype","seurat_clusters","Label")]
Labels <- Labels[ duplicated(Labels$seurat_clusters)==F,]
for (i in 1:nrow(Labels)) {
  Endomeso_R3@meta.data[ Endomeso_R3@meta.data$seurat_clusters %in% Labels[i,"seurat_clusters"],"Subtype"] <- Labels[i,"Subtype"]
  Endomeso_R3@meta.data[ Endomeso_R3@meta.data$seurat_clusters %in% Labels[i,"seurat_clusters"],"Label"] <- Labels[i,"Label"]
}

# Transfer all subtype labels to overall object
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Subtype == "EPC",]),"Subtype"] <- "EPC"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Subtype == "IMC",]),"Subtype"] <- "IMC"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Subtype == "LEC",]),"Subtype"] <- "LEC"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Subtype == "MC",]),"Subtype"] <- "MC"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Subtype == "VEC",]),"Subtype"] <- "VEC"

# Save the objects
saveRDS(Endomeso_R3, "eWAT_R3_Endothelial_Mesothelial.Rds") ## This object is downloadable from Open Science Framework
saveRDS(eWAT_R3, "eWAT_R3.Rds") ## This object is downloadable from Open Science Framework, with futher annotations as subsequent scripts add to the file

## DimPlots
DimPlot(Endomeso_R3, label=T, group.by="Subtype")
DimPlot(Endomeso_R3, label=T, group.by="Subtype", split.by="Diet")

## Cell type proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(Endomeso_R3$Subtype)), nrow=2))
colnames(Composition) <- unique(Endomeso_R3$Subtype)
for (i in 1:length(unique(Endomeso_R3$Subtype))) {
  Composition[1,i] <- nrow(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Dataset == "LFD_R3" & Endomeso_R3@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Dataset == "LFD_R3",])
  Composition[2,i] <- nrow(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Dataset == "HFD_R3" & Endomeso_R3@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Dataset == "HFD_R3",])
}
rownames(Composition) <- c("LFD","HFD")

# Plot it
barplot(as.matrix(Composition), beside=T, las=1, ylim = c(0, 1), xlab = "Subtypes", ylab = "Fraction of nuclei")

## Gene modules
# Extract marker gene results
Clust1 <- Markers[[1]] 
Clust2 <- Markers[[2]] 
Clust3 <- Markers[[3]] 
Clust4 <- Markers[[4]] 
Clust5 <- Markers[[5]] 

# Subset to genes that are enriched or exclusive (!= NS)
Clust1 <- Clust1[Clust1$Marker != "NS", ]
Clust2 <- Clust2[Clust2$Marker != "NS", ]
Clust3 <- Clust3[Clust3$Marker != "NS", ]
Clust4 <- Clust4[Clust4$Marker != "NS", ]
Clust5 <- Clust5[Clust5$Marker != "NS", ]

# Extracting the minimum logFC across replicates in pairwise cluster comparisons
# RATIONALE: The genes with the highest minimum fold change are the most specific ones for any given cluster
Clust1$logFC_OvO <- apply(Clust1[, grep("logFC_Cluster", colnames(Clust1))], 1, FUN = "min")
Clust2$logFC_OvO <- apply(Clust2[, grep("logFC_Cluster", colnames(Clust2))], 1, FUN = "min")
Clust3$logFC_OvO <- apply(Clust3[, grep("logFC_Cluster", colnames(Clust3))], 1, FUN = "min")
Clust4$logFC_OvO <- apply(Clust4[, grep("logFC_Cluster", colnames(Clust4))], 1, FUN = "min")
Clust5$logFC_OvO <- apply(Clust5[, grep("logFC_Cluster", colnames(Clust5))], 1, FUN = "min")

# Computing gene module scores using the top 50 most specific marker genes. NOTE: THIS IS NON-DETERMINISTIC. THE RESULTS VARY SLIGHTLY FROM RUN-TO-RUN.
Endomeso_R3 <- AddModuleScore(Endomeso_R3, features = list(Cluster1 = Clust1[order(-Clust1$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster2 = Clust2[order(-Clust2$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster3 = Clust3[order(-Clust3$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster4 = Clust4[order(-Clust4$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster5 = Clust5[order(-Clust5$logFC_OvO), ][1:50, ]$Symbol))

# Scale the module scores (across all nuclei)
Endomeso_R3$Cluster1 <- scale(Endomeso_R3$Cluster1)
Endomeso_R3$Cluster2 <- scale(Endomeso_R3$Cluster2)
Endomeso_R3$Cluster3 <- scale(Endomeso_R3$Cluster3)
Endomeso_R3$Cluster4 <- scale(Endomeso_R3$Cluster4)
Endomeso_R3$Cluster5 <- scale(Endomeso_R3$Cluster5)

# Computing mean scaled module scores in each cluster in each replicate
Averages <- as.data.frame(matrix(ncol=length(Markers), nrow=length(Markers)))

# Set column and row names
Labels <- Endomeso_R3@meta.data[, c("Subtype","Label")]
Labels <- Labels[ duplicated(Labels$Label)==F,]
Labels <- Labels[ order(Labels$Label),]
colnames(Averages) <- paste(Labels$Subtype,"module",sep="_")
rownames(Averages) <- Labels$Subtype

# Calculate averages
for (module in 1:length(Markers)) {
  for (label in 1:length(Markers)) {
      Averages[label,module] <- mean(Endomeso_R3@meta.data[ Endomeso_R3@meta.data$Label == label,paste("Cluster",module,sep="")])
  }
}

## Plot mean scaled module scores
col <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
Heatmap(as.matrix(Averages), cluster_columns = F, cluster_rows = F, col = col)
```
[Back to start](../README.md)<br>
