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
FAP <- readRDS("eWAT_FAP.Rds")
eWAT_R3 <- readRDS("eWAT_R3.Rds")
FAP_R3 <- subset(eWAT_R3, subset = Annotation == "FAP")
Amb <- readRDS("Ambient.Rds")
Markers <- readRDS("eWAT_FAP_Markers.Rds") # Generated during analysis of the main datasets
# FAP_R3 <- readRDS("eWAT_R3_FAP.Rds")

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
VariableFeatures(FAP_R3) <- VariableFeatures(FAP) 
FAP_R3 <- ScaleData(FAP_R3)

# Create a UMAP embedding of the main dataset. Identical to the one generated using Seurat, but the umap package saves the model to allow for future predictions. 
set.seed(42)
UModel <- umap::umap(d = FAP@reductions$harmony@cell.embeddings[, 1:20], config = custom_settings, method = "umap-learn") 

# Re-use the loadings of features in the main dataset to perform an identical PC rotation in R3
Loadings <- FAP@reductions$pca@feature.loadings
Loadings <- Loadings[order(rownames(Loadings)), ]
Counts <- FAP_R3@assays$RNA@scale.data
Counts <- Counts[order(rownames(Counts)), ]
PCA <- t(Counts) %*% Loadings
colnames(PCA) <- paste0("PCA_", 1:50)
FAP_R3[["pca"]] <- CreateDimReducObject(embeddings = PCA, key = "PCA_", assay = DefaultAssay(FAP_R3))  

# Harmonize the R3 PCs across diets. NOTE: This step is non-deterministic! Results vary from run-to-run!
FAP_R3 <- RunHarmony(object = FAP_R3, group.by.vars = "Diet")

# Predict UMAP coordinates for R3 using the model from the main dataset
tmp2 <- predict(UModel, Embeddings(FAP_R3, "harmony")[, c(1:20)])
colnames(tmp2) <- paste0("UMAP_", 1:2)
FAP_R3[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = DefaultAssay(FAP_R3))

# For each nuclei in R3, use the predicted UMAP coordinates to find the closest nuclei in the main data in UMAP space
IDs <- data.frame(ID = nn2(data = UModel$layout, query = tmp2, k = 1)$nn.idx)
IDs$cluster <- 10
for (i in 1:nrow(IDs)) {
  IDs[i, 2] <- as.character(FAP@meta.data[IDs[i, 1], "seurat_clusters"])
}

# Transfer the labels to the R3 object
FAP_R3$seurat_clusters <- IDs$cluster

# Set additional annotation
Labels <- FAP@meta.data[ ,c("Subtype","seurat_clusters","Label")]
Labels <- Labels[ duplicated(Labels$seurat_clusters)==F,]
for (i in 1:nrow(Labels)) {
  FAP_R3@meta.data[ FAP_R3@meta.data$seurat_clusters %in% Labels[i,"seurat_clusters"],"Subtype"] <- Labels[i,"Subtype"]
  FAP_R3@meta.data[ FAP_R3@meta.data$seurat_clusters %in% Labels[i,"seurat_clusters"],"Label"] <- Labels[i,"Label"]
}

# Transfer all subtype labels to overall object
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(FAP_R3@meta.data[ FAP_R3@meta.data$Subtype == "FAP1",]),"Subtype"] <- "FAP1"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(FAP_R3@meta.data[ FAP_R3@meta.data$Subtype == "FAP2",]),"Subtype"] <- "FAP2"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(FAP_R3@meta.data[ FAP_R3@meta.data$Subtype == "FAP3",]),"Subtype"] <- "FAP3"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(FAP_R3@meta.data[ FAP_R3@meta.data$Subtype == "FAP4",]),"Subtype"] <- "FAP4"

# Save the objects
saveRDS(FAP_R3, "eWAT_R3_FAP.Rds") ## This object is downloadable from Open Science Framework
saveRDS(eWAT_R3, "eWAT_R3.Rds") ## This object is downloadable from Open Science Framework, with futher annotations as subsequent scripts add to the file

## DimPlots
DimPlot(FAP_R3, group.by="Subtype", label=T)
DimPlot(FAP_R3, group.by="Subtype", label=T, split.by="Diet")

## Cell type proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(FAP_R3$Subtype)), nrow=2))
colnames(Composition) <- unique(FAP_R3$Subtype)
for (i in 1:length(unique(FAP_R3$Subtype))) {
  Composition[1,i] <- nrow(FAP_R3@meta.data[ FAP_R3@meta.data$Dataset == "LFD_R3" & FAP_R3@meta.data$Subtype == colnames(Composition)[i],]) / nrow(FAP_R3@meta.data[ FAP_R3@meta.data$Dataset == "LFD_R3",])
  Composition[2,i] <- nrow(FAP_R3@meta.data[ FAP_R3@meta.data$Dataset == "HFD_R3" & FAP_R3@meta.data$Subtype == colnames(Composition)[i],]) / nrow(FAP_R3@meta.data[ FAP_R3@meta.data$Dataset == "HFD_R3",])
}
rownames(Composition) <- c("LFD","HFD")

# Plot it
barplot(as.matrix(Composition), beside=T, las=1, ylim = c(0, 0.5), xlab = "Subtypes", ylab = "Fraction of nuclei")

## Gene modules
# Extract marker gene results
Clust1 <- Markers[[1]] 
Clust2 <- Markers[[2]] 
Clust3 <- Markers[[3]] 
Clust4 <- Markers[[4]] 

# Subset to genes that are enriched or exclusive (!= NS)
Clust1 <- Clust1[Clust1$Marker != "NS", ]
Clust2 <- Clust2[Clust2$Marker != "NS", ]
Clust3 <- Clust3[Clust3$Marker != "NS", ]
Clust4 <- Clust4[Clust4$Marker != "NS", ]

# Extracting the minimum logFC across replicates in pairwise cluster comparisons
# RATIONALE: The genes with the highest minimum fold change are the most specific ones for any given cluster
Clust1$logFC_OvO <- apply(Clust1[, grep("logFC_Cluster", colnames(Clust1))], 1, FUN = "min")
Clust2$logFC_OvO <- apply(Clust2[, grep("logFC_Cluster", colnames(Clust2))], 1, FUN = "min")
Clust3$logFC_OvO <- apply(Clust3[, grep("logFC_Cluster", colnames(Clust3))], 1, FUN = "min")
Clust4$logFC_OvO <- apply(Clust4[, grep("logFC_Cluster", colnames(Clust4))], 1, FUN = "min")

# Computing gene module scores using the top 50 most specific marker genes. NOTE: THIS IS NON-DETERMINISTIC. THE RESULTS VARY SLIGHTLY FROM RUN-TO-RUN.
FAP_R3 <- AddModuleScore(FAP_R3, features = list(Cluster1 = Clust1[order(-Clust1$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster2 = Clust2[order(-Clust2$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster3 = Clust3[order(-Clust3$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster4 = Clust4[order(-Clust4$logFC_OvO), ][1:50, ]$Symbol))

# Scale the module scores (across all nuclei)
FAP_R3$Cluster1 <- scale(FAP_R3$Cluster1)
FAP_R3$Cluster2 <- scale(FAP_R3$Cluster2)
FAP_R3$Cluster3 <- scale(FAP_R3$Cluster3)
FAP_R3$Cluster4 <- scale(FAP_R3$Cluster4)

# Computing mean scaled module scores in each cluster in each replicate
Averages <- as.data.frame(matrix(ncol=length(Markers), nrow=length(Markers)))

# Set column and row names
Labels <- FAP_R3@meta.data[, c("Subtype","Label")]
Labels <- Labels[ duplicated(Labels$Label)==F,]
Labels <- Labels[ order(Labels$Label),]
colnames(Averages) <- paste(Labels$Subtype,"module",sep="_")
rownames(Averages) <- Labels$Subtype

# Calculate averages
for (module in 1:length(Markers)) {
  for (label in 1:length(Markers)) {
      Averages[label,module] <- mean(FAP_R3@meta.data[ FAP_R3@meta.data$Label == label,paste("Cluster",module,sep="")])
  }
}

## Plot mean scaled module scores
col <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
Heatmap(as.matrix(Averages), cluster_columns = F, cluster_rows = F, col = col)

## Diet-dependent gene regulation
# Import clustering results
Clusters <- readRDS("eWAT_FAP_Diet_Clusters.Rds") # Generated during analysis of the main datasets

# Extract diet-regulated genes in all subtypes
Cluster0 <- Markers[[1]]
Cluster0 <- Cluster0[Cluster0$Diet != "NS", ]
Cluster1 <- Markers[[2]]
Cluster1 <- Cluster1[Cluster1$Diet != "NS", ]
Cluster2 <- Markers[[3]]
Cluster2 <- Cluster2[Cluster2$Diet != "NS", ]
Cluster3 <- Markers[[4]]
Cluster3 <- Cluster3[Cluster3$Diet != "NS", ]

# Calculate average expression for diet-regulated genes
Averages <- as.data.frame(matrix(ncol = 1, nrow = length(unique(sort(c(Cluster0[, 1], Cluster1[, 1], Cluster2[, 1], Cluster3[, 1]))))))
Averages[, 1] <- unique(sort(c(Cluster0[, 1], Cluster1[, 1], Cluster2[, 1], Cluster3[, 1])))
for (i in 1:nrow(Averages)) {
  Averages[i, 2] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "1" & FAP@meta.data$Diet == "LFD" & FAP@meta.data$Replicate == "R1", ])])))
  Averages[i, 3] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "1" & FAP@meta.data$Diet == "LFD" & FAP@meta.data$Replicate == "R2", ])])))
  Averages[i, 4] <- log1p(mean(expm1(FAP_R3@assays$RNA@data[rownames(FAP_R3@assays$RNA@data) %in% Averages[i, 1], colnames(FAP_R3@assays$RNA@data) %in% rownames(FAP_R3@meta.data[FAP_R3@meta.data$Label == "1" & FAP_R3@meta.data$Diet == "LFD", ])])))
  Averages[i, 5] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "1" & FAP@meta.data$Diet == "HFD" & FAP@meta.data$Replicate == "R1", ])])))
  Averages[i, 6] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "1" & FAP@meta.data$Diet == "HFD" & FAP@meta.data$Replicate == "R2", ])])))
  Averages[i, 7] <- log1p(mean(expm1(FAP_R3@assays$RNA@data[rownames(FAP_R3@assays$RNA@data) %in% Averages[i, 1], colnames(FAP_R3@assays$RNA@data) %in% rownames(FAP_R3@meta.data[FAP_R3@meta.data$Label == "1" & FAP_R3@meta.data$Diet == "HFD", ])])))
  Averages[i, 8] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "2" & FAP@meta.data$Diet == "LFD" & FAP@meta.data$Replicate == "R1", ])])))
  Averages[i, 9] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "2" & FAP@meta.data$Diet == "LFD" & FAP@meta.data$Replicate == "R2", ])])))
  Averages[i, 10] <- log1p(mean(expm1(FAP_R3@assays$RNA@data[rownames(FAP_R3@assays$RNA@data) %in% Averages[i, 1], colnames(FAP_R3@assays$RNA@data) %in% rownames(FAP_R3@meta.data[FAP_R3@meta.data$Label == "2" & FAP_R3@meta.data$Diet == "LFD", ])])))
  Averages[i, 11] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "2" & FAP@meta.data$Diet == "HFD" & FAP@meta.data$Replicate == "R1", ])])))
  Averages[i, 12] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "2" & FAP@meta.data$Diet == "HFD" & FAP@meta.data$Replicate == "R2", ])])))
  Averages[i, 13] <- log1p(mean(expm1(FAP_R3@assays$RNA@data[rownames(FAP_R3@assays$RNA@data) %in% Averages[i, 1], colnames(FAP_R3@assays$RNA@data) %in% rownames(FAP_R3@meta.data[FAP_R3@meta.data$Label == "2" & FAP_R3@meta.data$Diet == "HFD", ])])))
  Averages[i, 14] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "3" & FAP@meta.data$Diet == "LFD" & FAP@meta.data$Replicate == "R1", ])])))
  Averages[i, 15] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "3" & FAP@meta.data$Diet == "LFD" & FAP@meta.data$Replicate == "R2", ])])))
  Averages[i, 16] <- log1p(mean(expm1(FAP_R3@assays$RNA@data[rownames(FAP_R3@assays$RNA@data) %in% Averages[i, 1], colnames(FAP_R3@assays$RNA@data) %in% rownames(FAP_R3@meta.data[FAP_R3@meta.data$Label == "3" & FAP_R3@meta.data$Diet == "LFD", ])])))
  Averages[i, 17] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "3" & FAP@meta.data$Diet == "HFD" & FAP@meta.data$Replicate == "R1", ])])))
  Averages[i, 18] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "3" & FAP@meta.data$Diet == "HFD" & FAP@meta.data$Replicate == "R2", ])])))
  Averages[i, 19] <- log1p(mean(expm1(FAP_R3@assays$RNA@data[rownames(FAP_R3@assays$RNA@data) %in% Averages[i, 1], colnames(FAP_R3@assays$RNA@data) %in% rownames(FAP_R3@meta.data[FAP_R3@meta.data$Label == "3" & FAP_R3@meta.data$Diet == "HFD", ])])))
  Averages[i, 20] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "4" & FAP@meta.data$Diet == "LFD" & FAP@meta.data$Replicate == "R1", ])])))
  Averages[i, 21] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "4" & FAP@meta.data$Diet == "LFD" & FAP@meta.data$Replicate == "R2", ])])))
  Averages[i, 22] <- log1p(mean(expm1(FAP_R3@assays$RNA@data[rownames(FAP_R3@assays$RNA@data) %in% Averages[i, 1], colnames(FAP_R3@assays$RNA@data) %in% rownames(FAP_R3@meta.data[FAP_R3@meta.data$Label == "4" & FAP_R3@meta.data$Diet == "LFD", ])])))
  Averages[i, 23] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "4" & FAP@meta.data$Diet == "HFD" & FAP@meta.data$Replicate == "R1", ])])))
  Averages[i, 24] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[FAP@meta.data$Label == "4" & FAP@meta.data$Diet == "HFD" & FAP@meta.data$Replicate == "R2", ])])))
  Averages[i, 25] <- log1p(mean(expm1(FAP_R3@assays$RNA@data[rownames(FAP_R3@assays$RNA@data) %in% Averages[i, 1], colnames(FAP_R3@assays$RNA@data) %in% rownames(FAP_R3@meta.data[FAP_R3@meta.data$Label == "4" & FAP_R3@meta.data$Diet == "HFD", ])])))
}
rownames(Averages) <- Averages[, 1]

# Scale each experiment
Exp1 <- Averages[, c(2, 5, 8, 11, 14, 17, 20, 23)]
Exp1 <- t(scale(t(Exp1)))
Exp2 <- Averages[, c(3, 6, 9, 12, 15, 18, 21, 24)]
Exp2 <- t(scale(t(Exp2)))
Exp3 <- Averages[, c(4, 7, 10, 13, 16, 19, 22, 25)]
Exp3 <- t(scale(t(Exp3)))
Averages[, c(2, 5, 8, 11, 14, 17, 20, 23)] <- Exp1
Averages[, c(3, 6, 9, 12, 15, 18, 21, 24)] <- Exp2
Averages[, c(4, 7, 10, 13, 16, 19, 22, 25)] <- Exp3

# Set NAs to 0 to avoid errors
for (i in 2:ncol(Averages)) { Averages[ is.na(Averages[,i]),i] <- 0 }

# Plot the heatmap with order and split based on the original clusters
Lens <- unlist(lapply(Clusters, FUN="length"))
SplitVec <- c()
for (i in 1:length(Lens)) { SplitVec <- c(SplitVec, rep(paste("Cluster",i, sep=""), Lens[i])) }
Heatmap(as.matrix(Averages[unlist(Clusters), 2:ncol(Averages)]), show_row_names = F, cluster_columns = F, cluster_rows = F, row_split = SplitVec)
```
[Back to start](../README.md)<br>