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
Immune <- readRDS("eWAT_Immune.Rds")
eWAT_R3 <- readRDS("eWAT_R3.Rds")
Immune_R3 <- subset(eWAT_R3, subset = Annotation == "Immune")
Amb <- readRDS("Ambient.Rds")
Markers <- readRDS("eWAT_Immune_Markers.Rds") # Generated during analysis of the main datasets
# Immune_R3 <- readRDS("eWAT_R3_Immune.Rds")

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
VariableFeatures(Immune_R3) <- VariableFeatures(Immune) 
Immune_R3 <- ScaleData(Immune_R3)

# Create a UMAP embedding of the main dataset. Identical to the one generated using Seurat, but the umap package saves the model to allow for future predictions. 
set.seed(42)
UModel <- umap::umap(d = Immune@reductions$harmony@cell.embeddings[, 1:16], config = custom_settings, method = "umap-learn") 

# Re-use the loadings of features in the main dataset to perform an identical PC rotation in R3
Loadings <- Immune@reductions$pca@feature.loadings
Loadings <- Loadings[order(rownames(Loadings)), ]
Counts <- Immune_R3@assays$RNA@scale.data
Counts <- Counts[order(rownames(Counts)), ]
PCA <- t(Counts) %*% Loadings
colnames(PCA) <- paste0("PCA_", 1:50)
Immune_R3[["pca"]] <- CreateDimReducObject(embeddings = PCA, key = "PCA_", assay = DefaultAssay(Immune_R3))  

# Harmonize the R3 PCs across diets. NOTE: This step is non-deterministic! Results vary from run-to-run!
Immune_R3 <- RunHarmony(object = Immune_R3, group.by.vars = "Diet")

# Predict UMAP coordinates for R3 using the model from the main dataset
tmp2 <- predict(UModel, Embeddings(Immune_R3, "harmony")[, c(1:16)])
colnames(tmp2) <- paste0("UMAP_", 1:2)
Immune_R3[["umap"]] <- CreateDimReducObject(embeddings = tmp2, key = "UMAP_", assay = DefaultAssay(Immune_R3))

# For each nuclei in R3, use the predicted UMAP coordinates to find the closest nuclei in the main data in UMAP space
IDs <- data.frame(ID = nn2(data = UModel$layout, query = tmp2, k = 1)$nn.idx)
IDs$cluster <- 10
for (i in 1:nrow(IDs)) {
  IDs[i, 2] <- as.character(Immune@meta.data[IDs[i, 1], "seurat_clusters"])
}

# Transfer the labels to the R3 object
Immune_R3$seurat_clusters <- IDs$cluster

# Set additional annotation
Labels <- Immune@meta.data[ ,c("Subtype","seurat_clusters","Label")]
Labels <- Labels[ duplicated(Labels$seurat_clusters)==F,]
for (i in 1:nrow(Labels)) {
  Immune_R3@meta.data[ Immune_R3@meta.data$seurat_clusters %in% Labels[i,"seurat_clusters"],"Subtype"] <- Labels[i,"Subtype"]
  Immune_R3@meta.data[ Immune_R3@meta.data$seurat_clusters %in% Labels[i,"seurat_clusters"],"Label"] <- Labels[i,"Label"]
}

# Transfer all subtype labels to overall object
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Immune_R3@meta.data[ Immune_R3@meta.data$Subtype == "Macrophages",]),"Subtype"] <- "Macrophages"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Immune_R3@meta.data[ Immune_R3@meta.data$Subtype == "DCs",]),"Subtype"] <- "DCs"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Immune_R3@meta.data[ Immune_R3@meta.data$Subtype == "B-cells",]),"Subtype"] <- "B-cells"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Immune_R3@meta.data[ Immune_R3@meta.data$Subtype == "T-cells",]),"Subtype"] <- "T-cells"

# Save the objects
saveRDS(Immune_R3, "eWAT_R3_Immune.Rds") ## This object is downloadable from Open Science Framework
saveRDS(eWAT_R3, "eWAT_R3.Rds") ## This object is downloadable from Open Science Framework, with futher annotations as subsequent scripts add to the file

## Cell type proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(Immune_R3$Subtype)), nrow=2))
colnames(Composition) <- unique(Immune_R3$Subtype)
for (i in 1:length(unique(Immune_R3$Subtype))) {
  Composition[1,i] <- nrow(Immune_R3@meta.data[ Immune_R3@meta.data$Dataset == "LFD_R3" & Immune_R3@meta.data$Subtype == colnames(Composition)[i],]) / nrow(eWAT_R3@meta.data[ eWAT_R3@meta.data$Dataset == "LFD_R3",])
  Composition[2,i] <- nrow(Immune_R3@meta.data[ Immune_R3@meta.data$Dataset == "HFD_R3" & Immune_R3@meta.data$Subtype == colnames(Composition)[i],]) / nrow(eWAT_R3@meta.data[ eWAT_R3@meta.data$Dataset == "HFD_R3",])
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
Clust5 <- Markers[[5]] 
Clust6 <- Markers[[6]] 
Clust7 <- Markers[[7]]
Clust8 <- Markers[[8]] 
Clust9 <- Markers[[9]]

# Subset to genes that are enriched or exclusive (!= NS)
Clust1 <- Clust1[Clust1$Marker != "NS", ]
Clust2 <- Clust2[Clust2$Marker != "NS", ]
Clust3 <- Clust3[Clust3$Marker != "NS", ]
Clust4 <- Clust4[Clust4$Marker != "NS", ]
Clust5 <- Clust5[Clust5$Marker != "NS", ]
Clust6 <- Clust6[Clust6$Marker != "NS", ]
Clust7 <- Clust7[Clust7$Marker != "NS", ]
Clust8 <- Clust8[Clust6$Marker != "NS", ]
Clust9 <- Clust9[Clust7$Marker != "NS", ]

# Extracting the minimum logFC across replicates in pairwise cluster comparisons
# RATIONALE: The genes with the highest minimum fold change are the most specific ones for any given cluster
Clust1$logFC_OvO <- apply(Clust1[, grep("logFC_Cluster", colnames(Clust1))], 1, FUN = "min")
Clust2$logFC_OvO <- apply(Clust2[, grep("logFC_Cluster", colnames(Clust2))], 1, FUN = "min")
Clust3$logFC_OvO <- apply(Clust3[, grep("logFC_Cluster", colnames(Clust3))], 1, FUN = "min")
Clust4$logFC_OvO <- apply(Clust4[, grep("logFC_Cluster", colnames(Clust4))], 1, FUN = "min")
Clust5$logFC_OvO <- apply(Clust5[, grep("logFC_Cluster", colnames(Clust5))], 1, FUN = "min")
Clust6$logFC_OvO <- apply(Clust6[, grep("logFC_Cluster", colnames(Clust6))], 1, FUN = "min")
Clust7$logFC_OvO <- apply(Clust7[, grep("logFC_Cluster", colnames(Clust7))], 1, FUN = "min")
Clust8$logFC_OvO <- apply(Clust8[, grep("logFC_Cluster", colnames(Clust8))], 1, FUN = "min")
Clust9$logFC_OvO <- apply(Clust9[, grep("logFC_Cluster", colnames(Clust9))], 1, FUN = "min")

# Computing gene module scores using the top 50 most specific marker genes. NOTE: THIS IS NON-DETERMINISTIC. THE RESULTS VARY SLIGHTLY FROM RUN-TO-RUN.
Immune_R3 <- AddModuleScore(Immune_R3, features = list(Cluster1 = Clust1[order(-Clust1$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster2 = Clust2[order(-Clust2$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster3 = Clust3[order(-Clust3$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster4 = Clust4[order(-Clust4$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster5 = Clust5[order(-Clust5$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster6 = Clust6[order(-Clust6$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster7 = Clust7[order(-Clust7$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster8 = Clust8[order(-Clust8$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster9 = Clust9[order(-Clust9$logFC_OvO), ][1:50, ]$Symbol))

# Scale the module scores (across all nuclei)
Immune_R3$Cluster1 <- scale(Immune_R3$Cluster1)
Immune_R3$Cluster2 <- scale(Immune_R3$Cluster2)
Immune_R3$Cluster3 <- scale(Immune_R3$Cluster3)
Immune_R3$Cluster4 <- scale(Immune_R3$Cluster4)
Immune_R3$Cluster5 <- scale(Immune_R3$Cluster5)
Immune_R3$Cluster6 <- scale(Immune_R3$Cluster6)
Immune_R3$Cluster7 <- scale(Immune_R3$Cluster7)
Immune_R3$Cluster8 <- scale(Immune_R3$Cluster8)
Immune_R3$Cluster9 <- scale(Immune_R3$Cluster9)

# Computing mean scaled module scores in each cluster in each replicate
Averages <- as.data.frame(matrix(ncol=length(Markers), nrow=length(Markers)))

# Set column and row names
Labels <- Immune_R3@meta.data[, c("Subtype","Label")]
Labels <- Labels[ duplicated(Labels$Label)==F,]
Labels <- Labels[ order(Labels$Label),]
Labels[1:6,1] <- paste("Macrophages",seq(1,6,by=1), sep="")
colnames(Averages) <- paste(Labels$Subtype,"module",sep="_")
rownames(Averages) <- Labels$Subtype

# Calculate averages
for (module in 1:length(Markers)) {
  for (label in 1:length(Markers)) {
      Averages[label,module] <- mean(Immune_R3@meta.data[ Immune_R3@meta.data$Label == label,paste("Cluster",module,sep="")])
  }
}

## Plot mean scaled module scores
col <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(as.matrix(Averages), cluster_columns = F, cluster_rows = F, col = col)

### Analysis of macrophages
# Subset
Macrophages_R3 <- subset(Immune_R3, subset = Subtype == "Macrophages")

# Manually set macrophage subtypes (not present in the saved object)
Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label %in% 1, "Subtype"] <- "PVM"
Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label %in% 2, "Subtype"] <- "LAM"
Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label %in% 3, "Subtype"] <- "NPVM"
Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label %in% 4, "Subtype"] <- "CEM"
Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label %in% 5, "Subtype"] <- "P-LAM"
Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label %in% 6, "Subtype"] <- "RM"
Macrophages_R3 <- SetIdent(Macrophages_R3, value = Macrophages_R3$Subtype)

# Transfer detailed subtype labels to overall object
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Subtype == "PVM",]),"Subtype"] <- "PVM"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Subtype == "LAM",]),"Subtype"] <- "LAM"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Subtype == "NPVM",]),"Subtype"] <- "NPVM"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Subtype == "CEM",]),"Subtype"] <- "CEM"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Subtype == "P-LAM",]),"Subtype"] <- "P-LAM"
eWAT_R3@meta.data[ rownames(eWAT_R3@meta.data) %in% rownames(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Subtype == "RM",]),"Subtype"] <- "RM"

# Save the objects
saveRDS(eWAT_R3, "eWAT_R3.Rds") ## This object is downloadable from Open Science Framework, with futher annotations as subsequent scripts add to the file

## DimPlots
DimPlot(Macrophages_R3, group.by="Subtype", label=T)
DimPlot(Macrophages_R3, group.by="Subtype", label=T, split.by="Diet")

## Cell type proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(Macrophages_R3$Subtype)), nrow=2))
colnames(Composition) <- unique(Macrophages_R3$Subtype)
for (i in 1:length(unique(Macrophages_R3$Subtype))) {
  Composition[1,i] <- nrow(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Dataset == "LFD_R3" & Macrophages_R3@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Dataset == "LFD_R3",])
  Composition[2,i] <- nrow(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Dataset == "HFD_R3" & Macrophages_R3@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Macrophages_R3@meta.data[ Macrophages_R3@meta.data$Dataset == "HFD_R3",])
}
rownames(Composition) <- c("LFD","HFD")

# Plot it
barplot(as.matrix(Composition), beside=T, las=1, ylim = c(0, 0.8), xlab = "Subtypes", ylab = "Fraction of nuclei")

## Marker gene heatmap
# Select marker genes for each cell type
PVM <- data.frame(gene = c("Mrc1","Lyve1","Cd163"), type = "PVM")
LAM <- data.frame(gene = c("Lpl","Trem2","Cd9"), type = "LAM")
NPVM <- data.frame(gene = c("Cd74","Ear2","Fcrls"), type = "NPVM")
CEM <- data.frame(gene = c("Col5a2","Tgfbr3","Col3a1"), type = "CEM")
PLAM <- data.frame(gene = c("Kif15","Kif11","Pola1"), type = "PLAM") 
RM <- data.frame(gene = c("Prg4", "Tgfb2","Ltbp1"), type = "RM")

# Combine results
Averages <- rbind(PVM,LAM,NPVM,CEM, PLAM, RM)

# Calculate average expression for selected marker genes
for (i in 1:nrow(Averages)) {
  Averages[i, 3] <- log1p(mean(expm1(Macrophages_R3@assays$RNA@data[rownames(Macrophages_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages_R3@assays$RNA@data) %in% rownames(Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label == "1", ])])))
  Averages[i, 4] <- log1p(mean(expm1(Macrophages_R3@assays$RNA@data[rownames(Macrophages_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages_R3@assays$RNA@data) %in% rownames(Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label == "2", ])])))
  Averages[i, 5] <- log1p(mean(expm1(Macrophages_R3@assays$RNA@data[rownames(Macrophages_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages_R3@assays$RNA@data) %in% rownames(Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label == "3", ])])))
  Averages[i, 6] <- log1p(mean(expm1(Macrophages_R3@assays$RNA@data[rownames(Macrophages_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages_R3@assays$RNA@data) %in% rownames(Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label == "4", ])])))
  Averages[i, 7] <- log1p(mean(expm1(Macrophages_R3@assays$RNA@data[rownames(Macrophages_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages_R3@assays$RNA@data) %in% rownames(Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label == "5", ])])))
  Averages[i, 8] <- log1p(mean(expm1(Macrophages_R3@assays$RNA@data[rownames(Macrophages_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages_R3@assays$RNA@data) %in% rownames(Macrophages_R3@meta.data[Macrophages_R3@meta.data$Label == "6", ])])))
}
rownames(Averages) <- Averages[, 1]

# Plot the heatmap
Heatmap(t(scale(t(Averages[, 3:ncol(Averages)]))), cluster_columns = F, cluster_rows = F)

## Diet-dependent gene regulation
# Import clustering results
Clusters <- readRDS("eWAT_Macrophage_Diet_Clusters.Rds") # Generated during analysis of the main datasets

# Extract diet-regulated genes from PVMs and NPVMs
NPVM_Diet <- Markers[[1]]
NPVM_Diet <- NPVM_Diet[NPVM_Diet$Diet != "NS", ]
PVM_Diet <- Markers[[3]]
PVM_Diet <- PVM_Diet[PVM_Diet$Diet != "NS", ]

# Calculate average expression for diet-regulated genes
Averages <- as.data.frame(matrix(ncol = 1, nrow = length(unique(sort(c(NPVM_Diet[, 1], PVM_Diet[, 1]))))))
Averages[, 1] <- unique(sort(c(NPVM_Diet[, 1], PVM_Diet[, 1])))
for (i in 1:nrow(Averages)) {
  Averages[i, 2] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "1" & Immune@meta.data$Diet == "LFD" & Immune@meta.data$Replicate == "R1", ])])))
  Averages[i, 3] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "1" & Immune@meta.data$Diet == "LFD" & Immune@meta.data$Replicate == "R2", ])])))
  Averages[i, 4] <- log1p(mean(expm1(Immune_R3@assays$RNA@data[rownames(Immune_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Immune_R3@assays$RNA@data) %in% rownames(Immune_R3@meta.data[Immune_R3@meta.data$Label == "1" & Immune_R3@meta.data$Diet == "LFD", ])])))
  Averages[i, 5] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "1" & Immune@meta.data$Diet == "HFD" & Immune@meta.data$Replicate == "R1", ])])))
  Averages[i, 6] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "1" & Immune@meta.data$Diet == "HFD" & Immune@meta.data$Replicate == "R2", ])])))
  Averages[i, 7] <- log1p(mean(expm1(Immune_R3@assays$RNA@data[rownames(Immune_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Immune_R3@assays$RNA@data) %in% rownames(Immune_R3@meta.data[Immune_R3@meta.data$Label == "1" & Immune_R3@meta.data$Diet == "HFD", ])])))
  Averages[i, 8] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "3" & Immune@meta.data$Diet == "LFD" & Immune@meta.data$Replicate == "R1", ])])))
  Averages[i, 9] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "3" & Immune@meta.data$Diet == "LFD" & Immune@meta.data$Replicate == "R2", ])])))
  Averages[i, 10] <- log1p(mean(expm1(Immune_R3@assays$RNA@data[rownames(Immune_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Immune_R3@assays$RNA@data) %in% rownames(Immune_R3@meta.data[Immune_R3@meta.data$Label == "3" & Immune_R3@meta.data$Diet == "LFD", ])])))
  Averages[i, 11] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "3" & Immune@meta.data$Diet == "HFD" & Immune@meta.data$Replicate == "R1", ])])))
  Averages[i, 12] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "3" & Immune@meta.data$Diet == "HFD" & Immune@meta.data$Replicate == "R2", ])])))
  Averages[i, 13] <- log1p(mean(expm1(Immune_R3@assays$RNA@data[rownames(Immune_R3@assays$RNA@data) %in% Averages[i, 1], colnames(Immune_R3@assays$RNA@data) %in% rownames(Immune_R3@meta.data[Immune_R3@meta.data$Label == "3" & Immune_R3@meta.data$Diet == "HFD", ])])))
}
rownames(Averages) <- Averages[, 1]

# Scale each experiment
Exp1 <- Averages[, c(2, 5, 8, 11)]
Exp1 <- t(scale(t(Exp1)))
Exp2 <- Averages[, c(3, 6, 9, 12)]
Exp2 <- t(scale(t(Exp2)))
Exp3 <- Averages[, c(4, 7, 10, 13)]
Exp3 <- t(scale(t(Exp3)))
Averages[, c(2, 5, 8, 11)] <- Exp1
Averages[, c(3, 6, 9, 12)] <- Exp2
Averages[, c(4, 7, 10, 13)] <- Exp3

# Set NAs to 0 to avoid errors
for (i in 2:ncol(Averages)) { Averages[ is.na(Averages[,i]),i] <- 0 }

# Plot the heatmap with order and split based on the original clusters
Lens <- unlist(lapply(Clusters, FUN="length"))
SplitVec <- c()
for (i in 1:length(Lens)) { SplitVec <- c(SplitVec, rep(paste("Cluster",i, sep=""), Lens[i])) }
Heatmap(as.matrix(Averages[unlist(Clusters), 2:ncol(Averages)]), show_row_names = F, cluster_columns = F, cluster_rows = F, row_split = SplitVec)
```
[Back to start](../README.md)<br>