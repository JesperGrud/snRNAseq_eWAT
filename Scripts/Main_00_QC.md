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

## Import count objects generated using zUMIs
# Import data
LFD_R1 <- readRDS("~/LFD_R1.dgecounts.rds") # Equivalent to count matrices on GEO
LFD_R2 <- readRDS("~/LFD_R2.dgecounts.rds") # Equivalent to count matrices on GEO
HFD_R1 <- readRDS("~/HFD_R1.dgecounts.rds") # Equivalent to count matrices on GEO
HFD_R2 <- readRDS("~/HFD_R2.dgecounts.rds") # Equivalent to count matrices on GEO

# Extracting intron+exons counts 
LFD_R1 <- LFD_R1$umicount$inex$all
LFD_R2 <- LFD_R2$umicount$inex$all 
HFD_R1 <- HFD_R1$umicount$inex$all 
HFD_R2 <- HFD_R2$umicount$inex$all

## Fill in the matrices to give them all the same dimensions
## RATIONALE: Genes with 0 counts across all barcodes in a particular experiment are left out from zUMIs.
# Find all non-zero genes in all conditions
Genes <- unique(c(rownames(LFD_R1),rownames(LFD_R2),rownames(HFD_R1),rownames(HFD_R2)))

# Insert empty lines with missing genes to get same dimensions
Tmp <- as.sparse(matrix(ncol=ncol(LFD_R1), nrow=length(Genes[!(Genes %in% rownames(LFD_R1))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(LFD_R1))]
LFD_R1 <- Matrix::rbind2(LFD_R1, Tmp)
LFD_R1 <- LFD_R1[ order(rownames(LFD_R1)),]
dim(LFD_R1) # 29638 492758

Tmp <- as.sparse(matrix(ncol=ncol(LFD_R2), nrow=length(Genes[!(Genes %in% rownames(LFD_R2))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(LFD_R2))]
LFD_R2 <- Matrix::rbind2(LFD_R2, Tmp)
LFD_R2 <- LFD_R2[ order(rownames(LFD_R2)),]
dim(LFD_R2) # 29638 643085

Tmp <- as.sparse(matrix(ncol=ncol(HFD_R1), nrow=length(Genes[!(Genes %in% rownames(HFD_R1))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(HFD_R1))]
HFD_R1 <- Matrix::rbind2(HFD_R1, Tmp)
HFD_R1 <- HFD_R1[ order(rownames(HFD_R1)),]
dim(HFD_R1) # 29638 595784

Tmp <- as.sparse(matrix(ncol=ncol(HFD_R2), nrow=length(Genes[!(Genes %in% rownames(HFD_R2))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(HFD_R2))]
HFD_R2 <- Matrix::rbind2(HFD_R2, Tmp)
HFD_R2 <- HFD_R2[ order(rownames(HFD_R2)),]
dim(HFD_R2) # 29638 639691

## Paste in the experiment name into the column names to make barcodes/cell IDs unique.
colnames(LFD_R1) <- paste("LFD_R1_",colnames(LFD_R1), sep="")
colnames(LFD_R2) <- paste("LFD_R2_",colnames(LFD_R2), sep="")
colnames(HFD_R1) <- paste("HFD_R1_",colnames(HFD_R1), sep="")
colnames(HFD_R2) <- paste("HFD_R2_",colnames(HFD_R2), sep="")

## Convert Ensemble IDs to Gene Symbols
# Process gene list (downloaded from BioMart), keep only non-empty gene symbols and deduplicate.
Genes <- read.delim("Ensembl_Mouse.txt")
Genes <- Genes[ Genes$Gene.stable.ID %in% rownames(LFD_R1),]
Genes <- Genes[ Genes$MGI.symbol !=  "",]
Genes <- Genes[ duplicated(Genes$Gene.stable.ID) == F,]
Genes <- Genes[ duplicated(Genes$MGI.symbol) == F,]
Genes <- Genes[ order(Genes$Gene.stable.ID),]
dim(Genes) # 29286     7

# Subset count matrices and replace Ensemble IDs with gene symbols
LFD_R1 <- LFD_R1[ rownames(LFD_R1) %in% Genes$Gene.stable.ID,]
LFD_R1 <- LFD_R1[ order(rownames(LFD_R1)),]
rownames(LFD_R1) <- as.character(Genes$MGI.symbol)
dim(LFD_R1) # 29286 492758

LFD_R2 <- LFD_R2[ rownames(LFD_R2) %in% Genes$Gene.stable.ID,]
LFD_R2 <- LFD_R2[ order(rownames(LFD_R2)),]
rownames(LFD_R2) <- as.character(Genes$MGI.symbol)
dim(LFD_R2) # 29286 643085

HFD_R1 <- HFD_R1[ rownames(HFD_R1) %in% Genes$Gene.stable.ID,]
HFD_R1 <- HFD_R1[ order(rownames(HFD_R1)),]
rownames(HFD_R1) <- as.character(Genes$MGI.symbol)
dim(HFD_R1) # 29286 595784

HFD_R2 <- HFD_R2[ rownames(HFD_R2) %in% Genes$Gene.stable.ID,]
HFD_R2 <- HFD_R2[ order(rownames(HFD_R2)),]
rownames(HFD_R2) <- as.character(Genes$MGI.symbol)
dim(HFD_R2) # 29286 639691

## Create SingleCellExperiment objects to use DropUtils, scran, scater, etc.
LFD_R1_sce <- SingleCellExperiment(list(counts=LFD_R1))
LFD_R2_sce <- SingleCellExperiment(list(counts=LFD_R2))
HFD_R1_sce <- SingleCellExperiment(list(counts=HFD_R1))
HFD_R2_sce <- SingleCellExperiment(list(counts=HFD_R2))

## Find empty droplets. NOTE: THIS STEP IS NON-DETERMINSTIC - RESULTS VARY FROM RUN TO RUN
## Lower threshold were set based on being in the range of knee points and adjusted to match the expectation of recovering approximately 10.000 non-empty droplets, based on the loading of the 10X chip.
e.out <- emptyDrops(counts(LFD_R1_sce), lower=400) 
LFD_R1_sce <- LFD_R1_sce[,which(e.out$FDR <= 0.05)]
e.out <- emptyDrops(counts(LFD_R2_sce), lower=400) 
LFD_R2_sce <- LFD_R2_sce[,which(e.out$FDR <= 0.05)] 
e.out <- emptyDrops(counts(HFD_R1_sce), lower=700) 
HFD_R1_sce <- HFD_R1_sce[,which(e.out$FDR <= 0.05)] 
e.out <- emptyDrops(counts(HFD_R2_sce), lower=650) 
HFD_R2_sce <- HFD_R2_sce[,which(e.out$FDR <= 0.05)] 

## Calculate QC parameters (throws a warning, that can be ignored)
LFD_R1_sce <- calculateQCMetrics(LFD_R1_sce, feature_controls=list(Mito=grep(pattern = "^mt-", x = rownames(LFD_R1_sce), value = FALSE)))
LFD_R2_sce <- calculateQCMetrics(LFD_R2_sce, feature_controls=list(Mito=grep(pattern = "^mt-", x = rownames(LFD_R2_sce), value = FALSE)))
HFD_R1_sce <- calculateQCMetrics(HFD_R1_sce, feature_controls=list(Mito=grep(pattern = "^mt-", x = rownames(HFD_R1_sce), value = FALSE)))
HFD_R2_sce <- calculateQCMetrics(HFD_R2_sce, feature_controls=list(Mito=grep(pattern = "^mt-", x = rownames(HFD_R2_sce), value = FALSE)))

### Histograms of quality measures for each dataset
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R1", "log10_total_counts"], breaks = 50)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R2", "log10_total_counts"], breaks = 50)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R1", "log10_total_counts"], breaks = 50)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R2", "log10_total_counts"], breaks = 50)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R1", "total_features_by_counts"], breaks = 20)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R2", "total_features_by_counts"], breaks = 20)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R1", "total_features_by_counts"], breaks = 20)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R2", "total_features_by_counts"], breaks = 20)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R1", "pct_counts_Mito"], breaks = 20)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R2", "pct_counts_Mito"], breaks = 20)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R1", "pct_counts_Mito"], breaks = 20)
hist(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R2", "pct_counts_Mito"], breaks = 20)

## Threshold filtering of droplets in each dataset
## REMOVE: Droplets with more than 15% mitochondrial reads, less than 1000 UMIs, less than 500 genes or extremely high ratio between counts and genes (low complexity))
LFD_R1_sce <- LFD_R1_sce[,!(LFD_R1_sce$pct_counts_Mito > 15)] 
LFD_R1_sce <- LFD_R1_sce[,!(LFD_R1_sce$total_counts/LFD_R1_sce$total_features_by_counts > 2.5)] 
LFD_R1_sce <- LFD_R1_sce[,(LFD_R1_sce$total_counts >= 1000 & LFD_R1_sce$total_features_by_counts >= 500)] 

LFD_R2_sce <- LFD_R2_sce[,!(LFD_R2_sce$pct_counts_Mito > 15)] 
LFD_R2_sce <- LFD_R2_sce[,!(LFD_R2_sce$total_counts/LFD_R2_sce$total_features_by_counts > 2.5)] 
LFD_R2_sce <- LFD_R2_sce[,(LFD_R2_sce$total_counts >= 1000 & LFD_R2_sce$total_features_by_counts >= 500)] 

HFD_R1_sce <- HFD_R1_sce[,!(HFD_R1_sce$pct_counts_Mito > 15)]
HFD_R1_sce <- HFD_R1_sce[,!(HFD_R1_sce$total_counts/HFD_R1_sce$total_features_by_counts > 2.5)] 
HFD_R1_sce <- HFD_R1_sce[,(HFD_R1_sce$total_counts >= 1000 & HFD_R1_sce$total_features_by_counts >= 500)] 

HFD_R2_sce <- HFD_R2_sce[,!(HFD_R2_sce$pct_counts_Mito > 15)] 
HFD_R2_sce <- HFD_R2_sce[,!(HFD_R2_sce$total_counts/HFD_R2_sce$total_features_by_counts > 2.5)] 
HFD_R2_sce <- HFD_R2_sce[,(HFD_R2_sce$total_counts >= 1000 & HFD_R2_sce$total_features_by_counts >= 500)] 

## Automatic filtering droplets in each dataset using PCA across all QC metrics
# Calculate outliers
LFD_R1_sce <- runPCA(LFD_R1_sce, use_coldata=TRUE, detect_outliers=TRUE)
LFD_R2_sce <- runPCA(LFD_R2_sce, use_coldata=TRUE, detect_outliers=TRUE)
HFD_R1_sce <- runPCA(HFD_R1_sce, use_coldata=TRUE, detect_outliers=TRUE)
HFD_R2_sce <- runPCA(HFD_R2_sce, use_coldata=TRUE, detect_outliers=TRUE)

# Remove outliers with PCA
LFD_R1_sce <- LFD_R1_sce[ ,!LFD_R1_sce$outlier] 
LFD_R2_sce <- LFD_R2_sce[ ,!LFD_R2_sce$outlier] 
HFD_R1_sce <- HFD_R1_sce[ ,!HFD_R1_sce$outlier] 
HFD_R2_sce <- HFD_R2_sce[ ,!HFD_R2_sce$outlier] 

## Threshold filtering of genes in each dataset
## REMOVE: Genes expressed in less than 10 nuclei in all datasets
# Find lowly expressed genes and get the intersection
LFD_R1_low <- names(which(nexprs(LFD_R1_sce, byrow=T) <= 10))
LFD_R2_low <- names(which(nexprs(LFD_R2_sce, byrow=T) <= 10))
HFD_R1_low <- names(which(nexprs(HFD_R1_sce, byrow=T) <= 10))
HFD_R2_low <- names(which(nexprs(HFD_R2_sce, byrow=T) <= 10))
Low <- Reduce(intersect, list(LFD_R1_low,LFD_R2_low,HFD_R1_low,HFD_R2_low))

# Remove lowly expressed genes
LFD_R1_sce <- LFD_R1_sce[ which(!(rownames(LFD_R1_sce) %in% Low)),] 
LFD_R2_sce <- LFD_R2_sce[ which(!(rownames(LFD_R2_sce) %in% Low)),] 
HFD_R1_sce <- HFD_R1_sce[ which(!(rownames(HFD_R1_sce) %in% Low)),] 
HFD_R2_sce <- HFD_R2_sce[ which(!(rownames(HFD_R2_sce) %in% Low)),] 

## Detect ambient genes
## RATIONALE: Ambient genes are detected in droplets that do not contain nuclei. 
ExprsFun <- function(x) { sum(x > 0) }

Empty <- LFD_R1[, colnames(LFD_R1) %in% names(which(Matrix::colSums(LFD_R1) >= 1 & Matrix::colSums(LFD_R1) <= 10)),] # Find droplets that are empty, but has recovered genes (less than 10 UMI counts, but more than 1)
dim(Empty) # 29286 347222
Empty <- Empty[ rownames(Empty) %in% names(which(Matrix::rowSums(Empty) >= 50)),] # Find genes with at least 50 UMI counts in 'empty' droplets
dim(Empty) # 3561 347222
Empty <- as.matrix(t(Empty))
Empty <- apply(Empty,2,ExprsFun)/nrow(Empty)
LFD_R1_empty <- sort(names(which(Empty >= 0.001))) # Find genes expressed in at least 0.1% of the empty droplets. 

Empty <- LFD_R2[, colnames(LFD_R2) %in% names(which(Matrix::colSums(LFD_R2) >= 1 & Matrix::colSums(LFD_R2) <= 10)),]
dim(Empty) # 29286 326155
Empty <- Empty[ rownames(Empty) %in% names(which(Matrix::rowSums(Empty) >= 50)),]
dim(Empty) # 3088 326155
Empty <- as.matrix(t(Empty))
Empty <- apply(Empty,2,ExprsFun)/nrow(Empty)
LFD_R2_empty <- sort(names(which(Empty >= 0.001))) 

Empty <- HFD_R1[, colnames(HFD_R1) %in% names(which(Matrix::colSums(HFD_R1) >= 1 & Matrix::colSums(HFD_R1) <= 10)),]
dim(Empty) # 29286 454271
Empty <- Empty[ rownames(Empty) %in% names(which(Matrix::rowSums(Empty) >= 50)),]
dim(Empty) # 4205 454271
Empty <- as.matrix(t(Empty))
Empty <- apply(Empty,2,ExprsFun)/nrow(Empty)
HFD_R1_empty <- sort(names(which(Empty >= 0.001))) 

Empty <- HFD_R2[, colnames(HFD_R2) %in% names(which(Matrix::colSums(HFD_R2) >= 1 & Matrix::colSums(HFD_R2) <= 10)),]
dim(Empty) # 29286 347222
Empty <- Empty[ rownames(Empty) %in% names(which(Matrix::rowSums(Empty) >= 50)),]
dim(Empty) # 3364 347222
Empty <- as.matrix(t(Empty))
Empty <- apply(Empty,2,ExprsFun)/nrow(Empty)
HFD_R2_empty <- sort(names(which(Empty >= 0.001))) 

# Get a list of genes that are ambient in anyone dataset
Amb <- unique(c(LFD_R1_empty,LFD_R2_empty,HFD_R1_empty,HFD_R2_empty)) 

## Filtering genes based on biotype and transcript level support
## REMOVE: Non-protein coding genes, keep the one with high transcript level support (tsl)
Genes <- read.delim("Ensembl_Mouse.txt")
Genes <- Genes[ Genes$MGI.symbol %in% rownames(LFD_R2_sce),]
Genes <- Genes[ grep("tsl1", Genes$Transcript.support.level..TSL.),]
Genes <- Genes[ Genes$MGI.symbol !=  "",]
Genes <- Genes[ Genes$Transcript.type == "protein_coding",]
Genes <- Genes[ Genes$Gene.type == "protein_coding",]
Genes <- Genes[ !is.na(Genes$MGI.symbol),]
Genes <- Genes[ duplicated(Genes$Gene.stable.ID) == F,] 

# Filter the sce objects and ambient genes
LFD_R1_sce <- LFD_R1_sce[ which(rownames(LFD_R1_sce) %in% Genes$MGI.symbol),] # 13718 genes left
LFD_R2_sce <- LFD_R2_sce[ which(rownames(LFD_R2_sce) %in% Genes$MGI.symbol),] # 13718 genes left
HFD_R1_sce <- HFD_R1_sce[ which(rownames(HFD_R1_sce) %in% Genes$MGI.symbol),] # 13718 genes left
HFD_R2_sce <- HFD_R2_sce[ which(rownames(HFD_R2_sce) %in% Genes$MGI.symbol),] # 13718 genes left
Amb <- Amb[Amb %in% rownames(LFD_R1_sce)] 
 
# Save the ambient RNA object
saveRDS(Amb, "Ambient.Rds") # Available as Supplemental Table S1 or from the Open Science Framework

## Normalize the count matrices. 
# Cluster each data. NOTE: THIS STEP IS NON-DETERMINISTIC - RESULTS VARY FROM RUN TO RUN
LFD_R1_clusters <- quickCluster(LFD_R1_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
LFD_R2_clusters <- quickCluster(LFD_R2_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
HFD_R1_clusters <- quickCluster(HFD_R1_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
HFD_R2_clusters <- quickCluster(HFD_R2_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

# Compute scaling factors
LFD_R1_sce <- computeSumFactors(LFD_R1_sce, min.mean=0.1, cluster=LFD_R1_clusters)
LFD_R2_sce <- computeSumFactors(LFD_R2_sce, min.mean=0.1, cluster=LFD_R2_clusters)
HFD_R1_sce <- computeSumFactors(HFD_R1_sce, min.mean=0.1, cluster=HFD_R1_clusters)
HFD_R2_sce <- computeSumFactors(HFD_R2_sce, min.mean=0.1, cluster=HFD_R2_clusters)

# Normalize the counts
LFD_R1_sce <- normalize(LFD_R1_sce)
LFD_R2_sce <- normalize(LFD_R2_sce)
HFD_R1_sce <- normalize(HFD_R1_sce)
HFD_R2_sce <- normalize(HFD_R2_sce)

## Calculate doublet scores. NOTE: THIS STEP IS NON-DETERMINISTI - RESULTS VARY FROM RUN TO RUN
LFD_R1_sce$DoubletScore <- doubletCells(LFD_R1_sce, BSPARAM=IrlbaParam())
LFD_R2_sce$DoubletScore <- doubletCells(LFD_R2_sce, BSPARAM=IrlbaParam())
HFD_R1_sce$DoubletScore <- doubletCells(HFD_R1_sce, BSPARAM=IrlbaParam())
HFD_R2_sce$DoubletScore <- doubletCells(HFD_R2_sce, BSPARAM=IrlbaParam())

### QC by clustering - Individual datasets. NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
### RATIONaLE: Low quality nuclei may be included due to threshold effects. Deep clustering can help to reveal if there are groups of low quality nuclei that does not mix with the remaining nuclei, and thus can be removed.
## Create Seurat objects
LFD_R1_seurat <- as.Seurat(LFD_R1_sce, counts = "counts", data = "logcounts")
LFD_R2_seurat <- as.Seurat(LFD_R2_sce, counts = "counts", data = "logcounts")
HFD_R1_seurat <- as.Seurat(HFD_R1_sce, counts = "counts", data = "logcounts")
HFD_R2_seurat <- as.Seurat(HFD_R2_sce, counts = "counts", data = "logcounts")

## LFD_R1
# Iteration 1 - Use all genes
VariableFeatures(LFD_R1_seurat) <- rownames(LFD_R1_seurat)
LFD_R1_seurat <- ScaleData(LFD_R1_seurat)
LFD_R1_seurat <- RunPCA(LFD_R1_seurat)
LFD_R1_seurat <- RunUMAP(LFD_R1_seurat, dims=1:20, reduction="pca")
LFD_R1_seurat <- FindNeighbors(object = LFD_R1_seurat, dims = 1:20, reduction = "pca")
LFD_R1_seurat <- FindClusters(LFD_R1_seurat, resolution = 6, algorithm = 1)
VlnPlot(LFD_R1_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 37, 42, 46, 47 and 50 have high doublet scores, high mitochondrial percents or low counts/genes
LFD_R1_seurat <- subset(LFD_R1_seurat, cells=rownames(LFD_R1_seurat@meta.data[ !(LFD_R1_seurat@meta.data$seurat_clusters %in% c(37,42,46,47,50)),]))

# Iteration 2 - Exclude ambient genes
VariableFeatures(LFD_R1_seurat) <- rownames(LFD_R1_seurat)[!(rownames(LFD_R1_seurat) %in% Amb)]
LFD_R1_seurat <- ScaleData(LFD_R1_seurat)
LFD_R1_seurat <- RunPCA(LFD_R1_seurat)
LFD_R1_seurat <- RunUMAP(LFD_R1_seurat, dims=1:20, reduction="pca")
LFD_R1_seurat <- FindNeighbors(object = LFD_R1_seurat, dims = 1:20, reduction = "pca")
LFD_R1_seurat <- FindClusters(LFD_R1_seurat, resolution = 5, algorithm = 1)
VlnPlot(LFD_R1_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 39 and 41 have high doublet scores, high mitochondrial percents or low counts/genes
LFD_R1_seurat <- subset(LFD_R1_seurat, cells=rownames(LFD_R1_seurat@meta.data[ !(LFD_R1_seurat@meta.data$seurat_clusters %in% c(39,41)),]))

# Iteration 3
VariableFeatures(LFD_R1_seurat) <- rownames(LFD_R1_seurat)[!(rownames(LFD_R1_seurat) %in% Amb)]
LFD_R1_seurat <- ScaleData(LFD_R1_seurat)
LFD_R1_seurat <- RunPCA(LFD_R1_seurat)
LFD_R1_seurat <- RunUMAP(LFD_R1_seurat, dims=1:20, reduction="pca")
LFD_R1_seurat <- FindNeighbors(object = LFD_R1_seurat, dims = 1:20, reduction = "pca")
LFD_R1_seurat <- FindClusters(LFD_R1_seurat, resolution = 5, algorithm = 1)
VlnPlot(LFD_R1_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No more clusters to remove

## LFD_R2
# Iteration 1
VariableFeatures(LFD_R2_seurat) <- rownames(LFD_R2_seurat)
LFD_R2_seurat <- ScaleData(LFD_R2_seurat)
LFD_R2_seurat <- RunPCA(LFD_R2_seurat)
LFD_R2_seurat <- RunUMAP(LFD_R2_seurat, dims=1:20, reduction="pca")
LFD_R2_seurat <- FindNeighbors(object = LFD_R2_seurat, dims = 1:20, reduction = "pca")
LFD_R2_seurat <- FindClusters(LFD_R2_seurat, resolution = 7, algorithm = 1)
VlnPlot(LFD_R2_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 38, 48, 51, 53, 55 and 56 have high doublet scores, high mitochondrial percents or low counts/genes
LFD_R2_seurat <- subset(LFD_R2_seurat, cells=rownames(LFD_R2_seurat@meta.data[ !(LFD_R2_seurat@meta.data$seurat_clusters %in% c(38,48,51,53,55,56)),]))

# Iteration 2
VariableFeatures(LFD_R2_seurat) <- rownames(LFD_R2_seurat)[!(rownames(LFD_R2_seurat) %in% Amb)]
LFD_R2_seurat <- ScaleData(LFD_R2_seurat)
LFD_R2_seurat <- RunPCA(LFD_R2_seurat)
LFD_R2_seurat <- RunUMAP(LFD_R2_seurat, dims=1:20, reduction="pca")
LFD_R2_seurat <- FindNeighbors(object = LFD_R2_seurat, dims = 1:20, reduction = "pca")
LFD_R2_seurat <- FindClusters(LFD_R2_seurat, resolution = 5, algorithm = 1)
VlnPlot(LFD_R2_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 40 and 41 have high doublet scores, high mitochondrial percents or low counts/genes
LFD_R2_seurat <- subset(LFD_R2_seurat, cells=rownames(LFD_R2_seurat@meta.data[ !(LFD_R2_seurat@meta.data$seurat_clusters %in% c(40,41)),]))

# Iteration 3
VariableFeatures(LFD_R2_seurat) <- rownames(LFD_R2_seurat)[!(rownames(LFD_R2_seurat) %in% Amb)]
LFD_R2_seurat <- ScaleData(LFD_R2_seurat)
LFD_R2_seurat <- RunPCA(LFD_R2_seurat)
LFD_R2_seurat <- RunUMAP(LFD_R2_seurat, dims=1:20, reduction="pca")
LFD_R2_seurat <- FindNeighbors(object = LFD_R2_seurat, dims = 1:20, reduction = "pca")
LFD_R2_seurat <- FindClusters(LFD_R2_seurat, resolution = 5, algorithm = 1)
VlnPlot(LFD_R2_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No more clusters to remove

## HFD_R1
# Iteration 1
VariableFeatures(HFD_R1_seurat) <- rownames(HFD_R1_seurat)
HFD_R1_seurat <- ScaleData(HFD_R1_seurat)
HFD_R1_seurat <- RunPCA(HFD_R1_seurat)
HFD_R1_seurat <- RunUMAP(HFD_R1_seurat, dims=1:20, reduction="pca")
HFD_R1_seurat <- FindNeighbors(object = HFD_R1_seurat, dims = 1:20, reduction = "pca")
HFD_R1_seurat <- FindClusters(HFD_R1_seurat, resolution = 7, algorithm = 1)
VlnPlot(HFD_R1_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 4, 34 and 44 have high doublet scores, high mitochondrial percents or low counts/genes
HFD_R1_seurat <- subset(HFD_R1_seurat, cells=rownames(HFD_R1_seurat@meta.data[ !(HFD_R1_seurat@meta.data$seurat_clusters %in% c(4,34,44)),]))

# Iteration 2
VariableFeatures(HFD_R1_seurat) <- rownames(HFD_R1_seurat)[!(rownames(HFD_R1_seurat) %in% Amb)]
HFD_R1_seurat <- ScaleData(HFD_R1_seurat)
HFD_R1_seurat <- RunPCA(HFD_R1_seurat)
HFD_R1_seurat <- RunUMAP(HFD_R1_seurat, dims=1:20, reduction="pca")
HFD_R1_seurat <- FindNeighbors(object = HFD_R1_seurat, dims = 1:20, reduction = "pca")
HFD_R1_seurat <- FindClusters(HFD_R1_seurat, resolution = 5, algorithm = 1)
VlnPlot(HFD_R1_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 5 and 34 have high doublet scores, high mitochondrial percents or low counts/genes
HFD_R1_seurat <- subset(HFD_R1_seurat, cells=rownames(HFD_R1_seurat@meta.data[ !(HFD_R1_seurat@meta.data$seurat_clusters %in% c(5,34)),]))

# Iteration 3
VariableFeatures(HFD_R1_seurat) <- rownames(HFD_R1_seurat)
HFD_R1_seurat <- ScaleData(HFD_R1_seurat)
HFD_R1_seurat <- RunPCA(HFD_R1_seurat)
HFD_R1_seurat <- RunUMAP(HFD_R1_seurat, dims=1:20, reduction="pca")
HFD_R1_seurat <- FindNeighbors(object = HFD_R1_seurat, dims = 1:20, reduction = "pca")
HFD_R1_seurat <- FindClusters(HFD_R1_seurat, resolution = 5, algorithm = 1)
VlnPlot(HFD_R1_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 25 and 30 have high doublet scores, high mitochondrial percents or low counts/genes
HFD_R1_seurat <- subset(HFD_R1_seurat, cells=rownames(HFD_R1_seurat@meta.data[ !(HFD_R1_seurat@meta.data$seurat_clusters %in% c(25,30)),]))

# Iteration 4
HFD_R1_seurat <- ScaleData(HFD_R1_seurat)
HFD_R1_seurat <- RunPCA(HFD_R1_seurat)
HFD_R1_seurat <- RunUMAP(HFD_R1_seurat, dims=1:20, reduction="pca")
HFD_R1_seurat <- FindNeighbors(object = HFD_R1_seurat, dims = 1:20, reduction = "pca")
HFD_R1_seurat <- FindClusters(HFD_R1_seurat, resolution = 5, algorithm = 1)
VlnPlot(HFD_R1_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No more clusters to remove

## HFD_R2
# Iteration 1
VariableFeatures(HFD_R2_seurat) <- rownames(HFD_R2_seurat)
HFD_R2_seurat <- ScaleData(HFD_R2_seurat)
HFD_R2_seurat <- RunPCA(HFD_R2_seurat)
HFD_R2_seurat <- RunUMAP(HFD_R2_seurat, dims=1:20, reduction="pca")
HFD_R2_seurat <- FindNeighbors(object = HFD_R2_seurat, dims = 1:20, reduction = "pca")
HFD_R2_seurat <- FindClusters(HFD_R2_seurat, resolution = 7, algorithm = 1)
VlnPlot(HFD_R2_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 0, 23, 46, 44, 49 and 50 have high doublet scores, high mitochondrial percents or low counts/genes
HFD_R2_seurat <- subset(HFD_R2_seurat, cells=rownames(HFD_R2_seurat@meta.data[ !(HFD_R2_seurat@meta.data$seurat_clusters %in% c(0,23,44,46,49,50)),]))

# Iteration 2
VariableFeatures(HFD_R2_seurat) <- rownames(HFD_R2_seurat)[!(rownames(HFD_R2_seurat) %in% Amb)]
HFD_R2_seurat <- ScaleData(HFD_R2_seurat)
HFD_R2_seurat <- RunPCA(HFD_R2_seurat)
HFD_R2_seurat <- RunUMAP(HFD_R2_seurat, dims=1:20, reduction="pca")
HFD_R2_seurat <- FindNeighbors(object = HFD_R2_seurat, dims = 1:20, reduction = "pca")
HFD_R2_seurat <- FindClusters(HFD_R2_seurat, resolution = 5, algorithm = 1)
VlnPlot(HFD_R2_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 26, 35, 38, 39, 41 have high doublet scores, high mitochondrial percents or low counts/genes
HFD_R2_seurat <- subset(HFD_R2_seurat, cells=rownames(HFD_R2_seurat@meta.data[ !(HFD_R2_seurat@meta.data$seurat_clusters %in% c(26,35,38,39,41)),]))

# Iteration 3
VariableFeatures(HFD_R2_seurat) <- rownames(HFD_R2_seurat)[!(rownames(HFD_R2_seurat) %in% Amb)]
HFD_R2_seurat <- ScaleData(HFD_R2_seurat)
HFD_R2_seurat <- RunPCA(HFD_R2_seurat)
HFD_R2_seurat <- RunUMAP(HFD_R2_seurat, dims=1:20, reduction="pca")
HFD_R2_seurat <- FindNeighbors(object = HFD_R2_seurat, dims = 1:20, reduction = "pca")
HFD_R2_seurat <- FindClusters(HFD_R2_seurat, resolution = 5, algorithm = 1)
VlnPlot(HFD_R2_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No more clusters to remove

### QC by clustering - Replicates NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
# Merge the datasets
R1 <- merge(LFD_R1_seurat, y=HFD_R1_seurat)
R1$Dataset <- substr(rownames(R1@meta.data),0,11)

# Iteration 1. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
VariableFeatures(R1) <- rownames(R1)[!(rownames(R1) %in% Amb)]
R1 <- ScaleData(R1)
R1 <- RunPCA(R1)
R1 <- RunHarmony(R1, group.by.vars="Dataset", dims.use=1:20)
R1 <- RunUMAP(R1, dims=1:20, reduction="harmony")
R1 <- FindNeighbors(object = R1, dims = 1:20, reduction = "harmony")
R1 <- FindClusters(R1, resolution = 5, algorithm = 1)
VlnPlot(R1, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Cluster 48 have high doublet scores, high mitochondrial percents or low counts/genes
R1 <- subset(R1, cells=rownames(R1@meta.data[ !(R1@meta.data$seurat_clusters %in% c(48)),]))

# Iteration 2. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
R1 <- ScaleData(R1)
R1 <- RunPCA(R1)
R1 <- RunHarmony(R1, group.by.vars="Dataset", dims.use=1:20)
R1 <- RunUMAP(R1, dims=1:20, reduction="harmony")
R1 <- FindNeighbors(object = R1, dims = 1:20, reduction = "harmony")
R1 <- FindClusters(R1, resolution = 5, algorithm = 1)
VlnPlot(R1, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No more clusters to remove

## R2
# Merge the datasets
R2 <- merge(LFD_R2_seurat, y=HFD_R2_seurat)
R2$Dataset <- substr(rownames(R2@meta.data),0,11)

# Iteration 1. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
VariableFeatures(R2) <- rownames(R2)[!(rownames(R2) %in% Amb)]
R2 <- ScaleData(R2)
R2 <- RunPCA(R2)
R2 <- RunHarmony(R2, group.by.vars="Dataset", dims.use=1:20)
R2 <- RunUMAP(R2, dims=1:20, reduction="harmony")
R2 <- FindNeighbors(object = R2, dims = 1:20, reduction = "harmony")
R2 <- FindClusters(R2, resolution = 5, algorithm = 1)
VlnPlot(R2, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Cluster 41 have high doublet scores, high mitochondrial percents or low counts/genes
R2 <- subset(R2, cells=rownames(R2@meta.data[ !(R2@meta.data$seurat_clusters %in% c(41)),]))

# Iteration 2. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
R2 <- ScaleData(R2)
R2 <- RunPCA(R2)
R2 <- RunHarmony(R2, group.by.vars="Dataset", dims.use=1:20)
R2 <- RunUMAP(R2, dims=1:20, reduction="harmony")
R2 <- FindNeighbors(object = R2, dims = 1:20, reduction = "harmony")
R2 <- FindClusters(R2, resolution = 5, algorithm = 1)
VlnPlot(R2, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No more clusters to remove

### QC by clustering - Combined NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
# Merge the datasets
eWAT <- merge(R1, R2)
eWAT$Dataset <- substr(rownames(eWAT@meta.data),0,11)

# Iteration 1. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
VariableFeatures(eWAT) <- rownames(eWAT)
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
eWAT <- RunUMAP(eWAT, dims=1:20, reduction="harmony")
eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")
eWAT <- FindClusters(eWAT, resolution = 5, algorithm = 1)
VlnPlot(eWAT, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 54 and 55 have high doublet scores, high mitochondrial percents or low counts/genes
eWAT <- subset(eWAT, cells=rownames(eWAT@meta.data[ !(eWAT@meta.data$seurat_clusters %in% c(54,55)),]))

# Iteration 2. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
VariableFeatures(eWAT) <- rownames(eWAT)[!(rownames(eWAT) %in% Amb)]
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
eWAT <- RunUMAP(eWAT, dims=1:20, reduction="harmony")
eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")
eWAT <- FindClusters(eWAT, resolution = 5, algorithm = 1)
VlnPlot(eWAT, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 18, 24 and 50 have high doublet scores, high mitochondrial percents or low counts/genes
eWAT <- subset(eWAT, cells=rownames(eWAT@meta.data[ !(eWAT@meta.data$seurat_clusters %in% c(18,24,50)),]))

# Iteration 3. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
VariableFeatures(eWAT) <- rownames(eWAT)[!(rownames(eWAT) %in% Amb)]
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
eWAT <- RunUMAP(eWAT, dims=1:20, reduction="harmony")
eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")
eWAT <- FindClusters(eWAT, resolution = 5, algorithm = 1)
VlnPlot(eWAT, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Cluster 46 have high doublet scores, high mitochondrial percents or low counts/genes
eWAT <- subset(eWAT, cells=rownames(eWAT@meta.data[ !(eWAT@meta.data$seurat_clusters %in% c(46)),]))

# Iteration 4. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
eWAT <- RunUMAP(eWAT, dims=1:20, reduction="harmony")
eWAT <- FindNeighbors(object = eWAT, dims = 1:20, reduction = "harmony")
eWAT <- FindClusters(eWAT, resolution = 5, algorithm = 1)
VlnPlot(eWAT, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No more clusters to remove

### Subset each SCE object to the identified high-quality droplets. 
LFD_R1_sce <- LFD_R1_sce[,which(colnames(LFD_R1_sce) %in% colnames(eWAT))]
LFD_R2_sce <- LFD_R2_sce[,which(colnames(LFD_R2_sce) %in% colnames(eWAT))]
HFD_R1_sce <- HFD_R1_sce[,which(colnames(HFD_R1_sce) %in% colnames(eWAT))]
HFD_R2_sce <- HFD_R2_sce[,which(colnames(HFD_R2_sce) %in% colnames(eWAT))]

### Normalize the datasets. NOTE: THIS STEP IS NON-DETERMINISTIC - RESULTS VARY FROM RUN TO RUN
# Cluster the nuclei in each dataset
LFD_R1_clusters <- quickCluster(LFD_R1_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
LFD_R2_clusters <- quickCluster(LFD_R2_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
HFD_R1_clusters <- quickCluster(HFD_R1_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
HFD_R2_clusters <- quickCluster(HFD_R2_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

# Compute scaling factors
LFD_R1_sce <- computeSumFactors(LFD_R1_sce, min.mean=0.1, cluster=LFD_R1_clusters)
LFD_R2_sce <- computeSumFactors(LFD_R2_sce, min.mean=0.1, cluster=LFD_R2_clusters)
HFD_R1_sce <- computeSumFactors(HFD_R1_sce, min.mean=0.1, cluster=HFD_R1_clusters)
HFD_R2_sce <- computeSumFactors(HFD_R2_sce, min.mean=0.1, cluster=HFD_R2_clusters)

# Normalize the counts
LFD_R1_sce <- normalize(LFD_R1_sce)
LFD_R2_sce <- normalize(LFD_R2_sce)
HFD_R1_sce <- normalize(HFD_R1_sce)
HFD_R2_sce <- normalize(HFD_R2_sce)

### Reduce batch effects by rescaling across datasets
# Normalize across samples
rescaled <- batchelor::multiBatchNorm(
  LFD_R1_sce, 
  LFD_R2_sce,
  HFD_R1_sce, 
  HFD_R2_sce
)

### Merge all the data and embed
# Create seurat objects
LFD_R1_seurat <- as.Seurat(rescaled[[1]], counts = "counts", data = "logcounts")
LFD_R2_seurat <- as.Seurat(rescaled[[2]], counts = "counts", data = "logcounts")
HFD_R1_seurat <- as.Seurat(rescaled[[3]], counts = "counts", data = "logcounts")
HFD_R2_seurat <- as.Seurat(rescaled[[4]], counts = "counts", data = "logcounts")

# Merge the datasets
eWAT <- merge(LFD_R1_seurat, y=c(HFD_R1_seurat,LFD_R2_seurat,HFD_R2_seurat))
eWAT$Dataset <- substr(rownames(eWAT@meta.data),0,6)
eWAT$Replicate <- substr(rownames(eWAT@meta.data),5,6)
eWAT$Diet <- substr(rownames(eWAT@meta.data),0,3)

### Analyze effects of ambient gene removal and batch correction with Harmony integration
## Batch vector with numbers
Batch <- rep(1,nrow(eWAT@meta.data))
names(Batch) <- rownames(eWAT@meta.data)
Batch[ which(eWAT@meta.data$Replicate == "R2") ] <- 2

## Dont remove ambient genes and evaluate using PCA
# Decompose and embed
eWAT <- FindVariableFeatures(eWAT, nfeature=1000)
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunUMAP(eWAT, dims=1:15, reduction="pca")

# DimPlot (color by replicate)
DimPlot(eWAT, group.by="Replicate")

## Evaluate dataset integration
# Create a PC object from a mock object
testdata <- create_testset_multibatch(n.genes=1000,n.batch=3, plattform='any') # mock object
pca.data <- prcomp(testdata$data, center=TRUE)
pca.data$x <- Embeddings(eWAT, "pca")[,1:15]

# Calculate metrics
kBET.sd <- kBET(Embeddings(eWAT, "pca")[,1:15], batch = Batch, k0 = 20, do.pca = FALSE, plot=FALSE)
LISI.sd <- compute_lisi(Embeddings(eWAT, "pca")[,1:15], eWAT@meta.data, c('Replicate'), perplexity=40)
Sil.sd <- batch_sil(pca.data, Batch, nPCs=15)
PC.sd <- pcRegression(pca.data, Batch)

## Remove ambient genes and evaluate using PCA
# Decompose and embed
eWAT <- FindVariableFeatures(eWAT, nfeature=1000)
VariableFeatures(eWAT) <- VariableFeatures(eWAT)[!(VariableFeatures(eWAT) %in% Amb)]
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunUMAP(eWAT, dims=1:15, reduction="pca")

# DimPlot (color by replicate)
DimPlot(eWAT, group.by="Replicate")

## Evaluate dataset integration
# Create a PC object from a mock object
testdata <- create_testset_multibatch(n.genes=1000,n.batch=3, plattform='any') # mock object
pca.data <- prcomp(testdata$data, center=TRUE)
pca.data$x <- Embeddings(eWAT, "pca")[,1:15]

# Calculate metrics
kBET.amb <- kBET(Embeddings(eWAT, "pca")[,1:15], batch = Batch, k0 = 20, do.pca = FALSE, plot=FALSE)
LISI.amb <- compute_lisi(Embeddings(eWAT, "pca")[,1:15], eWAT@meta.data, c('Replicate'), perplexity=40)
Sil.amb <- batch_sil(pca.data, Batch, nPCs=15)
PC.amb <- pcRegression(pca.data, Batch)

## Remove ambient genes and evaluate using harmony. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
# Decompose and embed
eWAT <- FindVariableFeatures(eWAT, nfeature=1000)
VariableFeatures(eWAT) <- VariableFeatures(eWAT)[!(VariableFeatures(eWAT) %in% Amb)]
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
eWAT <- RunUMAP(eWAT, dims=1:15, reduction="harmony")

# DimPlot (color by replicate)
DimPlot(eWAT, group.by="Replicate")

## Evaluate dataset integration
# Create a PC object from a mock object
testdata <- create_testset_multibatch(n.genes=1000,n.batch=3, plattform='any') # mock object
pca.data <- prcomp(testdata$data, center=TRUE)
pca.data$x <- Embeddings(eWAT, "harmony")[,1:15]

# Calculate metrics
kBET.harm <- kBET(Embeddings(eWAT, "harmony")[,1:15], batch = Batch, k0 = 20, do.pca = FALSE, plot=FALSE)
LISI.harm <- compute_lisi(Embeddings(eWAT, "harmony")[,1:15], eWAT@meta.data, c('Replicate'), perplexity=40)
Sil.harm <- batch_sil(pca.data, Batch, nPCs=15)
PC.harm <- pcRegression(pca.data, Batch)

# Barplot the results
barplot(c(mean(LISI.sd[,1]),mean(LISI.amb[,1]),mean(LISI.harm[,1])), las=1, main="LISI", ylim=c(0,2))
barplot(c(mean(kBET.sd$stats$kBET.observed[,1]),mean(kBET.amb$stats$kBET.observed[,1]),mean(kBET.harm$stats$kBET.observed[,1])), las=1, main="kBET", ylim=c(0,1))
barplot(c(Sil.sd,Sil.amb,Sil.harm), las=1, main="Silhouette widt", ylim=c(0,0.06))
barplot(c(PC.sd$R2Var,PC.amb$R2Var,PC.harm$R2Var), las=1, main="PC Regression", ylim=c(0,0.06))

### Final embedding and clustering using ambient gene removal and Harmony (best performance in above test). NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
eWAT <- FindVariableFeatures(eWAT, nfeature=1000)
VariableFeatures(eWAT) <- VariableFeatures(eWAT)[!(VariableFeatures(eWAT) %in% Amb)]
eWAT <- ScaleData(eWAT)
eWAT <- RunPCA(eWAT)
eWAT <- RunHarmony(eWAT, group.by.vars="Dataset", dims.use=1:20)
eWAT <- RunUMAP(eWAT, dims=1:15, reduction="harmony")
eWAT <- FindNeighbors(object = eWAT, dims = 1:2, reduction = "umap")
eWAT <- FindClusters(eWAT, resolution = 0.005, algorithm = 1)

#### Quality control each detected cluster
### Cluster 0
## Process
# Subset
Subset <- subset(eWAT, subset = seurat_clusters == 0)

# Embed. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
Subset <- FindVariableFeatures(Subset, nfeatures=1000)
VariableFeatures(Subset) <- VariableFeatures(Subset)[!(VariableFeatures(Subset) %in% Amb)]
Subset <- ScaleData(Subset)
Subset <- RunPCA(Subset)
Subset <- RunHarmony(Subset, group.by.vars="Dataset")
Subset <- RunUMAP(Subset, dims=1:16, reduction="harmony")
Subset <- FindNeighbors(object = Subset, dims = 1:2, reduction = "umap", k.param=50)
Subset <- FindClusters(Subset, resolution = 0.1, algorithm = 1)

## Quality control
# Violin plots
VlnPlot(Subset, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No outliers / bad clusters found. No more filtering!

### Cluster 1 (Fibro-adipo progenitors, see Main_01_Overall.R)
## Process
# Subset
Subset <- subset(eWAT, subset = seurat_clusters == 1)

# Embed. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
Subset <- FindVariableFeatures(Subset, nfeatures=2000)
VariableFeatures(Subset) <- VariableFeatures(Subset)[!(VariableFeatures(Subset) %in% Amb)]
Subset <- ScaleData(Subset)
Subset <- RunPCA(Subset)
Subset <- RunHarmony(Subset, group.by.vars="Dataset")
Subset <- RunUMAP(Subset, dims=1:20, reduction="harmony")
Subset <- FindNeighbors(object = Subset, dims = 1:2, reduction = "umap", k.param=50)
Subset <- FindClusters(Subset, resolution = 0.1, algorithm = 1)

## Quality control
# Violin plots
VlnPlot(Subset, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)

# Remove high doublet cluster
Subset <- subset(Subset, subset = seurat_clusters != 4)

# Remove barcodes from the full dataset
eWAT <- subset(eWAT, cells = colnames(eWAT)[!(colnames(eWAT) %in% colnames(Subset))])

### Cluster 2  
## Process
# Subset
Subset <- subset(eWAT, subset = seurat_clusters == 2)

# Embed. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
Subset <- FindVariableFeatures(Subset, nfeatures=2000)
VariableFeatures(Subset) <- VariableFeatures(Subset)[!(VariableFeatures(Subset) %in% Amb)]
Subset <- ScaleData(Subset)
Subset <- RunPCA(Subset)
Subset <- RunHarmony(Subset, group.by.vars="Dataset")
Subset <- RunUMAP(Subset, dims=1:12, reduction="harmony")
Subset <- FindNeighbors(object = Subset, dims = 1:2, reduction = "umap", k.param=30)
Subset <- FindClusters(Subset, resolution = 0.1, algorithm = 1)

## Quality control
# Violin plots
VlnPlot(Subset, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)

# Remove high doublet clusters
Subset <- subset(Subset, subset = seurat_clusters != 4)
Subset <- subset(Subset, subset = seurat_clusters != 5)

# Remove barcodes from the full dataset
eWAT <- subset(eWAT, cells = colnames(eWAT)[!(colnames(eWAT) %in% colnames(Subset))])

### Cluster 3
## Process
# Subset
Subset <- subset(eWAT, subset = seurat_clusters == 3)

# Embed. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
Subset <- FindVariableFeatures(Subset, nfeatures=1000)
VariableFeatures(Subset) <- VariableFeatures(Subset)[!(VariableFeatures(Subset) %in% Amb)]
Subset <- ScaleData(Subset)
Subset <- RunPCA(Subset)
Subset <- RunHarmony(Subset, group.by.vars="Dataset")
Subset <- RunUMAP(Subset, dims=1:16, reduction="harmony")
Subset <- FindNeighbors(object = Subset, dims = 1:2, reduction = "umap", k.param=50)
Subset <- FindClusters(Subset, resolution = 0.1, algorithm = 1)

## Quality control
# Violin plots
VlnPlot(Subset, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No outliers / bad clusters found. No more filtering!

### Cluster 4
## Process
# Subset
Subset <- subset(eWAT, subset = seurat_clusters == 4)

# Embed. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
Subset <- FindVariableFeatures(Subset, nfeatures=1000)
VariableFeatures(Subset) <- VariableFeatures(Subset)[!(VariableFeatures(Subset) %in% Amb)]
Subset <- ScaleData(Subset)
Subset <- RunPCA(Subset)
Subset <- RunHarmony(Subset, group.by.vars="Dataset")
Subset <- RunUMAP(Subset, dims=1:16, reduction="harmony")
Subset <- FindNeighbors(object = Subset, dims = 1:2, reduction = "umap", k.param=50)
Subset <- FindClusters(Subset, resolution = 0.1, algorithm = 1)

## Quality control
# Violin plots
VlnPlot(Subset, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No outliers / bad clusters found. No more filtering!

### Cluster 5
## Process
# Subset
Subset <- subset(eWAT, subset = seurat_clusters == 5)

# Embed. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
Subset <- FindVariableFeatures(Subset, nfeatures=1000)
VariableFeatures(Subset) <- VariableFeatures(Subset)[!(VariableFeatures(Subset) %in% Amb)]
Subset <- ScaleData(Subset)
Subset <- RunPCA(Subset)
Subset <- RunHarmony(Subset, group.by.vars="Dataset")
Subset <- RunUMAP(Subset, dims=1:16, reduction="harmony")
Subset <- FindNeighbors(object = Subset, dims = 1:2, reduction = "umap", k.param=50)
Subset <- FindClusters(Subset, resolution = 0.1, algorithm = 1)

## Quality control
# Violin plots
VlnPlot(Subset, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No outliers / bad clusters found. No more filtering!

### Cluster 6
## Process
# Subset
Subset <- subset(eWAT, subset = seurat_clusters == 6)

# Embed. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
Subset <- FindVariableFeatures(Subset, nfeatures=1000)
VariableFeatures(Subset) <- VariableFeatures(Subset)[!(VariableFeatures(Subset) %in% Amb)]
Subset <- ScaleData(Subset)
Subset <- RunPCA(Subset)
Subset <- RunHarmony(Subset, group.by.vars="Dataset")
Subset <- RunUMAP(Subset, dims=1:16, reduction="harmony")
Subset <- FindNeighbors(object = Subset, dims = 1:2, reduction = "umap", k.param=50)
Subset <- FindClusters(Subset, resolution = 0.1, algorithm = 1)

## Quality control
# Violin plots
VlnPlot(Subset, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No outliers / bad clusters found. No more filtering!

### Cluster 7
## Process
# Subset
Subset <- subset(eWAT, subset = seurat_clusters == 7)

# Embed. NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN
Subset <- FindVariableFeatures(Subset, nfeatures=1000)
VariableFeatures(Subset) <- VariableFeatures(Subset)[!(VariableFeatures(Subset) %in% Amb)]
Subset <- ScaleData(Subset)
Subset <- RunPCA(Subset)
Subset <- RunHarmony(Subset, group.by.vars="Dataset")
Subset <- RunUMAP(Subset, dims=1:16, reduction="harmony")
Subset <- FindNeighbors(object = Subset, dims = 1:2, reduction = "umap", k.param=50)
Subset <- FindClusters(Subset, resolution = 0.1, algorithm = 1)

## Quality control
# Violin plots
VlnPlot(Subset, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No outliers / bad clusters found. No more filtering!

#### Save RDS object
saveRDS(eWAT, "eWAT_Annotated.Rds") ## This object is downloadable from Open Science Framework, with futher annotations as subsequent scripts add to the file

```
[Back to start](../README.md)<br>