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
LFD_R3 <- readRDS("~/LFD_R3.dgecounts.rds") # Equivalent to count matrices on GEO
HFD_R3 <- readRDS("~/HFD_R3.dgecounts.rds") # Equivalent to count matrices on GEO

# Extracting intron+exons counts 
LFD_R3 <- LFD_R3$umicount$inex$all
HFD_R3 <- HFD_R3$umicount$inex$all 

## Fill in the matrices to give them all the same dimensions
## RATIONALE: Genes with 0 counts across all barcodes in a particular experiment are left out from zUMIs.
# Find all non-zero genes in all conditions
Genes <- unique(c(rownames(LFD_R3),rownames(HFD_R3)))

# Insert empty lines with missing genes to get same dimensions
Tmp <- as.sparse(matrix(ncol=ncol(LFD_R3), nrow=length(Genes[!(Genes %in% rownames(LFD_R3))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(LFD_R3))]
LFD_R3 <- Matrix::rbind2(LFD_R3, Tmp)
LFD_R3 <- LFD_R3[ order(rownames(LFD_R3)),]

Tmp <- as.sparse(matrix(ncol=ncol(HFD_R3), nrow=length(Genes[!(Genes %in% rownames(HFD_R3))]), data=0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(HFD_R3))]
HFD_R3 <- Matrix::rbind2(HFD_R3, Tmp)
HFD_R3 <- HFD_R3[ order(rownames(HFD_R3)),]

## Paste in the experiment name into the column names to make barcodes/cell IDs unique.
colnames(LFD_R3) <- paste("LFD_R3_",colnames(LFD_R3), sep="")
colnames(HFD_R3) <- paste("HFD_R3_",colnames(HFD_R3), sep="")

## Convert Ensemble IDs to Gene Symbols
# Process gene list (downloaded from BioMart), keep only non-empty gene symbols and deduplicate.
Genes <- read.delim("Ensembl_Mouse.txt")
Genes <- Genes[ Genes$Gene.stable.ID %in% rownames(LFD_R3),]
Genes <- Genes[ Genes$MGI.symbol !=  "",]
Genes <- Genes[ duplicated(Genes$Gene.stable.ID) == F,]
Genes <- Genes[ duplicated(Genes$MGI.symbol) == F,]
Genes <- Genes[ order(Genes$Gene.stable.ID),]

# Subset count matrices and replace Ensemble IDs with gene symbols
LFD_R3 <- LFD_R3[ rownames(LFD_R3) %in% Genes$Gene.stable.ID,]
LFD_R3 <- LFD_R3[ order(rownames(LFD_R3)),]
rownames(LFD_R3) <- as.character(Genes$MGI.symbol)

HFD_R3 <- HFD_R3[ rownames(HFD_R3) %in% Genes$Gene.stable.ID,]
HFD_R3 <- HFD_R3[ order(rownames(HFD_R3)),]
rownames(HFD_R3) <- as.character(Genes$MGI.symbol)

## Match gene names with the main data
# Import the dataset
eWAT <- readRDS("eWAT_Annotated.Rds")

# Extract the genes
Genes <- rownames(eWAT)

# Insert empty lines when missing genes to get same dimensions across datasets
Tmp <- as.sparse(matrix(ncol = ncol(LFD_R3), nrow = length(Genes[!(Genes %in% rownames(LFD_R3))]), data = 0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(LFD_R3))]
LFD_R3 <- Matrix::rbind2(LFD_R3, Tmp)
LFD_R3 <- LFD_R3[order(rownames(LFD_R3)), ]

Tmp <- as.sparse(matrix(ncol = ncol(HFD_R3), nrow = length(Genes[!(Genes %in% rownames(HFD_R3))]), data = 0))
rownames(Tmp) <- Genes[!(Genes %in% rownames(HFD_R3))]
HFD_R3 <- Matrix::rbind2(HFD_R3, Tmp)
HFD_R3 <- HFD_R3[ order(rownames(HFD_R3)),]

## Create SingleCellExperiment objects to use DropUtils, scran, scater, etc.
LFD_R3_sce <- SingleCellExperiment(list(counts=LFD_R3))
HFD_R3_sce <- SingleCellExperiment(list(counts=HFD_R3))

## Find empty droplets. NOTE: THIS STEP IS NON-DETERMINSTIC - RESULTS VARY FROM RUN TO RUN
## Lower threshold were set based on being in the range of knee points and adjusted to match the expectation of recovering approximately 10.000 non-empty droplets, based on the loading of the 10X chip.
e.out <- emptyDrops(counts(LFD_R3_sce), lower=1000) 
LFD_R3_sce <- LFD_R3_sce[,which(e.out$FDR <= 0.05)]
e.out <- emptyDrops(counts(HFD_R3_sce), lower=1000) 
HFD_R3_sce <- HFD_R3_sce[,which(e.out$FDR <= 0.05)] 

## Calculate QC parameters (throws a warning, that can be ignored)
LFD_R3_sce <- calculateQCMetrics(LFD_R3_sce, feature_controls=list(Mito=grep(pattern = "^mt-", x = rownames(LFD_R3_sce), value = FALSE)))
HFD_R3_sce <- calculateQCMetrics(HFD_R3_sce, feature_controls=list(Mito=grep(pattern = "^mt-", x = rownames(HFD_R3_sce), value = FALSE)))

### Histograms of quality measures for each dataset
hist(R3@meta.data[ R3@meta.data$Dataset == "LFD_R3", "log10_total_counts"], breaks = 50)
hist(R3@meta.data[ R3@meta.data$Dataset == "HFD_R3", "log10_total_counts"], breaks = 50)
hist(R3@meta.data[ R3@meta.data$Dataset == "LFD_R3", "total_features_by_counts"], breaks = 20)
hist(R3@meta.data[ R3@meta.data$Dataset == "HFD_R3", "total_features_by_counts"], breaks = 20)
hist(R3@meta.data[ R3@meta.data$Dataset == "LFD_R3", "pct_counts_Mito"], breaks = 20)
hist(R3@meta.data[ R3@meta.data$Dataset == "HFD_R3", "pct_counts_Mito"], breaks = 20)

## Threshold filtering of droplets in each dataset
## REMOVE: Droplets with more than 15% mitochondrial reads, less than 1000 UMIs, less than 500 genes or extremely high ratio between counts and genes (low complexity))
LFD_R3_sce <- LFD_R3_sce[,!(LFD_R3_sce$pct_counts_Mito > 15)] 
LFD_R3_sce <- LFD_R3_sce[,!(LFD_R3_sce$total_counts/LFD_R3_sce$total_features_by_counts > 2.5)] 
LFD_R3_sce <- LFD_R3_sce[,(LFD_R3_sce$total_counts >= 1000 & LFD_R3_sce$total_features_by_counts >= 500)] 

HFD_R3_sce <- HFD_R3_sce[,!(HFD_R3_sce$pct_counts_Mito > 15)]
HFD_R3_sce <- HFD_R3_sce[,!(HFD_R3_sce$total_counts/HFD_R3_sce$total_features_by_counts > 2.5)] 
HFD_R3_sce <- HFD_R3_sce[,(HFD_R3_sce$total_counts >= 1000 & HFD_R3_sce$total_features_by_counts >= 500)] 

## Automatic filtering droplets in each dataset using PCA across all QC metrics
# Calculate outliers
LFD_R3_sce <- runPCA(LFD_R3_sce, use_coldata=TRUE, detect_outliers=TRUE)
HFD_R3_sce <- runPCA(HFD_R3_sce, use_coldata=TRUE, detect_outliers=TRUE)

# Remove outliers with PCA
LFD_R3_sce <- LFD_R3_sce[ ,!LFD_R3_sce$outlier] 
HFD_R3_sce <- HFD_R3_sce[ ,!HFD_R3_sce$outlier] 

# Filter the sce objects according to genes in the main dataset
LFD_R3_sce <- LFD_R3_sce[which(rownames(LFD_R3_sce) %in% rownames(eWAT)) ,]
HFD_R3_sce <- HFD_R3_sce[which(rownames(HFD_R3_sce) %in% rownames(eWAT)), ]

## Normalize the count matrices. 
# Cluster each data. NOTE: THIS STEP IS NON-DETERMINISTIC - RESULTS VARY FROM RUN TO RUN
LFD_R3_clusters <- quickCluster(LFD_R3_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
HFD_R3_clusters <- quickCluster(HFD_R3_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

# Compute scaling factors
LFD_R3_sce <- computeSumFactors(LFD_R3_sce, min.mean=0.1, cluster=LFD_R3_clusters)
HFD_R3_sce <- computeSumFactors(HFD_R3_sce, min.mean=0.1, cluster=HFD_R3_clusters)

# Normalize the counts
LFD_R3_sce <- normalize(LFD_R3_sce)
HFD_R3_sce <- normalize(HFD_R3_sce)

## Calculate doublet scores. NOTE: THIS STEP IS NON-DETERMINISTI - RESULTS VARY FROM RUN TO RUN
LFD_R3_sce$DoubletScore <- doubletCells(LFD_R3_sce, BSPARAM=IrlbaParam())
HFD_R3_sce$DoubletScore <- doubletCells(HFD_R3_sce, BSPARAM=IrlbaParam())

# Rescaling batch effects with batchelor
rescaled <- multiBatchNorm(LFD_R3_sce, HFD_R3_sce)

### QC by clustering - Individual datasets. NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
### RATIONaLE: Low quality nuclei may be included due to threshold effects. Deep clustering can help to reveal if there are groups of low quality nuclei that does not mix with the remaining nuclei, and thus can be removed.
## Create Seurat objects
LFD_R3_seurat <- as.Seurat(rescaled[[1]], counts = "counts", data = "logcounts")
HFD_R3_seurat <- as.Seurat(rescaled[[2]], counts = "counts", data = "logcounts")

## LFD_R3
# Iteration 1
LFD_R3_seurat <- ScaleData(LFD_R3_seurat)
LFD_R3_seurat <- RunPCA(LFD_R3_seurat)
LFD_R3_seurat <- RunUMAP(LFD_R3_seurat, dims = 1:20, reduction = "pca")
LFD_R3_seurat <- FindNeighbors(object = LFD_R3_seurat, dims = 1:20, reduction = "pca")
LFD_R3_seurat <- FindClusters(LFD_R3_seurat, resolution = 6, algorithm = 1)
VlnPlot(LFD_R3_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 17, 28, 30, 35, 36, 37, and 39 to 50 have high doublet scores, high mitochondrial percents or low counts/genes
LFD_R3_seurat <- subset(LFD_R3_seurat, cells = rownames(LFD_R3_seurat@meta.data[!(LFD_R3_seurat@meta.data$seurat_clusters %in% c(17, 28, 30, 35, 36:37, 39:50)), ]))

# Iteration 2
VariableFeatures(LFD_R3_seurat) <- rownames(LFD_R3_seurat)
LFD_R3_seurat <- ScaleData(LFD_R3_seurat)
LFD_R3_seurat <- RunPCA(LFD_R3_seurat)
LFD_R3_seurat <- RunUMAP(LFD_R3_seurat, dims = 1:20, reduction = "pca")
LFD_R3_seurat <- FindNeighbors(object = LFD_R3_seurat, dims = 1:20, reduction = "pca")
LFD_R3_seurat <- FindClusters(LFD_R3_seurat, resolution = 5, algorithm = 1)
VlnPlot(LFD_R3_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No more clusters to remove

## HFD_R3
# Iteration 1
VariableFeatures(HFD_R3_seurat) <- rownames(HFD_R3_seurat)
HFD_R3_seurat <- ScaleData(HFD_R3_seurat)
HFD_R3_seurat <- RunPCA(HFD_R3_seurat)
HFD_R3_seurat <- RunUMAP(HFD_R3_seurat, dims = 1:20, reduction = "pca")
HFD_R3_seurat <- FindNeighbors(object = HFD_R3_seurat, dims = 1:20, reduction = "pca")
HFD_R3_seurat <- FindClusters(HFD_R3_seurat, resolution = 6, algorithm = 1)
VlnPlot(HFD_R3_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 14, 25, 29, 44, 45, 49, and 51 to 53 have high doublet scores, high mitochondrial percents or low counts/genes
HFD_R3_seurat <- subset(HFD_R3_seurat, cells = rownames(HFD_R3_seurat@meta.data[!(HFD_R3_seurat@meta.data$seurat_clusters %in% c(14, 25, 29, 44:45, 49, 51:53)), ]))

# Iteration 2
HFD_R3_seurat <- ScaleData(HFD_R3_seurat)
HFD_R3_seurat <- RunPCA(HFD_R3_seurat)
HFD_R3_seurat <- RunUMAP(HFD_R3_seurat, dims = 1:20, reduction ="pca")
HFD_R3_seurat <- FindNeighbors(object = HFD_R3_seurat, dims = 1:20, reduction = "pca")
HFD_R3_seurat <- FindClusters(HFD_R3_seurat, resolution = 5, algorithm = 1)
VlnPlot(HFD_R3_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 18, 28, 29, 32, 33, 34 and 35 have high doublet scores, high mitochondrial percents or low counts/genes
HFD_R3_seurat <- subset(HFD_R3_seurat, cells = rownames(HFD_R3_seurat@meta.data[!(HFD_R3_seurat@meta.data$seurat_clusters %in% c(18, 28:29, 32:35)), ]))

# Iteration 3
HFD_R3_seurat <- ScaleData(HFD_R3_seurat)
HFD_R3_seurat <- RunPCA(HFD_R3_seurat)
HFD_R3_seurat <- RunUMAP(HFD_R3_seurat, dims = 1:20, reduction="pca")
HFD_R3_seurat <- FindNeighbors(object = HFD_R3_seurat, dims = 1:20, reduction = "pca")
HFD_R3_seurat <- FindClusters(HFD_R3_seurat, resolution = 4, algorithm = 1)
VlnPlot(HFD_R3_seurat, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# No more clusters to remove

### QC by clustering - Replicates NOTE: RERUNNING THE CODE MAY PRODUCE DIFFERENT CLUSTERS AND CLUSTER LABELS DUE TO NON-DETERMINISTIC STEPS
# Merge the datasets
R3 <- merge(LFD_R3_seurat, y = HFD_R3_seurat)
R3$Dataset <- substr(rownames(R3@meta.data), 0, 11)

# Iteration 1
VariableFeatures(R3) <- rownames(R3)
R3 <- ScaleData(R3)
R3 <- RunPCA(R3)
R3 <- RunHarmony(R3, group.by.vars = "Dataset", dims.use = 1:20)
R3 <- RunUMAP(R3, dims = 1:20, reduction = "harmony")
R3 <- FindNeighbors(object = R3, dims = 1:20, reduction = "harmony")
R3 <- FindClusters(R3, resolution = 6, algorithm = 1)
VlnPlot(R3, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 40, 43, 45, 46 and 47 have high doublet scores, high mitochondrial percents or low counts/genes
R3 <- subset(R3, cells = rownames(R3@meta.data[!(R3@meta.data$seurat_clusters %in% c(40, 43, 45:47)), ]))

# Iteration 2
R3 <- ScaleData(R3)
R3 <- RunPCA(R3)
R3 <- RunHarmony(R3, group.by.vars = "Dataset", dims.use = 1:20)
R3 <- RunUMAP(R3, dims = 1:20, reduction = "harmony")
R3 <- FindNeighbors(object = R3, dims = 1:20, reduction = "harmony")
R3 <- FindClusters(R3, resolution = 5, algorithm = 1)
VlnPlot(R3, c("nFeature_RNA", "nCount_RNA", "DoubletScore","pct_counts_Mito"), pt.size=0)
# Clusters 41, 42 and 43 have high doublet scores, high mitochondrial percents or low counts/genes
R3 <- subset(R3, cells = rownames(R3@meta.data[!(R3@meta.data$seurat_clusters %in% c(41:43)), ]))

# Iteration 3
R3 <- ScaleData(R3)
R3 <- RunPCA(R3)
R3 <- RunHarmony(R3, group.by.vars = "Dataset", dims.use = 1:20)
R3 <- RunUMAP(R3, dims = 1:20, reduction = "harmony")
R3 <- FindNeighbors(object = R3, dims = 1:20, reduction = "harmony")
R3 <- FindClusters(R3, resolution = 4, algorithm = 1)
# No more clusters to remove

### Subset each SCE object to the identified high-quality droplets. 
LFD_R3_sce <- LFD_R3_sce[,which(colnames(LFD_R3_sce) %in% colnames(R3))]
HFD_R3_sce <- HFD_R3_sce[,which(colnames(HFD_R3_sce) %in% colnames(R3))]

### Normalize the datasets. NOTE: THIS STEP IS NON-DETERMINISTIC - RESULTS VARY FROM RUN TO RUN
# Cluster the nuclei in each dataset
LFD_R3_clusters <- quickCluster(LFD_R3_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())
HFD_R3_clusters <- quickCluster(HFD_R3_sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

# Compute scaling factors
LFD_R3_sce <- computeSumFactors(LFD_R3_sce, min.mean=0.1, cluster=LFD_R3_clusters)
HFD_R3_sce <- computeSumFactors(HFD_R3_sce, min.mean=0.1, cluster=HFD_R3_clusters)

# Normalize the counts
LFD_R3_sce <- normalize(LFD_R3_sce)
HFD_R3_sce <- normalize(HFD_R3_sce)

### Reduce batch effects by rescaling across datasets
# Normalize across samples
rescaled <- batchelor::multiBatchNorm(
  LFD_R3_sce, 
  HFD_R3_sce
)

### Merge all the data and embed
# Create seurat objects
LFD_R3_seurat <- as.Seurat(rescaled[[1]], counts = "counts", data = "logcounts")
HFD_R3_seurat <- as.Seurat(rescaled[[2]], counts = "counts", data = "logcounts")

# Merge the datasets
R3 <- merge(LFD_R3_seurat, HFD_R3_seurat)
R3$Dataset <- substr(rownames(R3@meta.data),0,6)
R3$Replicate <- substr(rownames(R3@meta.data),5,6)
R3$Diet <- substr(rownames(R3@meta.data),0,3)
R3$Annotation <- "NA"
R3$Subtype <- "NA"
R3$Label <- "NA"
R3@project.name <- "Replicate 3"

# Save the object
saveRDS(R3, "eWAT_R3.Rds") ## This object is downloadable from Open Science Framework, with futher annotations as subsequent scripts add to the file
```
[Back to start](../README.md)<br>