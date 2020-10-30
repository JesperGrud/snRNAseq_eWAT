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

### Import data
Adipocytes <- readRDS("eWAT_Adipocytes.Rds")
FAP <- readRDS("eWAT_FAP.Rds")
Amb <- readRDS("Ambient.Rds")

### RNA Velocity in FAP compartment
# Importing raw zUMIs output
LFD_R1 <- readRDS("~/LFD_R1.dgecounts.rds")
LFD_R2 <- readRDS("~/LFD_R2.dgecounts.rds")
HFD_R1 <- readRDS("~/HFD_R1.dgecounts.rds")
HFD_R2 <- readRDS("~/HFD_R2.dgecounts.rds")

# Extract intron and exons counts 
LFD_R1_exon <- LFD_R1$umicount$exon$all
LFD_R2_exon <- LFD_R2$umicount$exon$all 
HFD_R1_exon <- HFD_R1$umicount$exon$all 
HFD_R2_exon <- HFD_R2$umicount$exon$all
LFD_R1_intron <- SM3869$umicount$intron$all
LFD_R2_intron <- SM4381$umicount$intron$all 
HFD_R1_intron <- SM3871$umicount$intron$all 
HFD_R2_intron <- SM4382$umicount$intron$all

# Subset to only QC-filtered cells
LFD_R1_exon <- LFD_R1_exon[, paste("LFD_R1_",colnames(LFD_R1_exon),sep="") %in% colnames(FAP)]
LFD_R2_exon <- LFD_R2_exon[, paste("LFD_R2_",colnames(LFD_R2_exon),sep="") %in% colnames(FAP)]
HFD_R1_exon <- HFD_R1_exon[, paste("HFD_R1_",colnames(HFD_R1_exon),sep="") %in% colnames(FAP)]
HFD_R2_exon <- HFD_R2_exon[, paste("HFD_R2_",colnames(HFD_R2_exon),sep="") %in% colnames(FAP)]
LFD_R1_intron <- LFD_R1_intron[, paste("LFD_R1_",colnames(LFD_R1_intron),sep="") %in% colnames(FAP)]
LFD_R2_intron <- LFD_R2_intron[, paste("LFD_R2_",colnames(LFD_R2_intron),sep="") %in% colnames(FAP)]
HFD_R1_intron <- HFD_R1_intron[, paste("HFD_R1_",colnames(HFD_R1_intron),sep="") %in% colnames(FAP)]
HFD_R2_intron <- HFD_R2_intron[, paste("HFD_R2_",colnames(HFD_R2_intron),sep="") %in% colnames(FAP)]

# Set column names
colnames(LFD_R1_exon) <- paste("LFD_R1_",colnames(LFD_R1_exon),sep="")
colnames(LFD_R2_exon) <- paste("LFD_R2_",colnames(LFD_R2_exon),sep="")
colnames(HFD_R1_exon) <- paste("HFD_R1_",colnames(HFD_R1_exon),sep="")
colnames(HFD_R2_exon) <- paste("HFD_R2_",colnames(HFD_R2_exon),sep="")
colnames(LFD_R1_intron) <- paste("LFD_R1_",colnames(LFD_R1_intron),sep="")
colnames(LFD_R2_intron) <- paste("LFD_R2_",colnames(LFD_R2_intron),sep="")
colnames(HFD_R1_intron) <- paste("HFD_R1_",colnames(HFD_R1_intron),sep="")
colnames(HFD_R2_intron) <- paste("HFD_R2_",colnames(HFD_R2_intron),sep="")

# Filter genes to only those present in all libraries and in genes list
Genes <- read.delim("Ensembl_Mouse.txt")
Genes <- Genes[ Genes$Gene.stable.ID %in% rownames(FAP),]
Common <- Reduce(intersect, list(rownames(LFD_R1_exon),rownames(LFD_R2_exon),rownames(HFD_R1_exon),rownames(HFD_R2_exon),rownames(LFD_R1_intron),rownames(LFD_R2_intron),rownames(HFD_R1_intron),rownames(HFD_R2_intron)))
LFD_R1_exon <- LFD_R1_exon[ rownames(LFD_R1_exon) %in% Common & rownames(LFD_R1_exon) %in% Genes[,1],]
LFD_R2_exon <- LFD_R2_exon[ rownames(LFD_R2_exon) %in% Common & rownames(LFD_R2_exon) %in% Genes[,1],]
HFD_R1_exon <- HFD_R1_exon[ rownames(HFD_R1_exon) %in% Common & rownames(HFD_R1_exon) %in% Genes[,1],]
HFD_R2_exon <- HFD_R2_exon[ rownames(HFD_R2_exon) %in% Common & rownames(HFD_R2_exon) %in% Genes[,1],]
LFD_R1_intron <- LFD_R1_intron[ rownames(LFD_R1_intron) %in% Common & rownames(LFD_R1_intron) %in% Genes[,1],]
LFD_R2_intron <- LFD_R2_intron[ rownames(LFD_R2_intron) %in% Common & rownames(LFD_R2_intron) %in% Genes[,1],]
HFD_R1_intron <- HFD_R1_intron[ rownames(HFD_R1_intron) %in% Common & rownames(HFD_R1_intron) %in% Genes[,1],]
HFD_R2_intron <- HFD_R2_intron[ rownames(HFD_R2_intron) %in% Common & rownames(HFD_R2_intron) %in% Genes[,1],]

# Set all objects to the same order of genes
LFD_R2_exon <- LFD_R2_exon[match(rownames(LFD_R1_exon),rownames(LFD_R2_exon)),]
HFD_R1_exon <- HFD_R1_exon[match(rownames(LFD_R1_exon),rownames(HFD_R1_exon)),]
HFD_R2_exon <- HFD_R2_exon[match(rownames(LFD_R1_exon),rownames(HFD_R2_exon)),]
LFD_R1_intron <- LFD_R1_intron[match(rownames(LFD_R1_exon),rownames(LFD_R1_intron)),]
LFD_R2_intron <- LFD_R2_intron[match(rownames(LFD_R1_exon),rownames(LFD_R2_intron)),]
HFD_R1_intron <- HFD_R1_intron[match(rownames(LFD_R1_exon),rownames(HFD_R1_intron)),]
HFD_R2_intron <- HFD_R2_intron[match(rownames(LFD_R1_exon),rownames(HFD_R2_intron)),]

# Generate exonic and intronic read matrix
emat <- cbind(LFD_R1_exon, LFD_R2_exon, HFD_R1_exon, HFD_R2_exon)
imat <- cbind(LFD_R1_intron, LFD_R2_intron, HFD_R1_intron, HFD_R2_intron)

# Order by colnames of FAP
emat <- emat[,match(rownames(FAP@meta.data),colnames(emat))]
imat <- imat[,match(rownames(FAP@meta.data),colnames(imat))]

# Setup the necessary info from Seurat objects
CellColors <- adegenet::fac2col(FAP@active.ident)
names(CellColors) <- names(FAP@active.ident)
CellEmb <- FAP@reductions$umap@cell.embeddings
CellDist <- FAP@reductions$harmony@cell.embeddings
CellDist <- as.dist(1-armaCor(t(CellDist)))

# Filter the genes
emat <- filter.genes.by.cluster.expression(emat, CellColors,0.1)
imat <- filter.genes.by.cluster.expression(imat, CellColors,0.1)

# Generate velocity estimates (smoothend with k=50)
fit.quantile <- 0.05
rval1 <- gene.relative.velocity.estimates(emat,imat,fit.quantile = fit.quantile, cell.dist = CellDist, verbose=T, n.cores=1, kCells=50)

# Plot the velocities on the existing embedding
show.velocity.on.embedding.cor(CellEmb,rval1,n=1000,corr.sigma=0.001,min.grid.cell.mass=20,scale='sqrt',cell.colors=ac(CellColors,alpha=0.4),cex=1,show.grid.flow=T,grid.n=30, n.cores=1)

### Trajectory analysis of FAP1/2 and adipocytes
## Trajectory
# Combine data
Adipogenesis <- merge(FAP, Adipocytes)

# Remove FAP3 and 4
Adipogenesis <- subset(Adipogenesis, subset = Subtype %in% c("LSA","SLSA","LGA","FAP1","FAP2"))

# PHATE embedding  of adipocytes and FAPs.
# NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN. 
# Thus, in order to reproduce downstream analyses using the same Harmony results, please download eWAT_Adipogenesis.Rds
Adipogenesis <- FindVariableFeatures(Adipogenesis, nfeatures=20000)
VariableFeatures(Adipogenesis) <- VariableFeatures(Adipogenesis)[!(VariableFeatures(Adipogenesis) %in% Amb)]
Adipogenesis <- ScaleData(Adipogenesis)
Adipogenesis <- RunPCA(Adipogenesis)
Adipogenesis <- RunHarmony(Adipogenesis, group.by.vars="Dataset")
PHATE <- phate(Embeddings(Tmp, "harmony"), knn = 3, gamma = 0, t = 8)
Adipogenesis[["phate"]] <- CreateDimReducObject(PHATE$embedding, key = "phate")
Adipogenesis@project.name <- "Adipogenesis"

# DimPlots
DimPlot(Adipogenesis, reduction = "phate")
DimPlot(Adipogenesis, reduction = "phate", split.by="Diet")
DimPlot(Adipogenesis, reduction = "phate", split.by="Replicate")

# Extract phate coordinates and train a tree
tree_data <- Embeddings(Adipogenesis, "phate")
TreeEPG <- computeElasticPrincipalTree(X = tree_data, NumNodes = 30, Lambda = .005, Mu = .001, drawAccuracyComplexity = FALSE, drawEnergy = FALSE)

# Plot the result
PlotPG(X = tree_data, TargetPG = TreeEPG[[1]], DimToPlot = 1:2, Do_PCA = FALSE)

# Infer pseudotime
Tree_Graph <- ConstructGraph(TreeEPG[[1]])
Tree_e2e <- GetSubGraph(Net = Tree_Graph, Structure = 'branches')
PartStruct <- PartitionData(X = tree_data, NodePositions = TreeEPG[[1]]$NodePositions)
ProjStruct <- project_point_onto_graph(X = tree_data,
                                       NodePositions = TreeEPG[[1]]$NodePositions,
                                       Edges = TreeEPG[[1]]$Edges$Edges,
                                       Partition = PartStruct$Partition)

# Get only pseudotime for branch 3 (the differentiation branch connecting FAPs and adipocytes)
Pt <- getPseudotime(ProjStruct = ProjStruct, NodeSeq = rev(Tree_e2e$Branch_3))
Adipogenesis$Pt <- Pt$Pt

# Save the object
saveRDS(Adipogenesis, "eWAT_Adipogenesis.Rds")

## Meta-stable states
# Subset to only the differentiation branch
Differentiation <- subset(Adipogenesis, Pt != "NA")

# Define the intermediary state
PHATE <- Embeddings(Differentiation, "phate")
Differentiation@meta.data$Type = "NA"
Differentiation@meta.data[ rownames(Differentiation@meta.data) %in% rownames(PHATE[ PHATE[,1] >= 0.02 & PHATE[,2] > 0,]), "Type"] <- "Early_Pre"
Differentiation@meta.data[ rownames(Differentiation@meta.data) %in% rownames(PHATE[ PHATE[,1] >= 0.0235 & PHATE[,2] < 0,]), "Type"] <- "Late_Pre"
Differentiation@meta.data[ rownames(Differentiation@meta.data) %in% rownames(PHATE[ PHATE[,1] < 0.0235 & PHATE[,1] > -0.0075 & PHATE[,2] < 0,]), "Type"] <- "Transitioning"
Differentiation@meta.data[ rownames(Differentiation@meta.data) %in% rownames(PHATE[ PHATE[,1] < -0.0075,]), "Type"] <- "Mature"

# DimPlots
DimPlot(Differentiation, reduction = "phate")
DimPlot(Differentiation, reduction = "phate", split.by="Diet")
DimPlot(Differentiation, reduction = "phate", split.by="Replicate")

# Histograms
H_R2 <- hist(Differentiation@meta.data[ grep("R2",Differentiation@meta.data$Dataset), "Pt"], breaks=seq(0,0.08, by=0.005))
H_R1 <- hist(Differentiation@meta.data[ grep("R1",Differentiation@meta.data$Dataset), "Pt"], breaks=seq(0,0.08, by=0.005))
Blt <- barplot(rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts))), ylim=c(0, 0.6), las = 1)
arrows(Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))+(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))-(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),code = 3, angle = 90, len=0.05)

## Gene modules
Overall_Markers <- readRDS("eWAT_Overall_Markers.Rds") # Generated during analysis of the main datasets

# Extract marker gene results
Clust1 <- Overall_Markers[[2]] # FAPs
Clust2 <- Overall_Markers[[3]] # Adipocytes

# Subset to genes that are enriched or exclusive (!= NS)
Clust1 <- Clust1[Clust1$Marker != "NS", ]
Clust2 <- Clust2[Clust2$Marker != "NS", ]

# Extracting the minimum logFC across replicates in pairwise cluster comparisons
# RATIONALE: The genes with the highest minimum fold change are the most specific ones for any given cluster
Clust1$logFC_OvO <- apply(Clust1[, grep("logFC_Cluster", colnames(Clust1))], 1, FUN = "min")
Clust2$logFC_OvO <- apply(Clust2[, grep("logFC_Cluster", colnames(Clust2))], 1, FUN = "min")

# Computing gene module scores using the top 50 most specific marker genes
Differentiation <- AddModuleScore(Differentiation, features = list(Cluster1 = Clust1[order(-Clust1$logFC_OvO), ][1:50, ]$Symbol,
                                                         Cluster2 = Clust2[order(-Clust2$logFC_OvO), ][1:50, ]$Symbol))

# Scale the module scores (across all nuclei)
Differentiation$Cluster1 <- scale(Differentiation$Cluster1)
Differentiation$Cluster2 <- scale(Differentiation$Cluster2)

# Boxplots
par(mfcol=c(1,2))
boxplot(
  Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre", "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Type == "Late_Pre", "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Type == "Transitioning", "Cluster1"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Type == "Mature", "Cluster1"][,1],
  outline=F, las=1, main = "FAP signature", names = c("Early_Pre", "Late_Pre","Transitioning","Mature"), ylab="Scaled score")
boxplot(
  Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre", "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Type == "Late_Pre", "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Type == "Transitioning", "Cluster2"][,1],
  Differentiation@meta.data[ Differentiation@meta.data$Type == "Mature", "Cluster2"][,1],
  outline=F, las=1, main = "Adipocyte signature", names = c("Early_Pre", "Late_Pre","Transitioning","Mature"), ylab="Scaled score")

### Trajectory validation using additional tools
## Trajectory
# Set labels (adipocytes set to a single label to focus on differentiation)
Differentiation$Label <- "NA"
Differentiation@meta.data[ Differentiation@meta.data$Annotation == "Adipocyte" & Differentiation@meta.data$seurat_clusters == "0", "Label"] <- 1
Differentiation@meta.data[ Differentiation@meta.data$Annotation == "Adipocyte" & Differentiation@meta.data$seurat_clusters == "1", "Label"] <- 1
Differentiation@meta.data[ Differentiation@meta.data$Annotation == "Adipocyte" & Differentiation@meta.data$seurat_clusters == "2", "Label"] <- 1
Differentiation@meta.data[ Differentiation@meta.data$Annotation == "FAP" & Differentiation@meta.data$seurat_clusters == "0", "Label"] <- 2
Differentiation@meta.data[ Differentiation@meta.data$Annotation == "FAP" & Differentiation@meta.data$seurat_clusters == "3", "Label"] <- 3

# Extract pseudo-time from the original object
ElPiTime <- as.data.frame(Differentiation@meta.data[,c("Label","Pt")])
ElPiTime$Cell <- rownames(ElPiTime)
colnames(ElPiTime)[2] <- "ElPiGraph"

## Test with different tools
# Slingshot
Test <- slingshot(Embeddings(Differentiation, "phate"), clusterLabels=Differentiation$Label)
SlingTime <- as.data.frame(slingPseudotime(Test))
SlingTime$Cell <- rownames(SlingTime)
colnames(SlingTime)[1] <- "Slingshot"

# TSCAN
Test <- exprmclust(as.data.frame(t(Embeddings(Differentiation, "phate"))), reduce=F,clusternum=2:50)
TSCANtime <- TSCANorder(Test, orderonly=F)
colnames(TSCANtime)[c(1,3)] <- c("Cell","TSCAN")

# DPT
DM <- DiffusionMap(Embeddings(Differentiation, "harmony"))
dpt <- DPT(DM)
DESTINY <- as.data.frame(DPT$DPT1942) # First cell in same trajectory as ElPiGraph.R
DESTINY$Cell <- rownames(Differentiation@meta.data)
colnames(DESTINY)[1] <- "Destiny"

# SCORPIUS
space <- reduce_dimensionality(Embeddings(Differentiation, "harmony"), dist = "pearson", ndim = 2)
traj <- infer_trajectory(space)
SCORPIUS <- as.data.frame(traj$time)
SCORPIUS$Cell <- rownames(SCORPIUS)
colnames(SCORPIUS)[1] <- "SCORPIUS"

# Monocle3
MD <- data.frame(gene_short_name = rownames(Differentiation))
rownames(MD) <- rownames(Differentiation)
cds <- new_cell_data_set(Differentiation@assays$RNA@counts, cell_metadata = Differentiation@meta.data, gene_metadata = MD)
cds <- preprocess_cds(cds, num_dim = 50, verbose=T, use_genes = VariableFeatures(Differentiation))
reducedDims(cds)[["Aligned"]] <- as.matrix(Embeddings(Differentiation, "harmony"))
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds@clusters[["UMAP"]]$clusters <- as.factor(Differentiation$Label)
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = names(which.min(Differentiation$Pt)))
Mcl3 <- as.data.frame(cds@principal_graph_aux[["UMAP"]]$pseudotime)
Mcl3$Cell <- rownames(Mcl3)
colnames(Mcl3)[1] <- "Monocle3"

## Combine results
Pt <- merge(ElPiTime, SlingTime, by="Cell")
Pt <- merge(Pt, TSCANtime[,c(1,3)], by="Cell")
Pt <- merge(Pt, SCORPIUS, by="Cell")
Pt <- merge(Pt, DESTINY, by="Cell")
Pt <- merge(Pt, Mcl3, by="Cell")
Pt <- Pt[ !is.na(Pt$ElPiGraph),]

# Calculate and plot the correlation
barplot(abs(cor(Pt[,3:8], method="spearman")[c(2,3,5,4,6),1]), ylim=c(0,1), las=1)

## Meta-stable states
# Scorpius
H_R2 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R1",Differentiation@meta.data$Dataset),]),"SCORPIUS"], breaks=seq(0,1.0,by=0.05))
H_R1 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R2",Differentiation@meta.data$Dataset),]),"SCORPIUS"], breaks=seq(0,1.0,by=0.05))
Blt <- barplot(rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts))), ylim=c(0, 0.2))
arrows(Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))+(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))-(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),code = 3, angle = 90, len=0.05)

# Destiny
H_R2 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R1",Differentiation@meta.data$Dataset),]),"Destiny"], breaks=seq(0,0.5,by=0.025))
H_R1 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R2",Differentiation@meta.data$Dataset),]),"Destiny"], breaks=seq(0,0.5,by=0.025))
Blt <- barplot(rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts))), ylim=c(0, 0.5))
arrows(Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))+(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))-(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),code = 3, angle = 90, len=0.05)

# Monocle3
H_R2 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R1",Differentiation@meta.data$Dataset),]),"Monocle3"], breaks=seq(0,30,by=0.75))
H_R1 <- hist(Pt[ Pt$Cell %in% rownames(Differentiation@meta.data[ grep("R2",Differentiation@meta.data$Dataset),]),"Monocle3"], breaks=seq(0,30,by=0.75))
Blt <- barplot(rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts))), ylim=c(0, 0.2))
arrows(Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))+(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),Blt, rowMeans(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)))-(apply(cbind(H_R1$counts/sum(H_R1$counts),H_R2$counts/sum(H_R2$counts)),1,FUN="sd")/sqrt(2)),code = 3, angle = 90, len=0.05)

### Differential expression along trajectory (with tradeSeq)
## Find expressed genes (in at least 10% of the nuclei)
sca <-  FromMatrix(as.matrix(Differentiation@assays$RNA@data), Differentiation@meta.data)
expressed_genes <- names(which(freq(sca) >= 0.1))

### Association tests
### NOTE: This step is non-deterministic! Results vary from run-to-run
## LFD
# Setup
Differentiation_LFD <- subset(Differentiation, subset = Diet == "LFD")
Differentiation_R1 <- subset(Differentiation_LFD, subset = Replicate == "R1")
Differentiation_R2 <- subset(Differentiation_LFD, subset = Replicate == "R2")

# Replicate 1
Pseudotime <- data.frame(curve1 = Differentiation_R1$Pt)
Weights <- data.frame(LFD = rep(1, ncol(Differentiation_R1)))
rownames(Weights) = colnames(Differentiation_R1)
Counts <- Differentiation_R1@assays$RNA@counts
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
LFD_sce_R1 <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)

# Replicate 2
Pseudotime <- data.frame(curve1 = Differentiation_R2$Pt)
Weights <- data.frame(LFD = rep(1, ncol(Differentiation_R2)))
rownames(Weights) = colnames(Differentiation_R2)
Counts <- Differentiation_R2@assays$RNA@counts
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
LFD_sce_R2 <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)

## HFD
# Setup
Differentiation_HFD <- subset(Differentiation, subset = Diet == "HFD")
Differentiation_R1 <- subset(Differentiation_HFD, subset = Replicate == "R1")
Differentiation_R2 <- subset(Differentiation_HFD, subset = Replicate == "R2")

# Replicate 1
Pseudotime <- data.frame(curve1 = Differentiation_R1$Pt)
Weights <- data.frame(LFD = rep(1, ncol(Differentiation_R1)))
rownames(Weights) = colnames(Differentiation_R1)
Counts <- Differentiation_R1@assays$RNA@counts
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
HFD_sce_R1 <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)

# Replicate 2
Pseudotime <- data.frame(curve1 = Differentiation_R2$Pt)
Weights <- data.frame(LFD = rep(1, ncol(Differentiation_R2)))
rownames(Weights) = colnames(Differentiation_R2)
Counts <- Differentiation_R2@assays$RNA@counts
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
HFD_sce_R2 <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)

## Tests and filtering
# Perform tests
Asso_LFD_R1 <- as.data.frame(associationTest(LFD_sce_R1))
Asso_LFD_R2 <- as.data.frame(associationTest(LFD_sce_R2))
Asso_HFD_R1 <- as.data.frame(associationTest(HFD_sce_R1))
Asso_HFD_R2 <- as.data.frame(associationTest(HFD_sce_R2))

# FDR-adjustment
Asso_LFD_R1$FDR <- p.adjust(Asso_LFD_R1$pvalue, method="fdr")
Asso_LFD_R2$FDR <- p.adjust(Asso_LFD_R2$pvalue, method="fdr")
Asso_HFD_R1$FDR <- p.adjust(Asso_HFD_R1$pvalue, method="fdr")
Asso_HFD_R2$FDR <- p.adjust(Asso_HFD_R2$pvalue, method="fdr")

# Filtering
Asso_LFD <- rownames(Asso_LFD_R1[ Asso_LFD_R1$FDR <= 0.05 & rownames(Asso_LFD_R1) %in% rownames(Asso_LFD_R2[Asso_LFD_R2$FDR <= 0.05,]),])
Asso_HFD <- rownames(Asso_HFD_R1[ Asso_HFD_R1$FDR <= 0.05 & rownames(Asso_HFD_R1) %in% rownames(Asso_HFD_R2[Asso_HFD_R2$FDR <= 0.05,]),])

# Combine
Ass <- unique(sort(c(Asso_LFD, Asso_HFD)))

### Heatmap of genes associated with pseudo-time in LFD
## Predict smoothend expression 
Smooth_R1 <- predictSmooth(LFD_sce_R1, gene = Asso_LFD, tidy=F, n=100)
Smooth_R2 <- predictSmooth(LFD_sce_R2, gene = Asso_LFD, tidy=F, n=100)

# Average across replicates and scale
Smooth <- Smooth_R1
for (i in 1:nrow(Smooth)) { Smooth[i,] <- colMeans(rbind(Smooth_R1[i,], Smooth_R2[i,])) }
Smooth <- t(scale(t(Smooth)))

# Seriate the results
Smooth <- Smooth[ get_order(seriate(Smooth, method="PCA_angle")),]

# Create heatmaps (in 4 visually defined groups)
HTM1 <- Heatmap(Smooth[1:210,], cluster_columns=F, cluster_rows=F)
HTM2 <- Heatmap(Smooth[211:320,], cluster_columns=F, cluster_rows=F)
HTM3 <- Heatmap(Smooth[321:800,], cluster_columns=F, cluster_rows=F)
HTM4 <- Heatmap(Smooth[801:1255,], cluster_columns=F, cluster_rows=F)

# Plot the heatmap
HTM3 %v% HTM4 %v% HTM1 %v% HTM2

## Pathway analysis
# Extract genes for each group
C1 <- rownames(Smooth[1:210,])
C2 <- rownames(Smooth[211:320,])
C3 <- rownames(Smooth[321:800,])
C4 <- rownames(Smooth[801:1255,])

# Setup
gmt <- rWikiPathways::downloadPathwayArchive(organism = "Mus musculus", format = "gmt")
wp2gene <- clusterProfiler::read.gmt(gmt)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name", "version", "wpid", "org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) # TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) # TERM2NAME

# Convert gene names
C1_Entrez <- clusterProfiler::bitr(C1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
C2_Entrez <- clusterProfiler::bitr(C2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
C3_Entrez <- clusterProfiler::bitr(C3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
C4_Entrez <- clusterProfiler::bitr(C4, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Do the pathway analysis
wiki_clust1 <- clusterProfiler::enricher(C1_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1,  TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
wiki_clust2 <- clusterProfiler::enricher(C2_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
wiki_clust3 <- clusterProfiler::enricher(C3_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
wiki_clust4 <- clusterProfiler::enricher(C4_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Adding gene symbols to the resulting pathway file
wiki_clust1 <- as.data.frame(DOSE::setReadable(wiki_clust1, org.Mm.eg.db, keyType = "ENTREZID"))
wiki_clust2 <- as.data.frame(DOSE::setReadable(wiki_clust2, org.Mm.eg.db, keyType = "ENTREZID"))
wiki_clust3 <- as.data.frame(DOSE::setReadable(wiki_clust3, org.Mm.eg.db, keyType = "ENTREZID"))
wiki_clust4 <- as.data.frame(DOSE::setReadable(wiki_clust4, org.Mm.eg.db, keyType = "ENTREZID"))

## Barplot of top 3 pathways
# Cluster 1: Fasn, Scd1, Acaca, Pparg, Fabp4, Lep, Lipe, Plin1, Cd36
pdf("Wikipathway_Cluster1.pdf", width=10, height=10, useDingbats=F)
barplot(-log10(wiki_clust1[1:3,"pvalue"]), names=wiki_clust1[1:3,"Description"], las=2)
dev.off()

# Cluster 2: Foxo1, Pten, Igf1, Akt2, Gsk3b, Ehd2
pdf("Wikipathway_Cluster2.pdf", width=10, height=10, useDingbats=F)
barplot(-log10(wiki_clust2[1:3,"pvalue"]), names=wiki_clust2[1:3,"Description"], las=2)
dev.off()

# Cluster 3: Atp5a1, Rpl37, Rps2, Cov5a, Ndufa4,  Ubcrb
pdf("Wikipathway_Cluster3.pdf", width=10, height=10, useDingbats=F)
barplot(-log10(wiki_clust3[1:3,"pvalue"]), names=wiki_clust3[1:3,"Description"], las=2)
dev.off()

# Cluster 4: Egfr, Pdgfrb, Pdgfra, Eps8, Eps15, Spry2
pdf("Wikipathway_Cluster4.pdf", width=10, height=10, useDingbats=F)
barplot(-log10(wiki_clust1[1:3,"pvalue"]), names=wiki_clust1[1:3,"Description"], las=2)
dev.off()

### Transcription factor waves
# Import list of TFs in mice
TF <- read.delim("Mus_musculus_TF.txt", header=T)

# Subset smoothend expression of genes associated with pseudo-time in LFD to only TFs
SmoothTF <- Smooth[ rownames(Smooth) %in% TF[,2],]

# Seriate the result
SmoothTF <- SmoothTF[ get_order(seriate(SmoothTF, method="PCA_angle")),]

# Create heatmaps (manual reordering to make it look nice)
HTM1 <- Heatmap(SmoothTF[c(1,3,4,9,6,7,2,11,5,10,8,12,13),], cluster_columns=F, cluster_rows=F)
HTM2 <- Heatmap(SmoothTF[c(15,14,18,17,23,16,22,19,21,20,24),], cluster_columns=F, cluster_rows=F)
HTM3 <- Heatmap(SmoothTF[c(25,29,30,31,32,33,27,34,28,26,35,36,37,38,39,40),], cluster_columns=F, cluster_rows=F)
HTM4 <- Heatmap(SmoothTF[c(41,42,48,72,52,43,44,46,45,47,50,61,49,69,51,55,53,54,60,59,58,65,62,68,56,71,64,57,74,75,63,70,66,67,73,76,77),], cluster_columns=F, cluster_rows=F)

# Plot the heatmap
HTM3 %v% HTM4 %v% HTM1 %v% HTM2

### Differential association with pseudo-time in LFD and HFD
## Pattern test
# Setup: R1
Differentiation_R1 <- subset(Differentiation, subset = Replicate == "R1")
Pseudotime <- data.frame(curve1 = Differentiation_R1$Pt, curve2 = Differentiation_R1$Pt)
Weights <- data.frame(LFD = rep(0, ncol(Differentiation_R1)),HFD = rep(0, ncol(Differentiation_R1)))
rownames(Weights) = colnames(Differentiation_R1)
Weights[ rownames(Weights) %in% rownames(Differentiation_R1@meta.data[Differentiation_R1@meta.data$Diet == "LFD",]),"LFD"] <- 1
Weights[ rownames(Weights) %in% rownames(Differentiation_R1@meta.data[Differentiation_R1@meta.data$Diet == "HFD",]),"HFD"] <- 1
Counts <- Differentiation_R1@assays$RNA@counts
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
R1_sce <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)

# Setup: R2
Differentiation_R1 <- subset(Differentiation, subset = Replicate == "R2")
Pseudotime <- data.frame(curve1 = Differentiation_R1$Pt, curve2 = Differentiation_R1$Pt)
Weights <- data.frame(LFD = rep(0, ncol(Differentiation_R1)),HFD = rep(0, ncol(Differentiation_R1)))
rownames(Weights) = colnames(Differentiation_R1)
Weights[ rownames(Weights) %in% rownames(Differentiation_R1@meta.data[Differentiation_R1@meta.data$Diet == "LFD",]),"LFD"] <- 1
Weights[ rownames(Weights) %in% rownames(Differentiation_R1@meta.data[Differentiation_R1@meta.data$Diet == "HFD",]),"HFD"] <- 1
Counts <- Differentiation_R1@assays$RNA@counts
Counts <- Counts[ rownames(Counts) %in% expressed_genes,]
R2_sce <- fitGAM(counts = Counts, pseudotime = Pseudotime, cellWeights = Weights, nknots = 4, verbose = TRUE)

## Tests and filtering
# Perform tests
Pat_R1 <- as.data.frame(patternTest(R1_sce))
Pat_R2 <- as.data.frame(patternTest(R2_sce))

# Subset to genes associated with pseudo-time
Pat_R2 <- Pat_R2[ rownames(Pat_R2) %in% Ass,]
Pat_R1 <- Pat_R1[ rownames(Pat_R1) %in% Ass,]

# FDR-adjustment
Pat_R1$FDR <- p.adjust(Pat_R1$pvalue, method="fdr")
Pat_R2$FDR <- p.adjust(Pat_R2$pvalue, method="fdr")

# Subset to significant ones
Pat_R1 <- Pat_R1[ Pat_R1$FDR <= 0.05,]
Pat_R2 <- Pat_R2[ Pat_R2$FDR <= 0.05,]

# Remove NAs
Pat_R1 <- Pat_R1[ grep("NA", rownames(Pat_R1), invert=T),]
Pat_R2 <- Pat_R2[ grep("NA", rownames(Pat_R2), invert=T),]

# Combine the results
Pat <- rownames(Pat_R1[ Pat_R1$FDR <= 0.05 & rownames(Pat_R1) %in% rownames(Pat_R2[ Pat_R2$FDR <= 0.05,]),])

### Heatmap of differentially associated TFs
# Import list of mouse TFs
TF <- read.delim("Mus_musculus_TF.txt", header=T)

# Subset results to only TFs
Pat_TFs <- Pat[ Pat %in% TF[,2]]

# Predict smoothend expression
Smooth_R1 <- predictSmooth(R1_sce, gene = Pat_TFs, tidy=F, n=100)
Smooth_R2 <- predictSmooth(R2_sce, gene = Pat_TFs, tidy=F, n=100)

# Average smoothend expression across replicates and scale
Smooth <- Smooth_R1
for (i in 1:nrow(Smooth)) { Smooth[i,] <- colMeans(rbind(Smooth_R1[i,], Smooth_R2[i,])) }
Smooth <- t(scale(t(Smooth)))

# Reorder the data (manually to make it look good visually)
Smooth <- Smooth[ c(5,17,25,8,3,9,13,6,28,10,24,19,2,4,27,7,18,23,20,1,14,11,22,12,16,15,29,26,21),]

# Plot the heatmap
Heatmap( Smooth, cluster_columns=F, cluster_rows=F, right_annotation = rowAnnotation(bar = c("U","I","U","U","I","U","I","U","I","U","U","I","A","I","U","I","U","A","U","A","A","A","I","U","U","U","I","U","U")))

### Cell-cell signaling
## Process data from CellPhoneDB
# Import gene lists and annotations (downloaded from here https://www.cellphonedb.org/downloads)
Genes_Annotations <- read.delim("gene_input.txt", header=T, sep=",")
Complex_Annotation <- read.delim("complex_input.txt", header=T, sep=",")
Interactions <- read.delim("interaction_curated.txt", sep = ",", header=T)
Complexes <- read.delim("complex_curated.txt", sep = ",")
Proteins <- read.delim("protein_curated.txt", sep = ",")

# Setup to loop
Pairs <- as.data.frame(matrix(ncol=7, nrow=1))
colnames(Pairs) <- c("ID","Receptor","Ligand","Receptor_Genes","Ligand_Genes", "Type_Receptor", "Type_Ligand")

# Loop through all interactions and save information about the partners and their interactions
Counter <- 1
for (i in 1:nrow(Interactions)) {
# Check if partner a is a secreted factor
if (nrow(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),]) == 1 | nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$secreted == "True",]) == 1) {
  # Check if partner b is a receptor
  if (nrow(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),]) == 1 | nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$receptor == "True",]) == 1) {
    # Save the interaction information
    Pairs[Counter,1] <- as.character(Interactions[i,1])
    if (nrow(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),]) == 1) {
      Pairs[Counter,"Ligand"] <- as.character(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),1])
      Pairs[Counter,"Ligand_Genes"] <- paste0(Genes_Annotations[ Genes_Annotations$uniprot %in% as.character(unlist(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),c(2,3,4,5)])),3], collapse = ",")
      Pairs[Counter,"Type_Ligand"] <- "Complex"
    }
    if (nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$secreted == "True",]) == 1) {
      Pairs[Counter,"Ligand"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$secreted == "True",1],3]), collapse = ",")
      Pairs[Counter,"Ligand_Genes"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$secreted == "True",1],3]), collapse = ",")
      Pairs[Counter,"Type_Ligand"] <- "Interaction"
    }
    
    if (nrow(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),]) == 1) {
      Pairs[Counter,"Receptor"] <- as.character(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),1])
      Pairs[Counter,"Receptor_Genes"] <- paste0(Genes_Annotations[ Genes_Annotations$uniprot %in% as.character(unlist(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),c(2,3,4,5)])),3], collapse = ",")
      Pairs[Counter,"Type_Receptor"] <- "Complex"
    }
    if (nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$receptor == "True",]) == 1) {
      Pairs[Counter,"Receptor"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$receptor == "True",1],3]), collapse = ",")
      Pairs[Counter,"Receptor_Genes"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$receptor == "True",1],3]), collapse = ",")
      Pairs[Counter,"Type_Receptor"] <- "Interaction"
    }
    Counter <- Counter + 1
  }
}
# Check if partner b is a secreted factor
if (nrow(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),]) == 1 | nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$secreted == "True",]) == 1) {
  # Check if partner b is a receptor
  if (nrow(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),]) == 1 | nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$receptor == "True",]) == 1) {
    # Save the interaction information
    Pairs[Counter,1] <- as.character(Interactions[i,1])
    if (nrow(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),]) == 1) {
      Pairs[Counter,"Ligand"] <- as.character(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),1])
      Pairs[Counter,"Ligand_Genes"] <- paste0(Genes_Annotations[ Genes_Annotations$uniprot %in% as.character(unlist(Complexes[ Complexes$secreted == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_b"]),c(2,3,4,5)])),3], collapse = ",")
      Pairs[Counter,"Type_Ligand"] <- "Complex"
    }
    if (nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$secreted == "True",]) == 1) {
      Pairs[Counter,"Ligand"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$secreted == "True",1],3]), collapse = ",")
      Pairs[Counter,"Ligand_Genes"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_b"]) & Proteins$secreted == "True",1],3]), collapse = ",")
      Pairs[Counter,"Type_Ligand"] <- "Interaction"
    }
    
    if (nrow(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),]) == 1) {
      Pairs[Counter,"Receptor"] <- as.character(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),1])
      Pairs[Counter,"Receptor_Genes"] <- paste0(Genes_Annotations[ Genes_Annotations$uniprot %in% as.character(unlist(Complexes[ Complexes$receptor == "TRUE" & Complexes[,1] == as.character(Interactions[i, "partner_a"]),c(2,3,4,5)])),3], collapse = ",")
      Pairs[Counter,"Type_Receptor"] <- "Complex"
    }
    if (nrow(Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$receptor == "True",]) == 1) {
      Pairs[Counter,"Receptor"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$receptor == "True",1],3]), collapse = ",")
      Pairs[Counter,"Receptor_Genes"] <- paste0(as.character(Genes_Annotations[ Genes_Annotations$uniprot %in% Proteins[ Proteins$uniprot == as.character(Interactions[i, "partner_a"]) & Proteins$receptor == "True",1],3]), collapse = ",")
      Pairs[Counter,"Type_Receptor"] <- "Interaction"
    }
    Counter <- Counter + 1
  }
}
}

## Import gene mappings between human and mouse (from Ensembl biomart)
Convert <- read.delim("Human_Mouse.txt", header=T)

## First analyze expression of the receptors in early pre-adipocytes
# Subset
Receptors <- Pairs[ duplicated(Pairs$V2) == F,]

# Split into wether its a complex or single receptor.
Complexes <- Receptors[ Receptors$V6 == "Complex",]
Interactions <- Receptors[ Receptors$V6 != "Complex",]

# Process complexes. For each complex, find the lowest expressed member in each diet.
for (i in 1:nrow(Complexes)) {
	Tmp <- Convert[ Convert[,1] %in% unlist(strsplit(as.character(Complexes[i,4]),",")),]
	for (q in 1:nrow(Tmp)) {
		if (sum(rownames(Differentiation@assays$RNA@data) %in% as.character(Tmp[q,2])) == 1) {
			Tmp[q,3] <- mean(Differentiation@assays$RNA@data[ rownames(Differentiation@assays$RNA@data) %in% as.character(Tmp[q,2]), colnames(Differentiation@assays$RNA@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "LFD",])])
			Tmp[q,4] <- mean(Differentiation@assays$RNA@data[ rownames(Differentiation@assays$RNA@data) %in% as.character(Tmp[q,2]), colnames(Differentiation@assays$RNA@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "HFD",])])
		} else {
			Tmp[q,3] <- 0
			Tmp[q,4] <- 0
		}
	}
	# LFD
	Tmp <- Tmp[ order(Tmp[,1], -Tmp[,3]),]
	Tmp2 <- Tmp[ duplicated(Tmp$Gene.name) == F,]
	Complexes[i,8] <- min(Tmp2[,3])
	# HFD
	Tmp <- Tmp[ order(Tmp[,1], -Tmp[,4]),]
	Tmp2 <- Tmp[ duplicated(Tmp$Gene.name) == F,]
	Complexes[i,9] <- min(Tmp2[,4])
}

# Process interactions.
for (i in 1:nrow(Interactions)) {
	if (sum(rownames(Differentiation@assays$RNA@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2])) == 0) {
	Interactions[i,8] <- 0
	Interactions[i,9] <- 0
	}
	if (sum(rownames(Differentiation@assays$RNA@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2])) == 1) {
	Interactions[i,8] <- mean(Differentiation@assays$RNA@data[ rownames(Differentiation@assays$RNA@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2]), colnames(Differentiation@assays$RNA@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "LFD",])])
	Interactions[i,9] <- mean(Differentiation@assays$RNA@data[ rownames(Differentiation@assays$RNA@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2]), colnames(Differentiation@assays$RNA@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "HFD",])])
	}
	if (sum(rownames(Differentiation@assays$RNA@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2])) > 1) {
	Interactions[i,8] <- Interactions[i,8] <- max(Matrix::rowMeans(Differentiation@assays$RNA@data[ rownames(Differentiation@assays$RNA@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2]), colnames(Differentiation@assays$RNA@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "LFD",])]))
	Interactions[i,9] <- Interactions[i,8] <- max(Matrix::rowMeans(Differentiation@assays$RNA@data[ rownames(Differentiation@assays$RNA@data) %in% as.character(Convert[ Convert[,1] %in% unlist(strsplit(as.character(Interactions[i,4]),",")),2]), colnames(Differentiation@assays$RNA@data) %in% rownames(Differentiation@meta.data[ Differentiation@meta.data$Type == "Early_Pre" & Differentiation$Diet == "HFD",])]))
	}
}

# Combine the information
Receptors <- rbind(Complexes, Interactions)

# Subset to only expressed receptors (higher than 0.1 in each LFD or HFD)
Receptors <- Receptors[ Receptors[,9] >= 0.1 | Receptors[,8] >= 0.1,] # 32 receptors

## Now analyze expression of the ligands in cell types endogenous to the adipose tissue
# Subset to ligands of expressed receptors
Ligands <- Pairs[ Pairs[,2] %in% Receptors[,2],]
Ligands <- Ligands[ duplicated(Ligands[,3]) == F,]
Ligands[,5] <- as.character(Ligands[,5])

# Loop through all ligands
AllLigands <- data.frame()
for (i in 2:nrow(Ligands)) {
	Tmp <- Convert[ Convert[,1] %in% unlist(strsplit(Ligands[i,5],",")),]
	Tmp$Adipocytes_LFD <- 0
	Tmp$Adipocytes_HFD <- 0
	Tmp$Mesothelial_LFD <- 0
	Tmp$Mesothelial_HFD <- 0
	Tmp$Endothelial_LFD <- 0
	Tmp$Endothelial_HFD <- 0
	Tmp$FAP_LFD <- 0
	Tmp$FAP_HFD <- 0
	Tmp$Immune_LFD <- 0
	Tmp$Immune_HFD <- 0
	# Loop over all genes homologous to the ligands and get their mean expression
	for (q in 1:nrow(Tmp)) {
		if (sum(rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2])) > 0) {
			Tmp[q,"Adipocytes_LFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Adipocyte" & Data$Diet == "LFD",])])
			Tmp[q,"Adipocytes_HFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Adipocyte" & Data$Diet == "HFD",])])
			Tmp[q,"Mesothelial_LFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Mesothelial" & Data$Diet == "LFD",])])
			Tmp[q,"Mesothelial_HFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Mesothelial" & Data$Diet == "HFD",])])
			Tmp[q,"Endothelial_LFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Endothelial" & Data$Diet == "LFD",])])
			Tmp[q,"Endothelial_HFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Endothelial" & Data$Diet == "HFD",])])
			Tmp[q,"FAP_LFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "FAP" & Data$Diet == "LFD",])])
			Tmp[q,"FAP_HFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "FAP" & Data$Diet == "HFD",])])
			Tmp[q,"Immune_LFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Immune" & Data$Diet == "LFD",])])
			Tmp[q,"Immune_HFD"] <- mean(eWAT@assays$RNA@counts[ rownames(eWAT@assays$RNA@counts) %in% as.character(Tmp[q,2]), colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation == "Immune" & Data$Diet == "HFD",])])
		}
	}
	AllLigands <- rbind(AllLigands, Tmp)
}

# Remove ligands without names
AllLigands <- AllLigands[ AllLigands[,2] != "",]

# Subset to expressed ligands (>= 0.1) across cell types and diets and merge results
AllLigands <- AllLigands[which(apply(AllLigands[,3:ncol(AllLigands)],1,FUN="max") >= 0.1),]
Ligands <- merge(Ligands, AllLigands[,c(1,3:ncol(AllLigands))], by.x=3, by.y=1) # 17 ligands

## Combine the result of the ligand and receptor analysis
Interactions <- Pairs[ Pairs[,2] %in% Receptors[,2] & Pairs[,3] %in% Ligands[,1],1:3]
colnames(Receptors)[8] <- "Receptor_LFD"
colnames(Receptors)[9] <- "Receptor_HFD"
Interactions <- merge(Interactions, Receptors[,c(2,8,9)], by.x=2, by.y=1)
Interactions <- merge(Interactions, Ligands[,c(1,8:17)], by.x=3, by.y=1) # 18 receptors, 17 ligands, 29 interactions

## Heatmap the expression of receptors in early pre-adipocytes
Receptors <- Interactions[ duplicated(Interactions[,2])==F,]
rownames(Receptors) <- Receptors[,2]
ComplexHeatmap::Heatmap(Receptors[c(5,6,3,4,1,7,8,9,10,11,15,12,13,14,2,17,16,18),4:5], cluster_rows=F, cluster_columns=F, col = colorRamp2(breaks = seq(0,1.65, length.out=11), rev(brewer.pal(11,"Spectral"))))

## Heatmap of the expression of ligands in cell types endogenous to the adipose tissue
Ligands <- Interactions[ duplicated(Interactions[,1])==F,]
rownames(Ligands) <- Ligands[,1]
ComplexHeatmap::Heatmap(Ligands[c(4,3,1,5,10,6,7,8,9,14,17,11,13,12,2,15,16),6:15], cluster_rows=F, cluster_columns=F, col = colorRamp2(breaks = seq(0,1.65, length.out=11), rev(brewer.pal(11,"Spectral"))))
```
[Back to start](../README.md)<br>




