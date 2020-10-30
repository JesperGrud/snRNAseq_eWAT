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
eWAT <- readRDS("eWAT_Annotated.Rds")
Amb <- readRDS("Ambient.Rds")
Endomeso <- subset(eWAT, subset = Annotation %in% c("Mesothelial", "Endothelial"))

# Embedding and clustering of mesothelial and endothelial cells into subpopulations
# NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN. 
# Thus, in order to reproduce downstream analyses using the same Harmony results, please download eWAT_Endothelial_Mesothelial.Rds
Endomeso <- FindVariableFeatures(Endomeso, nfeatures=3000)
VariableFeatures(Endomeso) <- VariableFeatures(Endomeso)[!(VariableFeatures(Endomeso) %in% Amb)]
Endomeso <- ScaleData(Endomeso)
Endomeso <- RunPCA(Endomeso)
Endomeso <- RunHarmony(Endomeso, group.by.vars="Dataset")
Endomeso <- RunUMAP(Endomeso, dims=1:12, reduction="harmony")
Endomeso <- FindNeighbors(object = Endomeso, dims = 1:2, reduction = "umap", k.param=10)
Endomeso <- FindClusters(Endomeso, resolution = 0.00001, algorithm = 1)
Endomeso@project.name <- "Endothelial_Mesothelial"
# Endomeso <- readRDS("eWAT_Endothelial_Mesothelial.Rds") 

## Adjust cluster labels from 0-indexed to 1-indexed
Endomeso$Label <- 0
Endomeso@meta.data[ Endomeso@meta.data$seurat_clusters == 0,"Label"] <- 1
Endomeso@meta.data[ Endomeso@meta.data$seurat_clusters == 1,"Label"] <- 2
Endomeso@meta.data[ Endomeso@meta.data$seurat_clusters == 2,"Label"] <- 3
Endomeso@meta.data[ Endomeso@meta.data$seurat_clusters == 3,"Label"] <- 4
Endomeso@meta.data[ Endomeso@meta.data$seurat_clusters == 4,"Label"] <- 5
Endomeso <- SetIdent(Endomeso, value = Endomeso$Label)

## Detection of marker genes and dietary effects
# Setup
Data <- Endomeso
Data <- SetIdent(Data, value=Data$Label)
Clusters <- sort(levels(Data@active.ident))
Result <- list()

# Loop through clusters
for (clust in Clusters) {
  ## Diet effects within the current cluster
  # Subset the dataset
  DataSubset <- subset(Data, subset = Label == clust)
  
  # Subset to expressed genes
  sca <-  FromMatrix(as.matrix(DataSubset@assays$RNA@data), DataSubset@meta.data)
  expressed_genes <- unique(sort(c(names(which(freq(sca[, sca$Diet == "LFD"]) >= 0.1)), names(which(freq(sca[, sca$Diet == "HFD"]) >= 0.1)))))  # Keeping only genes that are expressed in more than 10% of the nuclei in each diet
  
  # Extracting normalized data for expressed genes and split by replicate
  Exprs <- DataSubset@assays$RNA@data
  Exprs <- Exprs[rownames(Exprs) %in% expressed_genes, ]
  Exprs <- as.matrix(Exprs)
  Exprs_R1 <- Exprs[, colnames(Exprs) %in% rownames(DataSubset@meta.data[DataSubset@meta.data$Replicate == "R1", ])]
  Exprs_R2 <- Exprs[, colnames(Exprs) %in% rownames(DataSubset@meta.data[DataSubset@meta.data$Replicate == "R2", ])]
  
  # Setting up design matrices for differential testing of diet effects
  MD <- DataSubset@meta.data
  MD$Cell <- rownames(MD)
  Mat <- data.frame(Cell = colnames(Exprs))
  Mat <- merge(Mat, MD, by = "Cell", sort = F)
  Mat_R1 <- Mat[Mat$Cell %in% rownames(DataSubset@meta.data[DataSubset@meta.data$Replicate == "R1", ]), ]
  Mat_R2 <- Mat[Mat$Cell %in% rownames(DataSubset@meta.data[DataSubset@meta.data$Replicate == "R2", ]), ]
  Diet_R1 <- factor(as.character(as.numeric(Mat_R1$Diet == "HFD")), levels = c("0", "1"))
  Diet_R2 <- factor(as.character(as.numeric(Mat_R2$Diet == "HFD")), levels = c("0", "1"))
  Design_R1 <- model.matrix(~ Diet_R1)
  Design_R2 <- model.matrix(~ Diet_R2)
  
  # Fitting a linear model for each replicate
  Fit_R1 <- lmFit(Exprs_R1, Design_R1)
  Fit_R2 <- lmFit(Exprs_R2, Design_R2)
  
  # Computing Emperical Bayes statistics for each replicate
  Fit_R1 <- eBayes(Fit_R1, trend = T)
  Fit_R2 <- eBayes(Fit_R2, trend = T)
  
  # Differential testing in each replicate
  Result_R1 <- topTable(Fit_R1, coef = 2, adjust = "BH", number = nrow(Exprs), sort = "none")
  Result_R2 <- topTable(Fit_R2, coef = 2, adjust = "BH", number = nrow(Exprs), sort = "none")
  
  # Setting column names
  colnames(Result_R1)[c(1, 5)] <- c("logFC_Diet_R1", "Padj_Diet_R1")
  colnames(Result_R2)[c(1, 5)] <- c("logFC_Diet_R2", "Padj_Diet_R2")
  
  # Merging results from each replicate
  ResultFrame <- merge(Result_R1[, c(1,5)], Result_R2[, c(1,5)], by = 0)
  
  # Annotating diet effects
  ResultFrame$Diet <- "NS"
  ResultFrame[ResultFrame$logFC_Diet_R1 > 0 & ResultFrame$Padj_Diet_R1 <= 0.05 & ResultFrame$logFC_Diet_R2 > 0 & ResultFrame$Padj_Diet_R2 <= 0.05, "Diet"] <- "Induced"
  ResultFrame[ResultFrame$logFC_Diet_R1 < 0 & ResultFrame$Padj_Diet_R1 <= 0.05 & ResultFrame$logFC_Diet_R2 < 0 & ResultFrame$Padj_Diet_R2 <= 0.05, "Diet"] <- "Repressed"
  colnames(ResultFrame)[1] <- "Symbol"
  
  # Re-ordering columns
  ResultFrame <- ResultFrame[, c(1, 6, 2, 3, 4, 5)]
  
  ### Cluster marker genes
  # Expressed genes  
  expressed_genes <- names(which(freq(sca) >= 0.1))  # Keeping only genes that are expressed in more than 10% of the nuclei 
  
  ## One-vs-all cluster comparison
  # Extracting normalized data for expressed genes and split by replicate
  Exprs <- Data@assays$RNA@data
  Exprs <- Exprs[rownames(Exprs) %in% expressed_genes, ]
  Exprs <- as.matrix(Exprs)
  Exprs_R1 <- Exprs[, colnames(Exprs) %in% rownames(Data@meta.data[Data@meta.data$Replicate == "R1", ])]
  Exprs_R2 <- Exprs[, colnames(Exprs) %in% rownames(Data@meta.data[Data@meta.data$Replicate == "R2", ])]
  
  # Setting up design matrices for differential testing of one-vs-all cluster marker genes
  MD <- Data@meta.data
  MD$Cell <- rownames(MD)
  Mat <- data.frame(Cell = colnames(Exprs))
  Mat <- merge(Mat, MD, by = "Cell", sort = F)
  Mat_R1 <- Mat[Mat$Cell %in% rownames(Data@meta.data[Data@meta.data$Replicate == "R1", ]), ]
  Mat_R2 <- Mat[Mat$Cell %in% rownames(Data@meta.data[Data@meta.data$Replicate == "R2", ]), ]
  Clusters_R1 <- factor(as.character(as.numeric(Mat_R1$Label == clust)), levels=c("0", "1"))
  Clusters_R2 <- factor(as.character(as.numeric(Mat_R2$Label == clust)), levels=c("0", "1"))
  Design_R1 <- model.matrix(~ Clusters_R1)
  Design_R2 <- model.matrix(~ Clusters_R2)
  
  # Fitting a linear model for each replicate
  Fit_R1 <- lmFit(Exprs_R1, Design_R1)
  Fit_R2 <- lmFit(Exprs_R2, Design_R2)
  
  # Computing Emperical Bayes statistics for each replicate
  Fit_R1 <- eBayes(Fit_R1, trend = T)
  Fit_R2 <- eBayes(Fit_R2, trend = T)
  
  # Differential testing in each replicate
  Result_R1 <- topTable(Fit_R1, coef = 2, adjust = "BH", number = nrow(Exprs), sort = "none")
  Result_R2 <- topTable(Fit_R2, coef = 2, adjust = "BH", number = nrow(Exprs), sort = "none")
  
  # Setting column names
  colnames(Result_R1)[c(1, 5)] <- c("logFC_OvA_R1", "Padj_OvA_R1")
  colnames(Result_R2)[c(1, 5)] <- c("logFC_OvA_R2", "Padj_OvA_R2")
  
  # Merging results from each replicate
  OvA <- merge(Result_R1[, c(1, 5)], Result_R2[, c(1, 5)], by = 0)
  colnames(OvA)[1] <- "Symbol"
  
  # Merging with results from diet testing
  ResultFrame <- merge(ResultFrame, OvA, by = "Symbol", all = T)
  
  # Replacing NA values
  ResultFrame[is.na(ResultFrame$logFC_OvA_R1), "logFC_OvA_R1"] <- 0
  ResultFrame[is.na(ResultFrame$logFC_OvA_R2), "logFC_OvA_R2"] <- 0  
  ResultFrame[is.na(ResultFrame$Padj_OvA_R1), "Padj_OvA_R1"] <- 1
  ResultFrame[is.na(ResultFrame$Padj_OvA_R2), "Padj_OvA_R2"] <- 1
  
  ##  One-vs-one cluster comparison
  # Define comparisons
  Pairs <- data.frame(GroupA = clust, GroupB = Clusters[!(Clusters %in% clust)])
  
  # Process comparisons in a loop
  for (i in 1:nrow(Pairs)) {
    # Extracting normalized data for expressed genes and split by replicate
    Exprs <- Data@assays$RNA@data
    Exprs <- Exprs[rownames(Exprs) %in% expressed_genes, ]
    Exprs <- as.matrix(Exprs)
    Exprs_R1 <- Exprs[, colnames(Exprs) %in% rownames(Data@meta.data[Data@meta.data$Replicate == "R1" & Data@meta.data$Label %in% c(as.character(Pairs[i, 1]),as.character(Pairs[i, 2])), ])]
    Exprs_R2 <- Exprs[, colnames(Exprs) %in% rownames(Data@meta.data[Data@meta.data$Replicate == "R2" & Data@meta.data$Label %in% c(as.character(Pairs[i, 1]),as.character(Pairs[i, 2])), ])]
    
    # Setting up design matrices for differential testing of one-vs-one cluster marker genes
    MD <- Data@meta.data
    MD$Cell <- rownames(MD)
    MD <- MD[MD$Label %in% c(as.character(Pairs[i, 1]), as.character(Pairs[i, 2])), ]
    Mat <- data.frame(Cell = colnames(Exprs))
    Mat <- merge(Mat, MD, by = "Cell", sort = F)
    Mat_R1 <- Mat[Mat$Cell %in% rownames(Data@meta.data[Data@meta.data$Replicate == "R1", ]), ]
    Mat_R2 <- Mat[Mat$Cell %in% rownames(Data@meta.data[Data@meta.data$Replicate == "R2", ]), ]
    Clusters_R1 <- factor(as.character(as.numeric(Mat_R1$Label == clust)), levels = c("0", "1"))
    Clusters_R2 <- factor(as.character(as.numeric(Mat_R2$Label == clust)), levels = c("0", "1"))
    Design_R1 <- model.matrix(~ Clusters_R1)
    Design_R2 <- model.matrix(~ Clusters_R2)
    
    # Fitting a linear model for each replicate
    Fit_R1 <- lmFit(Exprs_R1, Design_R1)
    Fit_R2 <- lmFit(Exprs_R2, Design_R2)
    
    # Computing Emperical Bayes statistics for each replicate
    Fit_R1 <- eBayes(Fit_R1, trend = T)
    Fit_R2 <- eBayes(Fit_R2, trend = T)
    
    # Differential testing in each replicate
    Result_R1 <- topTable(Fit_R1, coef = 2, adjust = "BH", number = nrow(Exprs), sort = "none")
    Result_R2 <- topTable(Fit_R2, coef = 2, adjust = "BH", number = nrow(Exprs), sort = "none")
    
    # Setting column names
    colnames(Result_R1)[c(1, 5)] <- c(paste("logFC_Cluster", clust, "_vs_Cluster", as.character(Pairs[i, 2]), "_R1", sep = ""),paste("Padj_Cluster", clust, "_vs_Cluster", as.character(Pairs[i, 2]), "_R1", sep = ""))
    colnames(Result_R2)[c(1, 5)] <- c(paste("logFC_Cluster", clust, "_vs_Cluster", as.character(Pairs[i, 2]), "_R2", sep = ""),paste("Padj_Cluster", clust, "_vs_Cluster", as.character(Pairs[i, 2]), "_R2", sep = ""))
    
    # Merging results from each replicate
    OvO <- merge(Result_R1[ ,c(1,5)], Result_R2[, c(1,5)], by = 0)
    colnames(OvO)[1] <- "Symbol"
    
    # Annotating cluster marker genes
    OvO[6] <- 0
    colnames(OvO)[6] <- paste("Label_Cluster", clust, "_vs_Cluster", as.character(Pairs[i, 2]), sep = "")
    OvO[OvO[, 2] > 0 & OvO[, 3] <= 0.05 & OvO[, 4] > 0 & OvO[, 5] <= 0.05, 6] <- 1
    
    # Merging with results from diet testing and one-vs-all cluster marker testing
    ResultFrame <- merge(ResultFrame, OvO, by = "Symbol", all = T)
    
    # Replacing NA values
    ResultFrame[is.na(ResultFrame[, ncol(ResultFrame)]), ncol(ResultFrame)] <- 0
    ResultFrame[is.na(ResultFrame[, ncol(ResultFrame)-1]), ncol(ResultFrame)-1] <- 1
    ResultFrame[is.na(ResultFrame[, ncol(ResultFrame)-3]), ncol(ResultFrame)-3] <- 1
    ResultFrame[is.na(ResultFrame[, ncol(ResultFrame)-2]), ncol(ResultFrame)-2] <- 0
    ResultFrame[is.na(ResultFrame[, ncol(ResultFrame)-4]), ncol(ResultFrame)-4] <- 0
  }
  
  ## Finalize the results
  # Cluster annotations
  ResultFrame$Marker <- "NS"
  ResultFrame[ResultFrame[, 7] > 0 & ResultFrame[, 8] <= 0.05 & ResultFrame[, 9] > 0 & ResultFrame[, 10] <= 0.05, "Marker"] <- "Enriched" # Significant higher in a cluster compared to all other nuclei
  ResultFrame[which(apply(ResultFrame[, grep("Label", colnames(ResultFrame))], 1, FUN = "min") == 1), "Marker"] <- "Exclusive" # Significant higher in a cluster compared to all other clusters individually
  
  # Cleanup
  ResultFrame <- ResultFrame[, grep("Label", colnames(ResultFrame), invert = T)]
  ResultFrame <- ResultFrame[ResultFrame$Diet != "NS" | ResultFrame$Marker != "NS", ]
  ResultFrame <- ResultFrame[, c(1, ncol(ResultFrame), 2, 3:(ncol(ResultFrame)-1))]
  Result[[clust]] <- ResultFrame  ### These results are the basis for Supplementary Tables
}

# Save the results
saveRDS(Result, "eWAT_Endothelial_Mesothelial_Markers.Rds")

## Pathway analysis. NOTE: As WikiPathways is updated results may change!
# Prepare for WikiPathways analysis
gmt <- rWikiPathways::downloadPathwayArchive(organism = "Mus musculus", format = "gmt")
wp2gene <- clusterProfiler::read.gmt(gmt)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name", "version", "wpid", "org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) # TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) # TERM2NAME
Pathways <- list()

# Loop through all clusters
for (i in 1:length(Result)) {
  # Extract list
  Tmp <- Result[[i]]
  
  # Convert gene names
  Significant_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Marker != "NS",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  Enriched_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Marker == "Enriched",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  Exclusive_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Marker == "Exclusive",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  Induced_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Diet == "Induced",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  Repressed_Entrez <- clusterProfiler::bitr(Tmp[ Tmp$Diet == "Repressed",1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # Do the pathway analysis
  if (sum(wp2gene$gene %in% Significant_Entrez[[2]]) > 0) { wiki_significant <- clusterProfiler::enricher(Significant_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1,  TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_significant <- data.frame() }
  if (sum(wp2gene$gene %in% Enriched_Entrez[[2]]) > 0) { wiki_enriched <- clusterProfiler::enricher(Enriched_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_enriched <- data.frame() }
  if (sum(wp2gene$gene %in% Exclusive_Entrez[[2]]) > 0) { wiki_exclusive <- clusterProfiler::enricher(Exclusive_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_exclusive <- data.frame() }
  if (sum(wp2gene$gene %in% Induced_Entrez[[2]]) > 0) { wiki_induced <- clusterProfiler::enricher(Induced_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_induced <- data.frame() }
  if (sum(wp2gene$gene %in% Repressed_Entrez[[2]]) > 0) { wiki_repressed <- clusterProfiler::enricher(Repressed_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) } else { wiki_repressed <- data.frame() }
  
  # Adding gene symbols to the resulting pathway file
  if (nrow(wiki_significant) > 0) { wiki_significant <- as.data.frame(DOSE::setReadable(wiki_significant, org.Mm.eg.db, keyType = "ENTREZID")) }
  if (nrow(wiki_enriched) > 0) { wiki_enriched <- as.data.frame(DOSE::setReadable(wiki_enriched, org.Mm.eg.db, keyType = "ENTREZID")) }
  if (nrow(wiki_exclusive) > 0) { wiki_exclusive <- as.data.frame(DOSE::setReadable(wiki_exclusive, org.Mm.eg.db, keyType = "ENTREZID")) }
  if (nrow(wiki_induced) > 0) { wiki_induced <- as.data.frame(DOSE::setReadable(wiki_induced, org.Mm.eg.db, keyType = "ENTREZID")) }
  if (nrow(wiki_repressed) > 0) { wiki_repressed <- as.data.frame(DOSE::setReadable(wiki_repressed, org.Mm.eg.db, keyType = "ENTREZID")) }
  
  # Calculate gene ratios (i.e. the fraction of genes in the pathway)
  if (nrow(wiki_significant) > 0) { wiki_significant$Significant_GeneRatio <- as.numeric(substr(wiki_significant$GeneRatio, 1, regexpr("/", wiki_significant$GeneRatio)-1)) / as.numeric(substr(wiki_significant$GeneRatio, regexpr("/", wiki_significant$GeneRatio)+1, nchar(wiki_significant$GeneRatio))) }
  if (nrow(wiki_enriched) > 0) { wiki_enriched$Enriched_GeneRatio <- as.numeric(substr(wiki_enriched$GeneRatio, 1, regexpr("/", wiki_enriched$GeneRatio)-1)) / as.numeric(substr(wiki_enriched$GeneRatio, regexpr("/", wiki_enriched$GeneRatio)+1, nchar(wiki_enriched$GeneRatio))) }
  if (nrow(wiki_exclusive) > 0) { wiki_exclusive$Exclusive_GeneRatio <- as.numeric(substr(wiki_exclusive$GeneRatio, 1, regexpr("/", wiki_exclusive$GeneRatio)-1)) / as.numeric(substr(wiki_exclusive$GeneRatio, regexpr("/", wiki_exclusive$GeneRatio)+1, nchar(wiki_exclusive$GeneRatio))) }
  if (nrow(wiki_induced) > 0) { wiki_induced$Induced_GeneRatio <- as.numeric(substr(wiki_induced$GeneRatio, 1, regexpr("/", wiki_induced$GeneRatio)-1)) / as.numeric(substr(wiki_induced$GeneRatio, regexpr("/", wiki_induced$GeneRatio)+1, nchar(wiki_induced$GeneRatio))) }
  if (nrow(wiki_repressed) > 0) { wiki_repressed$Repressed_GeneRatio <- as.numeric(substr(wiki_repressed$GeneRatio, 1, regexpr("/", wiki_repressed$GeneRatio)-1)) / as.numeric(substr(wiki_repressed$GeneRatio, regexpr("/", wiki_repressed$GeneRatio)+1, nchar(wiki_repressed$GeneRatio))) }
  
  # Calculate pathway ratios (i.e. the fraction of the pathway in the genes)
  if (nrow(wiki_significant) > 0) { wiki_significant$Significant_PathRatio <- as.numeric(substr(wiki_significant$GeneRatio, 1, regexpr("/", wiki_significant$GeneRatio)-1)) / as.numeric(substr(wiki_significant$BgRatio, 1, regexpr("/", wiki_significant$BgRatio)-1)) }
  if (nrow(wiki_enriched) > 0) { wiki_enriched$Enriched_PathRatio <- as.numeric(substr(wiki_enriched$GeneRatio, 1, regexpr("/", wiki_enriched$GeneRatio)-1)) / as.numeric(substr(wiki_enriched$BgRatio, 1, regexpr("/", wiki_enriched$BgRatio)-1)) }
  if (nrow(wiki_exclusive) > 0) { wiki_exclusive$Exclusive_PathRatio <- as.numeric(substr(wiki_exclusive$GeneRatio, 1, regexpr("/", wiki_exclusive$GeneRatio)-1)) / as.numeric(substr(wiki_exclusive$BgRatio, 1, regexpr("/", wiki_exclusive$BgRatio)-1)) }
  if (nrow(wiki_induced) > 0) { wiki_induced$Induced_PathRatio <- as.numeric(substr(wiki_induced$GeneRatio, 1, regexpr("/", wiki_induced$GeneRatio)-1)) / as.numeric(substr(wiki_induced$BgRatio, 1, regexpr("/", wiki_induced$BgRatio)-1)) }
  if (nrow(wiki_repressed) > 0) { wiki_repressed$Repressed_PathRatio <- as.numeric(substr(wiki_repressed$GeneRatio, 1, regexpr("/", wiki_repressed$GeneRatio)-1)) / as.numeric(substr(wiki_repressed$BgRatio, 1, regexpr("/", wiki_repressed$BgRatio)-1)) }
  
  # Set column names
  if (nrow(wiki_significant) > 0) { colnames(wiki_significant)[c(5,7,8)] <- c("Significant_Pvalue","Significant_FDR","Significant_Genes") }
  if (nrow(wiki_enriched) > 0) { colnames(wiki_enriched)[c(5,7,8)] <- c("Enriched_Pvalue","Enriched_FDR","Enriched_Genes") }
  if (nrow(wiki_exclusive) > 0) { colnames(wiki_exclusive)[c(5,7,8)] <- c("Exclusive_Pvalue","Exclusive_FDR","Exclusive_Genes") }
  if (nrow(wiki_induced) > 0) { colnames(wiki_induced)[c(5,7,8)] <- c("Induced_Pvalue","Induced_FDR","Induced_Genes") }
  if (nrow(wiki_repressed) > 0) { colnames(wiki_repressed)[c(5,7,8)] <- c("Repressed_Pvalue","Repressed_FDR","Repressed_Genes") }
  
  # Combine the results
  significant <- c(wiki_significant[ wiki_significant$Significant_FDR <= 0.05,1],wiki_enriched[ wiki_enriched$Enriched_FDR <= 0.05,1],  wiki_exclusive[ wiki_exclusive$Exclusive_FDR <= 0.05,1],  wiki_induced[ wiki_induced$Induced_FDR <= 0.05,1],  wiki_repressed[ wiki_repressed$Repressed_FDR <= 0.05,1])
  if (length(significant) > 0) {
    # Subset from all pathways to only significant ones
    all <- wpid2name[ wpid2name[,1] %in% significant,]
    all <- all[ duplicated(all[,1])==F,]
    colnames(all) <- c("ID", "Description")
    
    # Merge statistics
    if (nrow(wiki_significant) > 0) { all <- merge(all, wiki_significant[,c(1,10,11,5,7,8)], all.x=T, by="ID") }
    if (nrow(wiki_enriched) > 0) { all <- merge(all, wiki_enriched[,c(1,10,11,5,7,8)], all.x=T, by="ID") }
    if (nrow(wiki_exclusive) > 0) { all <- merge(all, wiki_exclusive[,c(1,10,11,5,7,8)], all.x=T, by="ID") }
    if (nrow(wiki_induced) > 0) { all <- merge(all, wiki_induced[,c(1,10,11,5,7,8)], all.x=T, by="ID") }
    if (nrow(wiki_repressed) > 0) { all <- merge(all, wiki_repressed[,c(1,10,11,5,7,8)], all.x=T, by="ID") }
    
    # Handle NAs appropriately (ratios to 0, pvalues/FDRs to 1 and gene lists to "")
    for (m in grep("Ratio", colnames(all))) { all[ is.na(all[,m]),m] <- 0 }
    for (m in grep("Pvalue", colnames(all))) { all[ is.na(all[,m]),m] <- 1 }
    for (m in grep("FDR", colnames(all))) { all[ is.na(all[,m]),m] <- 1 }
    for (m in grep("Genes", colnames(all))) { all[ is.na(all[,m]),m] <- "" }
    
    # Store the results
    Pathways[[i]] <- all # The list is the basis for the Supplementary Tables
  }
}

## Label clusters based on marker genes and pathway analyses
# Labels
Endomeso$Subtype <- "NA"
Endomeso@meta.data[Endomeso@meta.data$Label %in% 1, "Subtype"] <- "MC"
Endomeso@meta.data[Endomeso@meta.data$Label %in% 2, "Subtype"] <- "IMC"
Endomeso@meta.data[Endomeso@meta.data$Label %in% 3, "Subtype"] <- "EPC"
Endomeso@meta.data[Endomeso@meta.data$Label %in% 4, "Subtype"] <- "VEC"
Endomeso@meta.data[Endomeso@meta.data$Label %in% 5, "Subtype"] <- "LEC"
Endomeso <- SetIdent(Endomeso, value = Endomeso$Subtype)

# Transfer labels to the eWAT object
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Subtype == "MC",]),"Subtype"] <- "MC"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Subtype == "IMC",]),"Subtype"] <- "IMC"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Subtype == "EPC",]),"Subtype"] <- "EPC"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Subtype == "VEC",]),"Subtype"] <- "VEC"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Subtype == "LEC",]),"Subtype"] <- "LEC"

# Save results
saveRDS(Endomeso, "eWAT_Endothelial_Mesothelial.Rds") ## This object is downloadable from Open Science Framework
saveRDS(eWAT, "eWAT_Annotated.Rds") ## This object is downloadable from Open Science Framework with additional annotations added by subclustering each major cell type

## DimPlots
DimPlot(Endomeso, label = T)
DimPlot(Endomeso, label = T, split.by="Diet")
DimPlot(Endomeso, label = T, split.by="Replicate")

## Marker gene heatmap
# Select marker genes for each cell type
MC <- data.frame(gene = c("Msln", "Upk3b", "Gpm6a"), type = "MC")
IMC <- data.frame(gene = c("Ccr2", "Alcam", "Cd200r1"), type = "IMC")
EPC <- data.frame(gene = c("Kdr", "Cdh5", "Flt1"), type = "EPC")
VEC <- data.frame(gene = c("Vegfc", "Vwf", "Vcam1"), type = "VEC")
LEC <- data.frame(gene = c("Lepr", "Lyve1", "Ccl21a"), type = "LEC")

# Combine results
Averages <- rbind(MC, IMC, EPC, VEC, LEC)

# Calculate average expression for selected marker genes
for (i in 1:nrow(Averages)) {
  Averages[i,3] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "1" & Endomeso@meta.data$Replicate == "R1",])])))
  Averages[i,4] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "1" & Endomeso@meta.data$Replicate == "R2",])])))
  Averages[i,5] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "2" & Endomeso@meta.data$Replicate == "R1",])])))
  Averages[i,6] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "2" & Endomeso@meta.data$Replicate == "R2",])])))
  Averages[i,7] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "3" & Endomeso@meta.data$Replicate == "R1",])])))
  Averages[i,8] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "3" & Endomeso@meta.data$Replicate == "R2",])])))
  Averages[i,9] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "4" & Endomeso@meta.data$Replicate == "R1",])])))
  Averages[i,10] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "4" & Endomeso@meta.data$Replicate == "R2",])])))
  Averages[i,11] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "5" & Endomeso@meta.data$Replicate == "R1",])])))
  Averages[i,12] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Label == "5" & Endomeso@meta.data$Replicate == "R2",])])))
}
rownames(Averages) <- Averages[,1]

# Scale each experiment
Exp1 <- Averages[,c(3,5,7,9,11)]
Exp1 <- t(scale(t(Exp1)))
Exp2 <- Averages[,c(4,6,8,10,12)]
Exp2 <- t(scale(t(Exp2)))
Averages[,c(3,5,7,9,11)] <- Exp1
Averages[,c(4,6,8,10,12)] <- Exp2

# Plot the heatmap
Heatmap(Averages[,c(3:12)], cluster_columns=F, cluster_rows=F)

## Gene module scoring
# Extract marker gene results
Clust1 <- Result[[1]] 
Clust2 <- Result[[2]] 
Clust3 <- Result[[3]] 
Clust4 <- Result[[4]] 
Clust5 <- Result[[5]] 

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
Endomeso <- AddModuleScore(Endomeso, features = list(Cluster1 = Clust1[order(-Clust1$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster2 = Clust2[order(-Clust2$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster3 = Clust3[order(-Clust3$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster4 = Clust4[order(-Clust4$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster5 = Clust5[order(-Clust5$logFC_OvO), ][1:50, ]$Symbol))

# Scale the module scores (across all nuclei)
Endomeso$Cluster1 <- scale(Endomeso$Cluster1)
Endomeso$Cluster2 <- scale(Endomeso$Cluster2)
Endomeso$Cluster3 <- scale(Endomeso$Cluster3)
Endomeso$Cluster4 <- scale(Endomeso$Cluster4)
Endomeso$Cluster5 <- scale(Endomeso$Cluster5)

# Computing mean scaled module scores in each cluster in each replicate
Averages <- as.data.frame(matrix(ncol=length(Result), nrow=length(Result)*2))

# Set column and row names
Labels <- Endomeso@meta.data[, c("Subtype","Label")]
Labels <- Labels[ duplicated(Labels$Label)==F,]
Labels <- Labels[ order(Labels$Label),]
colnames(Averages) <- paste(Labels$Subtype,"module",sep="_")
rownames(Averages)[seq(1,(nrow(Averages)-1),by=2)] <- paste(Labels$Subtype,"R1",sep="_")
rownames(Averages)[seq(2,nrow(Averages),by=2)] <- paste(Labels$Subtype,"R2",sep="_")

# Calculate averages
for (module in 1:length(Result)) {
  Counter <- 1
  for (label in 1:length(Result)) {
    for (rep in c("R1","R2")) {
      Averages[Counter,module] <- mean(Endomeso@meta.data[ Endomeso@meta.data$Label == label & Endomeso$Replicate == rep,paste("Cluster",module,sep="")])
      Counter <- Counter + 1
    }
  }
}

## Plot mean scaled module scores
col <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(as.matrix(Averages), cluster_columns = F, cluster_rows = F, col = col)

## Subtype proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(Endomeso$Subtype)), nrow=4))
colnames(Composition) <- unique(Endomeso$Subtype)
for (i in 1:length(unique(Endomeso$Subtype))) {
  Composition[1,i] <- nrow(Endomeso@meta.data[ Endomeso@meta.data$Dataset == "LFD_R1" & Endomeso@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Endomeso@meta.data[ Endomeso@meta.data$Dataset == "LFD_R1",])
  Composition[2,i] <- nrow(Endomeso@meta.data[ Endomeso@meta.data$Dataset == "LFD_R2" & Endomeso@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Endomeso@meta.data[ Endomeso@meta.data$Dataset == "LFD_R2",])
  Composition[3,i] <- nrow(Endomeso@meta.data[ Endomeso@meta.data$Dataset == "HFD_R1" & Endomeso@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Endomeso@meta.data[ Endomeso@meta.data$Dataset == "HFD_R1",])
  Composition[4,i] <- nrow(Endomeso@meta.data[ Endomeso@meta.data$Dataset == "HFD_R2" & Endomeso@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Endomeso@meta.data[ Endomeso@meta.data$Dataset == "HFD_R2",])
}

# Calculate stats (mean, SEM, P)
Stats <- data.frame(matrix(ncol=5, nrow=length(unique(Endomeso$Subtype))))
colnames(Stats) <- c("LFD_mean","LFD_SEM","HFD_mean","HFD_SEM","Pvalue")
rownames(Stats) <- colnames(Composition)
for (i in 1:length(unique(Endomeso$Subtype))) {
  Stats[i,1] <- mean(as.matrix(Composition[1:2,i]))
  Stats[i,2] <- sd(as.matrix(Composition[1:2,i])) / sqrt(2)
  Stats[i,3] <- mean(as.matrix(Composition[3:4,i]))
  Stats[i,4] <- sd(as.matrix(Composition[3:4,i])) / sqrt(2)
  Stats[i,5] <- t.test(as.matrix(Composition[1:2,i], as.matrix(Composition[3:4,i])))$p.value
}

# Plot the composition
Blt <- barplot(t(as.matrix(Stats[,c("LFD_mean","HFD_mean")])), beside=T, las=1, ylim=c(0,1))
arrows(Blt[1,], Stats$LFD_mean - Stats$LFD_SEM, Blt[1,], Stats$LFD_mean + Stats$LFD_SEM, angle=90, code=3, length=0.05)
arrows(Blt[2,], Stats$HFD_mean - Stats$HFD_SEM, Blt[2,], Stats$HFD_mean + Stats$HFD_SEM, angle=90, code=3, length=0.05)

### Expression of cytokines 
## Extract gene names for chemokines using biomaRt
ensembl <- useMart("ensembl" ,dataset="mmusculus_gene_ensembl")
res <- getBM(c("external_gene_name"), filters = "go", values = "GO:0005125", mart = ensembl)
Chemokines <- res[,1]

## Process counts
# Calculate average expression of genes 
Averages <- data.frame(gene = c(Chemokines))
for (i in 1:nrow(Averages)) {
  Averages[i,2] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Replicate == "R1" & Endomeso@meta.data$Diet == "LFD",])])))
  Averages[i,3] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Replicate == "R1" & Endomeso@meta.data$Diet == "HFD",])])))
  Averages[i,4] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Replicate == "R2" & Endomeso@meta.data$Diet == "LFD",])])))
  Averages[i,5] <- log1p(mean(expm1(Endomeso@assays$RNA@data[ rownames(Endomeso@assays$RNA@data) %in% Averages[i,1], colnames(Endomeso@assays$RNA@data) %in% rownames(Endomeso@meta.data[ Endomeso@meta.data$Replicate == "R2" & Endomeso@meta.data$Diet == "HFD",])])))
}
rownames(Averages) <- Averages[,1]

# Remove genes without signal
Averages <- Averages[ !is.na(Averages[,2]),]

# Calculate log fold changes between LFD and HFD in each replicate
Averages[,6] <- Averages[,3] - Averages[,2]
Averages[,7] <- Averages[,5] - Averages[,4]

# Calculate the average logFC and the SEM
Averages[,8] <- apply(Averages[,c(6,7)],1,FUN="mean")
Averages[,9] <- apply(Averages[,c(6,7)],1,FUN="sd")/sqrt(2)

# Get the maximum expression across replicates
Averages[,10] <- apply(cbind(rowMeans(Averages[,c(2,4)]),rowMeans(Averages[,c(3,5)])),1,FUN="max")

## Scatterplot the results
plot(Averages[,10], Averages[,8], xlab="Maximum expression", ylab="log2FC", pch=16, las=1, main="Chemokines")
abline(h=log2(1.1), lty=2)
abline(h=0)
abline(h=-log2(1.1), lty=2)
```
[Back to start](../README.md)<br>