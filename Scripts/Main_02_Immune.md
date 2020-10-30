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
Immune <- subset(eWAT, subset = Annotation == "Immune")

# Embedding and clustering of immune cells into cell types and subpopulations  
# NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN. 
# Thus, in order to reproduce downstream analyses using the same Harmony results, please download eWAT_Immune.Rds
Immune <- FindVariableFeatures(Immune, nfeatures = 1000)
VariableFeatures(Immune) <- VariableFeatures(Immune)[!(VariableFeatures(Immune) %in% Amb)]
Immune <- ScaleData(Immune)
Immune <- RunPCA(Immune)
Immune <- RunHarmony(Immune, group.by.vars = "Dataset") 
Immune <- RunUMAP(Immune, dims = 1:16, reduction = "harmony")
Immune <- FindNeighbors(object = Immune, dims = 1:2, reduction = "umap", k.param = 50)
Immune <- FindClusters(Immune, resolution = 0.1, algorithm = 1)
Immune@project.name <- "Immune"
# Immune <- readRDS("eWAT_Immune.Rds") 

## Adjust cluster labels from 0-indexed to 1-indexed
Immune$Label <- 0
Immune@meta.data[Immune@meta.data$seurat_clusters %in% 0, "Label"] <- 1
Immune@meta.data[Immune@meta.data$seurat_clusters %in% 1, "Label"] <- 2
Immune@meta.data[Immune@meta.data$seurat_clusters %in% 2, "Label"] <- 3
Immune@meta.data[Immune@meta.data$seurat_clusters %in% 3, "Label"] <- 5
Immune@meta.data[Immune@meta.data$seurat_clusters %in% 4, "Label"] <- 8
Immune@meta.data[Immune@meta.data$seurat_clusters %in% 5, "Label"] <- 9
Immune@meta.data[Immune@meta.data$seurat_clusters %in% 6, "Label"] <- 7
Immune@meta.data[Immune@meta.data$seurat_clusters %in% 7, "Label"] <- 6
Immune@meta.data[Immune@meta.data$seurat_clusters %in% 8, "Label"] <- 4
Immune <- SetIdent(Immune, value = Immune$Label)

## Detection of marker genes and dietary effects
# Setup
Data <- Immune
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
saveRDS(Result, "eWAT_Immune_Markers.Rds")

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

## Label clusters based on marker genes and pathway analyses (there are several cell types, and macrophages have additional subtypes. Start byt looking at cell types)
# Labels
Immune$Subtype <- "NA"
Immune@meta.data[Immune@meta.data$Label %in% c(1,2,3,4,5,6), "Subtype"] <- "Macrophages"
Immune@meta.data[Immune@meta.data$Label %in% 7, "Subtype"] <- "B-cells"
Immune@meta.data[Immune@meta.data$Label %in% 8, "Subtype"] <- "T-cells"
Immune@meta.data[Immune@meta.data$Label %in% 9, "Subtype"] <- "DCs"
Immune <- SetIdent(Immune, value = Immune$Subtype)

# Transfer labels to the eWAT object
eWAT$Subtype <- eWAT$Annotation
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Immune@meta.data[ Immune@meta.data$Subtype == "Macrophages",]),"Subtype"] <- "Macrophages"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Immune@meta.data[ Immune@meta.data$Subtype == "B-cells",]),"Subtype"] <- "B-cells"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Immune@meta.data[ Immune@meta.data$Subtype == "T-cells",]),"Subtype"] <- "T-cells"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Immune@meta.data[ Immune@meta.data$Subtype == "DCs",]),"Subtype"] <- "DCs"

# Save results
saveRDS(Immune, "eWAT_Immune.Rds") ## This object is downloadable from Open Science Framework

## DimPlots (using labels to not group macrophages in the plot)
DimPlot(Immune, label = T, group.by = "Label")

## Marker gene heatmap
# Select marker genes for each cell type
Macrophages <- data.frame(gene = c("Adgre1", "Lyz2", "Ccl6"), type = "Macrophage")
Tcells <- data.frame(gene = c("Skap1", "Lck", "Themis"), type = "Tcells")
Bcells <- data.frame(gene = c("Cd79a", "Cd79b", "Ms4a1"), type = "Bcells")
DCs <- data.frame(gene = c("Cd209a", "Flt3", "Cd244a"), type = "DCs")

# Combine results
Averages <- rbind(Macrophages, Tcells, Bcells, DCs)

# Calculating fractional expressions of marker genes for each replicate
# RATIONALE: There are multiple subtypes of macrophages that express macrophage marker to different levels, but with the same proportions. Thus, showing the fraction of nuclei expressing these markers highlights the cell types rather than differences between macrophages (which is explored later).
for (i in 1:nrow(Averages)) {
  for (q in 1:9) {
    Averages[i, q+2] <- sum(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) == Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == q & Immune@meta.data$Replicate == "R1", ])] > 0)/nrow(Immune@meta.data[Immune@meta.data$Label == q & Immune@meta.data$Replicate == "R1", ])
  }
}
for (i in 1:nrow(Averages)) {
  for (q in 1:9) {
    Averages[i, q+11] <- sum(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) == Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == q & Immune@meta.data$Replicate == "R2", ])] > 0)/nrow(Immune@meta.data[Immune@meta.data$Label == q & Immune@meta.data$Replicate == "R2", ])
  }
}
rownames(Averages) <- Averages$gene

# Plot the heatmap (scaling the percentages)
Heatmap(t(scale(t(Averages[, c(3:20)]))), cluster_columns = T, cluster_rows = F)

## Gene module scoring
# Extract marker gene results
Clust1 <- Result[[1]] 
Clust2 <- Result[[2]] 
Clust3 <- Result[[3]] 
Clust4 <- Result[[4]] 
Clust5 <- Result[[5]] 
Clust6 <- Result[[6]] 
Clust7 <- Result[[7]]
Clust8 <- Result[[8]] 
Clust9 <- Result[[9]]

# Subset to genes that are enriched or exclusive (!= NS)
Clust1 <- Clust1[Clust1$Marker != "NS", ]
Clust2 <- Clust2[Clust2$Marker != "NS", ]
Clust3 <- Clust3[Clust3$Marker != "NS", ]
Clust4 <- Clust4[Clust4$Marker != "NS", ]
Clust5 <- Clust5[Clust5$Marker != "NS", ]
Clust6 <- Clust6[Clust6$Marker != "NS", ]
Clust7 <- Clust7[Clust7$Marker != "NS", ]
Clust8 <- Clust8[Clust8$Marker != "NS", ]
Clust9 <- Clust9[Clust9$Marker != "NS", ]

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
Immune <- AddModuleScore(Immune, features = list(Cluster1 = Clust1[order(-Clust1$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster2 = Clust2[order(-Clust2$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster3 = Clust3[order(-Clust3$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster4 = Clust4[order(-Clust4$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster5 = Clust5[order(-Clust5$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster6 = Clust6[order(-Clust6$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster7 = Clust7[order(-Clust7$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster8 = Clust8[order(-Clust8$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster9 = Clust9[order(-Clust9$logFC_OvO), ][1:50, ]$Symbol))

# Scale the module scores (across all nuclei)
Immune$Cluster1 <- scale(Immune$Cluster1)
Immune$Cluster2 <- scale(Immune$Cluster2)
Immune$Cluster3 <- scale(Immune$Cluster3)
Immune$Cluster4 <- scale(Immune$Cluster4)
Immune$Cluster5 <- scale(Immune$Cluster5)
Immune$Cluster6 <- scale(Immune$Cluster6)
Immune$Cluster7 <- scale(Immune$Cluster7)
Immune$Cluster8 <- scale(Immune$Cluster8)
Immune$Cluster9 <- scale(Immune$Cluster9)

# Computing mean scaled module scores in each cluster in each replicate
Averages <- as.data.frame(matrix(ncol=length(Result), nrow=length(Result)*2))

# Set column and row names
Labels <- Immune@meta.data[, c("Subtype","Label")]
Labels <- Labels[ duplicated(Labels$Label)==F,]
Labels <- Labels[ order(Labels$Label),]
Labels[1:6,1] <- paste("Macrophage", seq(1,6), sep="")
colnames(Averages) <- paste(Labels$Subtype,"module",sep="_")
rownames(Averages)[seq(1,(nrow(Averages)-1),by=2)] <- paste(Labels$Subtype,"R1",sep="_")
rownames(Averages)[seq(2,nrow(Averages),by=2)] <- paste(Labels$Subtype,"R2",sep="_")

# Calculate averages
for (module in 1:length(Result)) {
  Counter <- 1
  for (label in 1:length(Result)) {
    for (rep in c("R1","R2")) {
      Averages[Counter,module] <- mean(Immune@meta.data[ Immune@meta.data$Label == label & Immune$Replicate == rep,paste("Cluster",module,sep="")])
      Counter <- Counter + 1
    }
  }
}

## Plot mean scaled module scores
col <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(as.matrix(Averages), cluster_columns = F, cluster_rows = F, col = col)

## Subtype proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(Immune$Subtype)), nrow=4))
colnames(Composition) <- unique(Immune$Subtype)
for (i in 1:length(unique(Immune$Subtype))) {
  Composition[1,i] <- nrow(Immune@meta.data[ Immune@meta.data$Dataset == "LFD_R1" & Immune@meta.data$Subtype == colnames(Composition)[i],]) / nrow(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R1",])
  Composition[2,i] <- nrow(Immune@meta.data[ Immune@meta.data$Dataset == "LFD_R2" & Immune@meta.data$Subtype == colnames(Composition)[i],]) / nrow(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R2",])
  Composition[3,i] <- nrow(Immune@meta.data[ Immune@meta.data$Dataset == "HFD_R1" & Immune@meta.data$Subtype == colnames(Composition)[i],]) / nrow(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R1",])
  Composition[4,i] <- nrow(Immune@meta.data[ Immune@meta.data$Dataset == "HFD_R2" & Immune@meta.data$Subtype == colnames(Composition)[i],]) / nrow(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R2",])
}

# Calculate stats (mean, SEM, P)
Stats <- data.frame(matrix(ncol=5, nrow=length(unique(Immune$Subtype))))
colnames(Stats) <- c("LFD_mean","LFD_SEM","HFD_mean","HFD_SEM","Pvalue")
rownames(Stats) <- colnames(Composition)
for (i in 1:length(unique(Immune$Subtype))) {
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

### Targeted analysis of macrophages
## Subset to only macrophages
Macrophages <- subset(Immune, subset = Subtype == "Macrophages")

## Cross-study comparison (Chakarov et al. and Jaitin et al.)
# Add module scores using marker genes from the supplementary and main figures of Chakarov et al.
Macrophages <- AddModuleScore(Macrophages, features= list(High = c("Sepp1","F13a1","Ctsb","Itm2b","Folr2","Fcna","Cbr2","Grn","Fcgrt","Hmox1","Wfdc17","Pf4","Gas6","Ms4a4a","Cd209f","Pltp","Pmp22","Lyve1","Ccl8","Ccl7","Cfp","Tgfbi","Rnase4","Cd163","Ninj1"), Low = c("Ear2","Fcrls","Coro1a","Slamf9","Axl","Tmem176b","H2-K1","Mpeg1","Apoe","Ctsz","Cd14","Cx3cr1","Tmsb4x","Hexb","H2-DMb1","H2-DMa","Cts8","Lyz1","Psap","H2-Aa","H2-Eb1","H2-Ab1","Cd74","Chil3","Plbd1")))
Macrophages$Cluster1 <- scale(Macrophages$Cluster1)
Macrophages$Cluster2 <- scale(Macrophages$Cluster2)

# Boxplots and FeaturePlots
boxplot(Cluster1 ~ Label, data = Macrophages@meta.data, outline=F)
boxplot(Cluster2 ~ Label, data = Macrophages@meta.data, outline=F)
FeaturePlot(Macrophages, "Cluster1", min.cutoff = -1)
FeaturePlot(Macrophages, "Cluster2", min.cutoff = -1)

# Add module scores using marker genes from the supplementary and main figures of Jaitun et al.
LAM_markers <- read.delim("LAM.tx", header=T, dec=",")
LAM_markers <- LAM_markers[ LAM_markers[,2] <= 0.05,]
Macrophages <- AddModuleScore(Macrophages, features= list(High = LAM_markers[ LAM_markers[,3] >= 1,1], Low = LAM_markers[ LAM_markers[,3] <= -1,1]))
Macrophages$Cluster1 <- scale(Macrophages$Cluster1)
Macrophages$Cluster2 <- scale(Macrophages$Cluster2)

# Boxplots and FeaturePlots
boxplot(Cluster1 ~ Label, data = Macrophages@meta.data, outline=F)
boxplot(Cluster2 ~ Label, data = Macrophages@meta.data, outline=F)
FeaturePlot(Macrophages, "Cluster1", min.cutoff = -1)
FeaturePlot(Macrophages, "Cluster2", min.cutoff = -1)

## Label clusters based on marker genes, pathway analyses and cross-study comparisons
# Labels
Macrophages$Subtype <- "NA"
Macrophages@meta.data[Macrophages@meta.data$Label %in% 1, "Subtype"] <- "PVM"
Macrophages@meta.data[Macrophages@meta.data$Label %in% 2, "Subtype"] <- "LAM"
Macrophages@meta.data[Macrophages@meta.data$Label %in% 3, "Subtype"] <- "NPVM"
Macrophages@meta.data[Macrophages@meta.data$Label %in% 4, "Subtype"] <- "CEM"
Macrophages@meta.data[Macrophages@meta.data$Label %in% 5, "Subtype"] <- "P-LAM"
Macrophages@meta.data[Macrophages@meta.data$Label %in% 6, "Subtype"] <- "RM"
Macrophages <- SetIdent(Macrophages, value = Macrophages$Subtype)

# Transfer labels to the eWAT object
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Subtype == "PVM",]),"Subtype"] <- "PVM"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Subtype == "LAM",]),"Subtype"] <- "LAM"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Subtype == "NPVM",]),"Subtype"] <- "NPVM"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Subtype == "CEM",]),"Subtype"] <- "CEM"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Subtype == "P-LAM",]),"Subtype"] <- "P-LAM"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Subtype == "RM",]),"Subtype"] <- "RM"

# Save results
saveRDS(eWAT, "eWAT_Annotated.Rds") ## This object is downloadable from Open Science Framework with additional annotations added by subclustering each major cell type

## DimPlots
DimPlot(Macrophages, label=T)
DimPlot(Macrophages, label=T, split.by="Diet")
DimPlot(Macrophages, label=T, split.by="Replicate")

## Marker gene heatmap
# Select marker genes for each cell type
PVM <- data.frame(gene = c("Mrc1","Lyve1","Cd163"), type = "PVM")
LAM <- data.frame(gene = c("Lpl","Trem2","Cd9"), type = "LAM")
NPVM <- data.frame(gene = c("Ear2","Fcrls","Cd74"), type = "NPVM")
CEM <- data.frame(gene = c("Col5a2","Tgfbr3","Col3a1"), type = "CEM")
PLAM <- data.frame(gene = c("Kif15","Kif11","Pola1"), type = "PLAM") 
RM <- data.frame(gene = c("Prg4", "Tgfb2","Ltbp1"), type = "RM")

# Combine results
Averages <- rbind(PVM,LAM,NPVM,CEM, PLAM, RM)

# Calculate average expression for selected marker genes
for (i in 1:nrow(Averages)) {
  Averages[i,3] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "1" & Macrophages@meta.data$Replicate == "R1",])])))
  Averages[i,4] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "1" & Macrophages@meta.data$Replicate == "R2",])])))
  Averages[i,5] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "2" & Macrophages@meta.data$Replicate == "R1",])])))
  Averages[i,6] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "2" & Macrophages@meta.data$Replicate == "R2",])])))
  Averages[i,7] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "3" & Macrophages@meta.data$Replicate == "R1",])])))
  Averages[i,8] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "3" & Macrophages@meta.data$Replicate == "R2",])])))
  Averages[i,9] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "4" & Macrophages@meta.data$Replicate == "R1",])])))
  Averages[i,10] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "4" & Macrophages@meta.data$Replicate == "R2",])])))
  Averages[i,11] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "5" & Macrophages@meta.data$Replicate == "R1",])])))
  Averages[i,12] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "5" & Macrophages@meta.data$Replicate == "R2",])])))
  Averages[i,13] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "6" & Macrophages@meta.data$Replicate == "R1",])])))
  Averages[i,14] <- log1p(mean(expm1(Macrophages@assays$RNA@data[ rownames(Macrophages@assays$RNA@data) %in% Averages[i,1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[ Macrophages@meta.data$Label == "6" & Macrophages@meta.data$Replicate == "R2",])])))
}
rownames(Averages) <- Averages[,1]

# Scale each experiment
Exp1 <- Averages[,c(3,5,7,9,11,13)]
Exp1 <- t(scale(t(Exp1)))
Exp2 <- Averages[,c(4,6,8,10,12,14)]
Exp2 <- t(scale(t(Exp2)))
Averages[,c(3,5,7,9,11,13)] <- Exp1
Averages[,c(4,6,8,10,12,14)] <- Exp2

# Plot the heatmap
Heatmap(Averages[,c(3:14)], cluster_columns=F)

## Subtype proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(Macrophages$Subtype)), nrow=4))
colnames(Composition) <- unique(Macrophages$Subtype)
for (i in 1:length(unique(Macrophages$Subtype))) {
  Composition[1,i] <- nrow(Macrophages@meta.data[ Macrophages@meta.data$Dataset == "LFD_R1" & Macrophages@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Macrophages@meta.data[ Macrophages@meta.data$Dataset == "LFD_R1",])
  Composition[2,i] <- nrow(Macrophages@meta.data[ Macrophages@meta.data$Dataset == "LFD_R2" & Macrophages@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Macrophages@meta.data[ Macrophages@meta.data$Dataset == "LFD_R2",])
  Composition[3,i] <- nrow(Macrophages@meta.data[ Macrophages@meta.data$Dataset == "HFD_R1" & Macrophages@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Macrophages@meta.data[ Macrophages@meta.data$Dataset == "HFD_R1",])
  Composition[4,i] <- nrow(Macrophages@meta.data[ Macrophages@meta.data$Dataset == "HFD_R2" & Macrophages@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Macrophages@meta.data[ Macrophages@meta.data$Dataset == "HFD_R2",])
}

# Calculate stats (mean, SEM, P)
Stats <- data.frame(matrix(ncol=5, nrow=length(unique(Macrophages$Subtype))))
colnames(Stats) <- c("LFD_mean","LFD_SEM","HFD_mean","HFD_SEM","Pvalue")
rownames(Stats) <- colnames(Composition)
for (i in 1:length(unique(Macrophages$Subtype))) {
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

### MacSpectrum analysis
## Find genes expressed in at least 10% of nuclei in any given macrophage subtype
sca <-  FromMatrix(as.matrix(Macrophages@assays$RNA@data), Macrophages@meta.data)
colData(sca)$Replicate <- substr(rownames(colData(sca)),10,11)
C1 <- freq(sca[ ,sca$Cluster == 1]) >= 0.1
C2 <- freq(sca[ ,sca$Cluster == 2]) >= 0.1
C3 <- freq(sca[ ,sca$Cluster == 3]) >= 0.1
C4 <- freq(sca[ ,sca$Cluster == 4]) >= 0.1
C5 <- freq(sca[ ,sca$Cluster == 5]) >= 0.1
C6 <- freq(sca[ ,sca$Cluster == 6]) >= 0.1
expressed_genes <- sort(unique(c(names(which(C1 == TRUE)),names(which(C2 == TRUE)),names(which(C3 == TRUE)),names(which(C4 == TRUE)),names(which(C5 == TRUE)),names(which(C6 == TRUE)))))

## Setup results for export to MacSpectrum (and convert to ENSEMBL IDs)
Convert1 <- as.data.frame(org.Mm.egSYMBOL)
Convert2 <- as.data.frame(org.Mm.egENSEMBL)
Convert <- merge(Convert1, Convert2, by="gene_id")
MacroSC <- Macrophages@assays$RNA@data
MacroSC <- as.data.frame(MacroSC[ rownames(MacroSC) %in% expressed_genes,])
MacroSC$Gene <- rownames(MacroSC)
MacroSC <- merge(MacroSC, Convert[,c("ensembl_id","symbol")], by.x="Gene", by.y="symbol")
MacroSC <- MacroSC[ duplicated(MacroSC$Gene) == F,]
MacroSC <- MacroSC[ duplicated(MacroSC$ensembl_id) == F,]
MacroSC <- MacroSC[,c(ncol(MacroSC),2:(ncol(MacroSC)-1))]
write.table(MacroSC, "Macro_SC.csv", sep=",", quote=F, col.names=F, row.names=F) # Submit this file to the MacSpectrum server
Features <- data.frame(Cell = colnames(MacroSC)[2:ncol(MacroSC)])
Features <- merge(Features, Macrophages@meta.data[,c("Diet","Replicate")], by.x="Cell", by.y=0, sort=F)
write.table(Features$Tissue_Diet, "Macro_SC_Features.csv", sep=",", col.names=F, row.names=F, quote=F) # Submit this file to the MacSpectrum server

## Process the results
# Read data
MPI <- read.csv(file="MPI_AMDI_table.csv")
MD <- Macrophages@meta.data
MD$Cell <- rownames(MD)
Features$MPI<- MPI$MPI
Features$AMDI <- MPI$AMDI
MD <- merge(MD, Features[,c("Cell","MPI","AMDI")], by="Cell", sort=F)
rownames(MD) <- MD$Cell
Macrophages@meta.data <- MD

# Capture stats (estimated median and confidence intervals) with wilcox.test
Stats <- as.data.frame(matrix(ncol=6, nrow=6))
for (i in 1:6) {
  Stats[i,1] <- wilcox.test(MD[ MD$Cluster == i, "AMDI"], conf.int=T)$estimate
  Stats[i,2] <- wilcox.test(MD[ MD$Cluster == i, "AMDI"], conf.int=T)$conf.int[1]
  Stats[i,3] <- wilcox.test(MD[ MD$Cluster == i, "AMDI"], conf.int=T)$conf.int[2]
  Stats[i,4] <- wilcox.test(MD[ MD$Cluster == i, "MPI"], conf.int=T)$estimate
  Stats[i,5] <- wilcox.test(MD[ MD$Cluster == i, "MPI"], conf.int=T)$conf.int[1]
  Stats[i,6] <- wilcox.test(MD[ MD$Cluster == i, "MPI"], conf.int=T)$conf.int[2]
}
  
# Plot the result
plot(Stats[,1], Stats[,4], xlab="AMDI", ylab="MPI", xlim=c(-18,6), pch=16, ylim=c(-5,5))
abline(h=0)
abline(v=0)
for (i in 1:6) { arrows(Stats[i,2], Stats[i,4],Stats[i,3],Stats[i,4], code=3, length = 0.05, angle=90) }
for (i in 1:6) { arrows(Stats[i,1], Stats[i,5],Stats[i,1],Stats[i,6], code=3, length = 0.05, angle=90) }

### Diet-dependent gene regulation 
## Heatmap in PVM and NPVMs
# Extract diet-regulated genes from PVMs and NPVMs
NPVM_Diet <- Result[[1]]
NPVM_Diet <- NPVM_Diet[NPVM_Diet$Diet != "NS", ]
PVM_Diet <- Result[[3]]
PVM_Diet <- PVM_Diet[PVM_Diet$Diet != "NS", ]

# Calculate average expression for diet-regulated genes
Averages <- as.data.frame(matrix(ncol = 1, nrow = length(unique(sort(c(NPVM_Diet[, 1], PVM_Diet[, 1]))))))
Averages[, 1] <- unique(sort(c(NPVM_Diet[, 1], PVM_Diet[, 1])))
for (i in 1:nrow(Averages)) {
  Averages[i, 2] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "1" & Immune@meta.data$Diet == "LFD" & Immune@meta.data$Replicate == "R1", ])])))
  Averages[i, 3] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "1" & Immune@meta.data$Diet == "LFD" & Immune@meta.data$Replicate == "R2", ])])))
  Averages[i, 4] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "1" & Immune@meta.data$Diet == "HFD" & Immune@meta.data$Replicate == "R1", ])])))
  Averages[i, 5] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "1" & Immune@meta.data$Diet == "HFD" & Immune@meta.data$Replicate == "R2", ])])))
  Averages[i, 6] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "3" & Immune@meta.data$Diet == "LFD" & Immune@meta.data$Replicate == "R1", ])])))
  Averages[i, 7] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "3" & Immune@meta.data$Diet == "LFD" & Immune@meta.data$Replicate == "R2", ])])))
  Averages[i, 8] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "3" & Immune@meta.data$Diet == "HFD" & Immune@meta.data$Replicate == "R1", ])])))
  Averages[i, 9] <- log1p(mean(expm1(Immune@assays$RNA@data[rownames(Immune@assays$RNA@data) %in% Averages[i, 1], colnames(Immune@assays$RNA@data) %in% rownames(Immune@meta.data[Immune@meta.data$Label == "3" & Immune@meta.data$Diet == "HFD" & Immune@meta.data$Replicate == "R2", ])])))
}
rownames(Averages) <- Averages[, 1]

# Scale each experiment
Exp1 <- Averages[, c(2, 4, 6, 8)]
Exp1 <- t(scale(t(Exp1)))
Exp2 <- Averages[, c(3, 5, 7, 9)]
Exp2 <- t(scale(t(Exp2)))
Averages[, c(2, 4, 6, 8)] <- Exp1
Averages[, c(3, 5, 7, 9)] <- Exp2

# Cluster the data and highlight some genes (set seed to make it reproducible)
ToHighlight <- c("Pparg","Lpl","Cd74","Dpep2","Cd84","Cd63","Cd9","Rel","Adgre1","Cd300a","Ccr1","Ccr2","Fn1","Lyn","Retnla","Acaca","Lifr","Insr","Pten","Fasn","Cd163","Csf1r","Nrp1","Il16","Cd209d")
set.seed(100)
HTM <- Heatmap(Averages[,c(2:9)], cluster_columns=F, km=5, clustering_distance_rows = "manhattan", row_km_repeats = 100, show_row_names=F, right_annotation = rowAnnotation(foo = anno_mark(at = which(rownames(Averages) %in% ToHighlight), labels = rownames(Averages)[which(rownames(Averages) %in% ToHighlight)])))

# Extract the cluster order and save it 
set.seed(100)
Clusters <- row_order(HTM)
saveRDS(Clusters, "eWAT_Macrophage_Diet_Clusters.Rds")

# Plot the heatmap
set.seed(100)
HTM

## Overlap between LAM marker genes and diet-regulated cluster 1
# Get groups of marker genes
C1_Markers <- rownames(Averages[ Clusters[[1]],])
LAM_Markers <- Result[[2]]
LAM_Markers <- LAM_Markers[ LAM_Markers$Marker != "NS" ,]

# Calculate the overlap
Overlap <- sum(LAM_Markers$Symbol %in% C1_Markers)/nrow(LAM_Markers)

# Get list of tested genes in PVM and NPVM
Genes <- unique(sort(c(Result[[1]]$Symbol, Result[[3]]$Symbol)))

# Randomized selection
Background <- data.frame()
for (i in 1:100) {
  Random <- Genes[sample(1:length(Genes), length(C1_Markers))]
  Background[i,1] <- sum(LAM_Markers$Symbol %in% Random)/nrow(LAM_Markers)
  }

# Plot
Blt <- barplot(c(Overlap, mean(Background[,1])), las=1, ylim=c(0, 0.35), names = c("Cluster 1", "Random"))
arrows(Blt[2],mean(Background[,1])+(sd(Background[,1])/sqrt(100)),Blt[2],mean(Background[,1])-(sd(Background[,1])/sqrt(100)), length = 0.05, angle = 90, code = 3)

### Cytokine expression
# Extracting gene names of cytokines using biomaRt
res <- getBM(c("external_gene_name"), filters = "go", values = "GO:0005125", mart = ensembl)
Cytokines <- res[, 1]

# Extracting mean cytokine expression values across replicates and diets in macrophages
Averages <- data.frame(gene = Cytokines)
for (i in 1:nrow(Averages)) {
  Averages[i, 2] <- log1p(mean(expm1(Macrophages@assays$RNA@data[rownames(Macrophages@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[Macrophages@meta.data$Replicate == "R1" & Macrophages@meta.data$Diet == "LFD", ])])))
  Averages[i, 3] <- log1p(mean(expm1(Macrophages@assays$RNA@data[rownames(Macrophages@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[Macrophages@meta.data$Replicate == "R1" & Macrophages@meta.data$Diet == "HFD", ])])))
  Averages[i, 4] <- log1p(mean(expm1(Macrophages@assays$RNA@data[rownames(Macrophages@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[Macrophages@meta.data$Replicate == "R2" & Macrophages@meta.data$Diet == "LFD", ])])))
  Averages[i, 5] <- log1p(mean(expm1(Macrophages@assays$RNA@data[rownames(Macrophages@assays$RNA@data) %in% Averages[i, 1], colnames(Macrophages@assays$RNA@data) %in% rownames(Macrophages@meta.data[Macrophages@meta.data$Replicate == "R2" & Macrophages@meta.data$Diet == "HFD", ])])))
}
rownames(Averages) <- Averages[, 1]
names(Averages)[2:5] <- c("LFD_R1", "HFD_R1", "LFD_R2", "HFD_R2")

# Remove genes without signal
Averages <- Averages[!is.na(Averages[, 2]), ]

# Add diet logFC for each replicate
Averages[, "logFC_R1"] <- Averages[, 3] - Averages[, 2]
Averages[, "logFC_R2"] <- Averages[, 5] - Averages[, 4]

# Add mean logFC and SEM
Averages[, "Mean_logFC"] <- apply(Averages[, c(6,7)], 1, FUN = "mean")
Averages[, "sem"] <- apply(Averages[, c(6,7)], 1, FUN = "sd") / sqrt(2)

# Add maximum expression values across diets
Averages[, "Max"] <- apply(cbind(rowMeans(Averages[, c(2, 4)]), rowMeans(Averages[, c(3, 5)])), 1, FUN = "max")

# Plot logFC against maximum expression
plot(Averages[, 10], Averages[, 8], xlab = "Maximum log expression", ylab = "logFC", pch = 16, las = 1, main = "Cytokines")
abline(h = 0)
```
[Back to start](../README.md)<br>