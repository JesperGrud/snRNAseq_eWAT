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
Amb <- readRDS("Ambient.Rds")

## Adjust cluster labels from 0-indexed to 1-indexed
eWAT$Label <- 1
eWAT@meta.data[ eWAT@meta.data$seurat_clusters == "0","Label"] <- 1
eWAT@meta.data[ eWAT@meta.data$seurat_clusters == "1","Label"] <- 2
eWAT@meta.data[ eWAT@meta.data$seurat_clusters == "2","Label"] <- 3
eWAT@meta.data[ eWAT@meta.data$seurat_clusters == "3","Label"] <- 4
eWAT@meta.data[ eWAT@meta.data$seurat_clusters == "4","Label"] <- 5
eWAT@meta.data[ eWAT@meta.data$seurat_clusters == "5","Label"] <- 6
eWAT@meta.data[ eWAT@meta.data$seurat_clusters == "6","Label"] <- 7
eWAT@meta.data[ eWAT@meta.data$seurat_clusters == "7","Label"] <- 8
eWAT <- SetIdent(eWAT, value = eWAT$Label)

## Detection of marker genes and dietary effects
# Setup
Data <- eWAT
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
saveRDS(Result, "eWAT_Overall_Markers.Rds")

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
eWAT$Annotation <- "NA"
eWAT@meta.data[ eWAT@meta.data$Label %in% 1,"Annotation"] <- "Immune"
eWAT@meta.data[ eWAT@meta.data$Label %in% 2,"Annotation"] <- "FAP"
eWAT@meta.data[ eWAT@meta.data$Label %in% 3,"Annotation"] <- "Adipocyte"
eWAT@meta.data[ eWAT@meta.data$Label %in% 4,"Annotation"] <- "Mesothelial"
eWAT@meta.data[ eWAT@meta.data$Label %in% 5,"Annotation"] <- "Epididymal"
eWAT@meta.data[ eWAT@meta.data$Label %in% 6,"Annotation"] <- "Spermatozoa"
eWAT@meta.data[ eWAT@meta.data$Label %in% 7,"Annotation"] <- "Endothelial"
eWAT <- SetIdent(eWAT, value = eWAT$Annotation)

## Save the annotated object
saveRDS(eWAT, "eWAT_Annotated.Rds") ## This object is downloadable from Open Science Framework with additional annotations added by subclustering each major cell type

## DimPlots
DimPlot(eWAT, label=T)
DimPlot(eWAT, label=T, split.by="Diet")
DimPlot(eWAT, label=T, split.by="Replicate")

## FeaturePlots for the largest WAT intrinsic cell types
FeaturePlot(eWAT, features = c("Adgre1","Lipe","Pdgfra","Upk3b"), split.by="Replicate", max.cutoff = 3, pt.size = 0.001)

## Marker gene heatmap
# Select marker genes for each cell type
Immune <- data.frame(gene = c("Adgre1","Tec","Btk"), type = "Immune")
FAP <- data.frame(gene = c("Dcn","Col1a1","Pdgfra"), type = "FAP")
Adipocytes <- data.frame(gene = c("Pparg","Plin4","Lipe"), type = "Adipocytes")
Mesothelial <- data.frame(gene = c("Upk3b","Gpm6a","Msln"), type = "Mesothelial")
Epididymal <- data.frame(gene = c("Abcb5","Adam7","Ces5a"), type = "Epididymal") 
Spermatozoa <- data.frame(gene = c("Hydin", "Dnah12","Spef2"), type = "Spermatozoa")
Endothelial <- data.frame(gene = c("Vwf","Flt1","Esam"), type = "Endothelial")

# Combine results
Averages <- rbind(Immune,FAP,Adipocytes,Epididymal, Mesothelial, Spermatozoa,Endothelial)

# Calculate average expression for selected marker genes
for (i in 1:nrow(Averages)) {
  Averages[i,3] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "1" & eWAT@meta.data$Replicate == "R1",])])))
  Averages[i,4] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "1" & eWAT@meta.data$Replicate == "R2",])])))
  Averages[i,5] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "2" & eWAT@meta.data$Replicate == "R1",])])))
  Averages[i,6] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "2" & eWAT@meta.data$Replicate == "R2",])])))
  Averages[i,7] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "3" & eWAT@meta.data$Replicate == "R1",])])))
  Averages[i,8] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "3" & eWAT@meta.data$Replicate == "R2",])])))
  Averages[i,9] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "4" & eWAT@meta.data$Replicate == "R1",])])))
  Averages[i,10] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "4" & eWAT@meta.data$Replicate == "R2",])])))
  Averages[i,11] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "5" & eWAT@meta.data$Replicate == "R1",])])))
  Averages[i,12] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "5" & eWAT@meta.data$Replicate == "R2",])])))
  Averages[i,13] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "6" & eWAT@meta.data$Replicate == "R1",])])))
  Averages[i,14] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "6" & eWAT@meta.data$Replicate == "R2",])])))
  Averages[i,15] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "7" & eWAT@meta.data$Replicate == "R1",])])))
  Averages[i,16] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Label == "7" & eWAT@meta.data$Replicate == "R2",])])))
}
rownames(Averages) <- Averages[,1]

# Scale each experiment
Exp1 <- Averages[,c(3,5,7,9,11,13,15)]
Exp1 <- t(scale(t(Exp1)))
Exp2 <- Averages[,c(4,6,8,10,12,14,16)]
Exp2 <- t(scale(t(Exp2)))
Averages[,c(3,5,7,9,11,13,15)] <- Exp1
Averages[,c(4,6,8,10,12,14,16)] <- Exp2

# Plot the heatmap
Heatmap(Averages[,c(3:16)], cluster_columns=F)

## Gene module scoring
# Extract marker gene results
Clust1 <- Result[[1]] 
Clust2 <- Result[[2]] 
Clust3 <- Result[[3]] 
Clust4 <- Result[[4]] 
Clust5 <- Result[[5]] 
Clust6 <- Result[[6]] 
Clust7 <- Result[[7]]

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
Averages <- as.data.frame(matrix(ncol=length(Result), nrow=length(Result)*2))

# Set column and row names
Labels <- eWAT@meta.data[, c("Annotation","Label")]
Labels <- Labels[ duplicated(Labels$Label)==F,]
Labels <- Labels[ order(Labels$Label),]
colnames(Averages) <- paste(Labels$Annotation,"module",sep="_")
rownames(Averages)[seq(1,(nrow(Averages)-1),by=2)] <- paste(Labels$Annotation,"R1",sep="_")
rownames(Averages)[seq(2,nrow(Averages),by=2)] <- paste(Labels$Annotation,"R2",sep="_")

# Calculate averages
for (module in 1:length(Result)) {
  Counter <- 1
  for (label in 1:length(Result)) {
    for (rep in c("R1","R2")) {
      Averages[Counter,module] <- mean(eWAT@meta.data[ eWAT@meta.data$Label == label & eWAT$Replicate == rep,paste("Cluster",module,sep="")])
      Counter <- Counter + 1
    }
  }
}

## Plot mean scaled module scores
col <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(as.matrix(Averages), cluster_columns = F, cluster_rows = F, col = col)

## Cell type proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=7, nrow=4))
colnames(Composition) <- c("Immune","FAP","Adipocyte","Mesothelial","Epididymal","Spermatozoa","Endothelial")
for (i in 1:7) {
  Composition[1,i] <- nrow(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R1" & eWAT@meta.data$Annotation == colnames(Composition)[i],])
  Composition[2,i] <- nrow(eWAT@meta.data[ eWAT@meta.data$Dataset == "LFD_R2" & eWAT@meta.data$Annotation == colnames(Composition)[i],])
  Composition[3,i] <- nrow(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R1" & eWAT@meta.data$Annotation == colnames(Composition)[i],])
  Composition[4,i] <- nrow(eWAT@meta.data[ eWAT@meta.data$Dataset == "HFD_R2" & eWAT@meta.data$Annotation == colnames(Composition)[i],])
}
Composition <- Composition/rowSums(Composition)

# Calculate stats (mean, SEM, P)
Stats <- data.frame(matrix(ncol=5, nrow=7))
colnames(Stats) <- c("LFD_mean","LFD_SEM","HFD_mean","HFD_SEM","Pvalue")
rownames(Stats) <- colnames(Composition)
for (i in 1:7) {
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
```
[Back to start](../README.md)<br>