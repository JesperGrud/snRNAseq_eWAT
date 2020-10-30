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
FAP <- subset(eWAT, subset = Annotation == "FAP")

# Embedding and clustering of FAPs into subpopulations  
# NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN. 
# Thus, in order to reproduce downstream analyses using the same Harmony results, please download eWAT_FAP.Rds
FAP <- FindVariableFeatures(FAP, nfeatures=2000)
VariableFeatures(FAP) <- VariableFeatures(FAP)[!(VariableFeatures(FAP) %in% Amb)]
FAP <- ScaleData(FAP)
FAP <- RunPCA(FAP)
FAP <- RunHarmony(FAP, group.by.vars="Dataset")
FAP <- RunUMAP(FAP, dims=1:20, reduction="harmony")
FAP <- FindNeighbors(object = FAP, dims = 1:2, reduction = "umap", k.param=50)
FAP <- FindClusters(FAP, resolution = 0.1, algorithm = 1)
FAP@project.name <- "FAP"
# FAP <- readRDS("eWAT_FAP.Rds") 

## Adjust cluster labels from 0-indexed to 1-indexed
FAP$Label <- 0
FAP@meta.data[ FAP@meta.data$seurat_clusters == 0,"Label"] <- 2
FAP@meta.data[ FAP@meta.data$seurat_clusters == 1,"Label"] <- 3
FAP@meta.data[ FAP@meta.data$seurat_clusters == 2,"Label"] <- 4
FAP@meta.data[ FAP@meta.data$seurat_clusters == 3,"Label"] <- 1
FAP <- SetIdent(FAP, value=FAP$Label)

## Detection of marker genes and dietary effects
# Setup
Data <- FAP
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
saveRDS(Result, "eWAT_FAP_Markers.Rds")

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
FAP$Subtype <- "NA"
FAP@meta.data[FAP@meta.data$Label %in% 1, "Subtype"] <- "FAP1"
FAP@meta.data[FAP@meta.data$Label %in% 2, "Subtype"] <- "FAP2"
FAP@meta.data[FAP@meta.data$Label %in% 3, "Subtype"] <- "FAP3"
FAP@meta.data[FAP@meta.data$Label %in% 4, "Subtype"] <- "FAP4"
FAP <- SetIdent(FAP, value = FAP$Subtype)

# Transfer labels to the eWAT object
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(FAP@meta.data[ FAP@meta.data$Subtype == "FAP1",]),"Subtype"] <- "FAP1"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(FAP@meta.data[ FAP@meta.data$Subtype == "FAP2",]),"Subtype"] <- "FAP2"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(FAP@meta.data[ FAP@meta.data$Subtype == "FAP3",]),"Subtype"] <- "FAP3"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(FAP@meta.data[ FAP@meta.data$Subtype == "FAP4",]),"Subtype"] <- "FAP4"

# Save results
saveRDS(FAP, "eWAT_FAP.Rds") ## This object is downloadable from Open Science Framework
saveRDS(eWAT, "eWAT_Annotated.Rds") ## This object is downloadable from Open Science Framework with additional annotations added by subclustering each major cell type

## DimPlots
DimPlot(FAP, label = T)
DimPlot(FAP, label = T, split.by="Diet")
DimPlot(FAP, label = T, split.by="Replicate")

## Marker gene heatmap
# Select marker genes for each cell type
FAP1 <- data.frame(gene = c("Hes1", "Lox", "Igf1", "Lbp", "Foxp2"), type = "FAP1")
FAP2 <- data.frame(gene = c("Cd36", "Lpl", "Gata4","Pparg","Fgf10"), type = "FAP2")
FAP3 <- data.frame(gene = c("Ebf1", "Rock2", "Zfp521","Bmp6","Ebf2"), type = "FAP3")
FAP4 <- data.frame(gene = c("Fbn1", "Klf4", "Klf2","Fn1","Loxl1"), type = "FAP4")

# Combine results
Averages <- rbind(FAP1, FAP2, FAP3, FAP4)

# Calculate average expression for selected marker genes
for (i in 1:nrow(Averages)) {
  Averages[i,3] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Averages[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "1" & FAP@meta.data$Replicate == "R1",])])))
  Averages[i,4] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Averages[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "1" & FAP@meta.data$Replicate == "R2",])])))
  Averages[i,5] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Averages[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "2" & FAP@meta.data$Replicate == "R1",])])))
  Averages[i,6] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Averages[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "2" & FAP@meta.data$Replicate == "R2",])])))
  Averages[i,7] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Averages[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "3" & FAP@meta.data$Replicate == "R1",])])))
  Averages[i,8] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Averages[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "3" & FAP@meta.data$Replicate == "R2",])])))
  Averages[i,9] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Averages[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "4" & FAP@meta.data$Replicate == "R1",])])))
  Averages[i,10] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Averages[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "4" & FAP@meta.data$Replicate == "R2",])])))
}
rownames(Averages) <- Averages[,1]

# Scale each experiment
Exp1 <- Averages[,c(3,5,7,9)]
Exp1 <- t(scale(t(Exp1)))
Exp2 <- Averages[,c(4,6,8,10)]
Exp2 <- t(scale(t(Exp2)))
Averages[,c(3,5,7,9)] <- Exp1
Averages[,c(4,6,8,10)] <- Exp2

# Plot the heatmap
Heatmap(Averages[,c(3:10)], cluster_columns=F)

## Gene module scoring
# Extract marker gene results
Clust1 <- Result[[1]] 
Clust2 <- Result[[2]] 
Clust3 <- Result[[3]] 
Clust4 <- Result[[4]] 

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

# Computing gene module scores using the top 50 most specific marker genes
FAP <- AddModuleScore(FAP, features = list(Cluster1 = Clust1[order(-Clust1$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster2 = Clust2[order(-Clust2$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster3 = Clust3[order(-Clust3$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster4 = Clust4[order(-Clust4$logFC_OvO), ][1:50, ]$Symbol))

# Scale the module scores (across all nuclei)
FAP$Cluster1 <- scale(FAP$Cluster1)
FAP$Cluster2 <- scale(FAP$Cluster2)
FAP$Cluster3 <- scale(FAP$Cluster3)
FAP$Cluster4 <- scale(FAP$Cluster4)

# Computing mean scaled module scores in each cluster in each replicate
Averages <- as.data.frame(matrix(ncol=length(Result), nrow=length(Result)*2))

# Set column and row names
Labels <- FAP@meta.data[, c("Subtype","Label")]
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
      Averages[Counter,module] <- mean(FAP@meta.data[ FAP@meta.data$Label == label & FAP$Replicate == rep,paste("Cluster",module,sep="")])
      Counter <- Counter + 1
    }
  }
}

## Plot mean scaled module scores
col <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Heatmap(as.matrix(Averages), cluster_columns = F, cluster_rows = F, col = col)

## Subtype proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(FAP$Subtype)), nrow=4))
colnames(Composition) <- unique(FAP$Subtype)
for (i in 1:length(unique(FAP$Subtype))) {
  Composition[1,i] <- nrow(FAP@meta.data[ FAP@meta.data$Dataset == "LFD_R1" & FAP@meta.data$Subtype == colnames(Composition)[i],]) / nrow(FAP@meta.data[ FAP@meta.data$Dataset == "LFD_R1",])
  Composition[2,i] <- nrow(FAP@meta.data[ FAP@meta.data$Dataset == "LFD_R2" & FAP@meta.data$Subtype == colnames(Composition)[i],]) / nrow(FAP@meta.data[ FAP@meta.data$Dataset == "LFD_R2",])
  Composition[3,i] <- nrow(FAP@meta.data[ FAP@meta.data$Dataset == "HFD_R1" & FAP@meta.data$Subtype == colnames(Composition)[i],]) / nrow(FAP@meta.data[ FAP@meta.data$Dataset == "HFD_R1",])
  Composition[4,i] <- nrow(FAP@meta.data[ FAP@meta.data$Dataset == "HFD_R2" & FAP@meta.data$Subtype == colnames(Composition)[i],]) / nrow(FAP@meta.data[ FAP@meta.data$Dataset == "HFD_R2",])
}

# Calculate stats (mean, SEM, P)
Stats <- data.frame(matrix(ncol=5, nrow=length(unique(FAP$Subtype))))
colnames(Stats) <- c("LFD_mean","LFD_SEM","HFD_mean","HFD_SEM","Pvalue")
rownames(Stats) <- colnames(Composition)
for (i in 1:length(unique(FAP$Subtype))) {
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

### Pre-adipocyte marker genes
# Define data.frame with markers
Markers <- data.frame(Gene = c("Ly6a","Pdgfra","Pdgfrb","Myh11","Acta2","Cd24a","Cd34","Zfp423"))

# In each subpopulation in each replicate, calculate average expression
for (i in 1:nrow(Markers)) {
  Markers[i,2] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Markers[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "1" & FAP@meta.data$Replicate == "R1",])])))
  Markers[i,3] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Markers[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "1" & FAP@meta.data$Replicate == "R2",])])))
  Markers[i,4] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Markers[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "2" & FAP@meta.data$Replicate == "R1",])])))
  Markers[i,5] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Markers[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "2" & FAP@meta.data$Replicate == "R2",])])))
  Markers[i,6] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Markers[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "3" & FAP@meta.data$Replicate == "R1",])])))
  Markers[i,7] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Markers[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "3" & FAP@meta.data$Replicate == "R2",])])))
  Markers[i,8] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Markers[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "4" & FAP@meta.data$Replicate == "R1",])])))
  Markers[i,9] <- log1p(mean(expm1(FAP@assays$RNA@data[ rownames(FAP@assays$RNA@data) %in% Markers[i,1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "4" & FAP@meta.data$Replicate == "R2",])])))
}

# Calculate mean and SEM across replicates for each gene in each dietary condition
Stats <- data.frame()
for (i in 1:nrow(Markers)) {
  Stats[i,1] <- mean(as.numeric(Markers[i,c(2,3)]))
  Stats[i,2] <- sd(as.numeric(Markers[i,c(2,3)]))/sqrt(2)
  Stats[i,3] <- mean(as.numeric(Markers[i,c(4,5)]))
  Stats[i,4] <- sd(as.numeric(Markers[i,c(4,5)]))/sqrt(2)
  Stats[i,5] <- mean(as.numeric(Markers[i,c(6,7)]))
  Stats[i,6] <- sd(as.numeric(Markers[i,c(6,7)]))/sqrt(2)
  Stats[i,7] <- mean(as.numeric(Markers[i,c(8,9)]))
  Stats[i,8] <- sd(as.numeric(Markers[i,c(8,9)]))/sqrt(2)
  rownames(Stats)[i] <- as.character(Markers[i,1])
}

# Barplot in a loop
par(mfcol=c(2,4))
for (i in 1:nrow(Markers)) {
  # Find ylimit
  ylim <- max(c(Stats[i,1]+Stats[i,2],Stats[i,3]+Stats[i,4],Stats[i,5]+Stats[i,6], Stats[i,7]+Stats[i,8]))
  ylim <- ylim * 1.2
  
  # Barplot
  Blt <- barplot(c(Stats[i,1],Stats[i,3], Stats[i,5], Stats[i,7]), las=1, ylim=c(0,ylim), col = c("purple", "orange","green","blue"), main = rownames(Stats)[i])
  
  # Arrows in another loop
  for (m in 1:4) { arrows(Blt[m], Stats[i,(((m-1)*2) + 1)]+Stats[i,(((m-1)*2) + 2)], Blt[m], Stats[i,(((m-1)*2) + 1)]-Stats[i,(((m-1)*2) + 2)], length = 0.05, angle = 90, code = 3 )}
}

### Comparison with other studies (Granneman lab: Burl et al., Deplanke lab: Schwalie et al., Gupta lab: Hepler et al., Seale lab: Merrick et al.)
## Disregard mesothelial cells from Hepler et al., as our mesothelial cells do not cluster with FAPs
# Import all gene lists and combine (Extracted from the papers, main/supplementary figures or tables)
Granneman_ASC1 <- read.delim("Granneman_ASC1.txt",header=F)
Granneman_ASC2 <- read.delim("Granneman_ASC2.txt",header=F)
Deplancke_P1 <- read.delim("Deplancke_P1.txt",header=F)
Deplancke_P2 <- read.delim("Deplancke_P2.txt",header=F)
Deplancke_P3 <- read.delim("Deplancke_P3.txt",header=F)
Gupta_APC <- read.delim("Gupta_APC.txt",header=F)
Gupta_FIP <- read.delim("/Gupta_FIP.txt",header=F)
Gupta_CPA <- read.delim("Gupta_CPA.txt",header=F)
Seale_DPP4 <- read.delim("Seale_DPP4.txt",header=F)
Seale_ICAM1 <- read.delim("Seale_ICAM1.txt",header=F)
Seale_CD142 <- read.delim("Seale_CD142.txt",header=F)
GeneLists <- list(as.character(Granneman_ASC1[,1]), as.character(Granneman_ASC2[,1]), as.character(Deplancke_P1[,1]), as.character(Deplancke_P2[,1]), as.character(Deplancke_P3[,1]), as.character(Gupta_FIP[,1]), as.character(Gupta_APC[,1]), as.character(Gupta_CPA[,1]), as.character(Seale_DPP4[,1]), as.character(Seale_ICAM1[,1]), as.character(Seale_CD142[,1]))
names(GeneLists) <- c("Granneman_ASC1", "Granneman_ASC2","Deplancke_P1","Deplancke_P2","Deplancke_P3", "Gupta_FIP","Gupta_APC","Gupta_CPA","Seale_DPP4","Seale_ICAM1","Seale_CD142")

# Calculate module scores. NOTE: THIS IS NON-DETERMINISTIC. THE RESULTS VARY SLIGHTLY FROM RUN-TO-RUN.
FAP <- AddModuleScore(FAP, features = GeneLists, name = "Comparison")

# Fix the names of the module scores
colnames(FAP@meta.data)[50:60] <- names(GeneLists)

# Scale the module scores
for (i in 50:60) { FAP@meta.data[,i] <- scale(FAP@meta.data[,i]) }

## Calculate label (binary classification)
# Granneman lab
FAP$Granneman <- "None"
FAP@meta.data[ FAP@meta.data$Granneman_ASC1 > FAP@meta.data$Granneman_ASC2, "Granneman"] <- "ASC1"
FAP@meta.data[ FAP@meta.data$Granneman_ASC2 > FAP@meta.data$Granneman_ASC1, "Granneman"] <- "ASC2"

# Deplancke lab
FAP$Deplancke <- "None"
FAP@meta.data[ FAP@meta.data$Deplancke_P1 > FAP@meta.data$Deplancke_P2 & FAP@meta.data$Deplancke_P1 > FAP@meta.data$Deplancke_P3, "Deplancke"] <- "P1"
FAP@meta.data[ FAP@meta.data$Deplancke_P2 > FAP@meta.data$Deplancke_P1 & FAP@meta.data$Deplancke_P2 > FAP@meta.data$Deplancke_P3, "Deplancke"] <- "P2"
FAP@meta.data[ FAP@meta.data$Deplancke_P3 > FAP@meta.data$Deplancke_P2 & FAP@meta.data$Deplancke_P3 > FAP@meta.data$Deplancke_P1, "Deplancke"] <- "P3"

# Gupta lab
FAP$Gupta <- "None"
FAP@meta.data[ FAP@meta.data$Gupta_FIP > FAP@meta.data$Gupta_APC & FAP@meta.data$Gupta_FIP > FAP@meta.data$Gupta_CPA, "Gupta"] <- "FIP"
FAP@meta.data[ FAP@meta.data$Gupta_APC > FAP@meta.data$Gupta_FIP & FAP@meta.data$Gupta_APC > FAP@meta.data$Gupta_CPA, "Gupta"] <- "APC"
FAP@meta.data[ FAP@meta.data$Gupta_CPA > FAP@meta.data$Gupta_FIP & FAP@meta.data$Gupta_CPA > FAP@meta.data$Gupta_APC, "Gupta"] <- "CPA"

# Seale lab
FAP$Seale <- "None"
FAP@meta.data[ FAP@meta.data$Seale_DPP4 > FAP@meta.data$Seale_ICAM1 & FAP@meta.data$Seale_DPP4 > FAP@meta.data$Seale_CD142, "Seale"] <- "DPP4"
FAP@meta.data[ FAP@meta.data$Seale_ICAM1 > FAP@meta.data$Seale_DPP4 & FAP@meta.data$Seale_ICAM1 > FAP@meta.data$Seale_CD142, "Seale"] <- "ICAM1"
FAP@meta.data[ FAP@meta.data$Seale_CD142 > FAP@meta.data$Seale_DPP4 & FAP@meta.data$Seale_CD142 > FAP@meta.data$Seale_ICAM1, "Seale"] <- "CD142"

## Barplots of fractions (using binary labels)
barplot(t(as.matrix(table(FAP$Subtype, FAP$Granneman))/rowSums(as.matrix(table(FAP$Subtype, FAP$Granneman)))), beside=T, las=1, ylab="Fraction of cluster", names = c("1","2","3","4"), ylim=c(0,1), main="Burl et al.")
barplot(t(as.matrix(table(FAP$Subtype, FAP$Deplancke))/rowSums(as.matrix(table(FAP$Subtype, FAP$Deplancke)))), beside=T, las=1, ylab="Fraction of cluster", names = c("1","2","3","4"), ylim=c(0,1), main="Schwalie et al.")
barplot(t(as.matrix(table(FAP$Subtype, FAP$Gupta))/rowSums(as.matrix(table(FAP$Subtype, FAP$Gupta)))), beside=T, las=1, ylab="Fraction of cluster", names = c("1","2","3","4"), ylim=c(0,1), main="Hepler et al.")
barplot(t(as.matrix(table(FAP$Subtype, FAP$Seale))/rowSums(as.matrix(table(FAP$Subtype, FAP$Seale)))), beside=T, las=1, ylab="Fraction of cluster", names = c("1","2","3","4"), ylim=c(0,1), main="Merrick et al.")

## DimPlots (binary classifications)
FAP <- SetIdent(FAP, value = FAP$Granneman)
DimPlot(FAP, label=T)
FAP <- SetIdent(FAP, value = FAP$Deplancke)
DimPlot(FAP, label=T)
FAP <- SetIdent(FAP, value = FAP$Gupta)
DimPlot(FAP, label=T)
FAP <- SetIdent(FAP, value = FAP$Seale)
DimPlot(FAP, label=T)
FAP <- SetIdent(FAP, value = FAP$Subtype)

## FeaturePlots (continuous scores)
FeaturePlot(FAP, "Granneman_ASC1")
FeaturePlot(FAP, "Granneman_ASC2")
FeaturePlot(FAP, "Deplancke_P1")
FeaturePlot(FAP, "Deplancke_P2")
FeaturePlot(FAP, "Deplancke_P3")
FeaturePlot(FAP, "Gupta_FIP")
FeaturePlot(FAP, "Gupta_APC")
FeaturePlot(FAP, "Gupta_CPA")
FeaturePlot(FAP, "Seale_DPP4")
FeaturePlot(FAP, "Seale_ICAM1")
FeaturePlot(FAP, "Seale_CD142")

## Boxplots (continuous scores)
boxplot(
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Granneman_ASC1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Granneman_ASC2"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Granneman_ASC1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Granneman_ASC2"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Granneman_ASC1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Granneman_ASC2"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Granneman_ASC1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Granneman_ASC2"][,1], outline=F, main = "Burl et al.", las=1)

boxplot(
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Deplancke_P1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Deplancke_P2"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Deplancke_P3"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Deplancke_P1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Deplancke_P2"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Deplancke_P3"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Deplancke_P1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Deplancke_P2"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Deplancke_P3"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Deplancke_P1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Deplancke_P2"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Deplancke_P3"][,1], outline=F, main = "Schwalie et al.", las=1)

boxplot(
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Gupta_FIP"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Gupta_APC"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Gupta_CPA"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Gupta_FIP"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Gupta_APC"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Gupta_CPA"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Gupta_FIP"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Gupta_APC"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Gupta_CPA"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Gupta_FIP"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Gupta_APC"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Gupta_CPA"][,1], outline=F, main = "Hepler et al.", las=1)

boxplot(
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Seale_DPP4"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Seale_ICAM1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP1", "Seale_CD142"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Seale_DPP4"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Seale_ICAM1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP2", "Seale_CD142"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Seale_DPP4"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Seale_ICAM1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP3", "Seale_CD142"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Seale_DPP4"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Seale_ICAM1"][,1],
  FAP@meta.data[ FAP@meta.data$Subtype == "FAP4", "Seale_CD142"][,1], outline=F, main = "Merrick et al.", las=1)

### Diet-dependent regulation of gene expression in FAP subtypes
# Extract diet-regulated genes from all FAP clusters
Cluster0 <- Result[[1]]
Cluster0 <- Cluster0[Cluster0$Diet != "NS", ]
Cluster1 <- Result[[2]]
Cluster1 <- Cluster1[Cluster1$Diet != "NS", ]
Cluster2 <- Result[[3]]
Cluster2 <- Cluster2[Cluster2$Diet != "NS", ]
Cluster3 <- Result[[4]]
Cluster3 <- Cluster3[Cluster3$Diet != "NS", ]

# Calculate average expression for diet-regulated genes
Averages <- as.data.frame(matrix(ncol = 1, nrow = length(unique(sort(c(Cluster0[, 1], Cluster1[, 1], Cluster2[, 1], Cluster3[, 1]))))))
Averages[, 1] <- unique(sort(c(Cluster0[, 1], Cluster1[, 1], Cluster2[, 1], Cluster3[, 1])))
for (i in 1:nrow(Averages)) {
  Averages[i, 2] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "1" & FAP@meta.data$Replicate == "R1" & FAP@meta.data$Diet == "LFD", ])])))
  Averages[i, 3] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "1" & FAP@meta.data$Replicate == "R1" & FAP@meta.data$Diet == "HFD", ])])))
  Averages[i, 4] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "2" & FAP@meta.data$Replicate == "R1" & FAP@meta.data$Diet == "LFD", ])])))
  Averages[i, 5] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "2" & FAP@meta.data$Replicate == "R1" & FAP@meta.data$Diet == "HFD", ])])))
  Averages[i, 6] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "3" & FAP@meta.data$Replicate == "R1" & FAP@meta.data$Diet == "LFD", ])])))
  Averages[i, 7] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "3" & FAP@meta.data$Replicate == "R1" & FAP@meta.data$Diet == "HFD", ])])))
  Averages[i, 8] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "4" & FAP@meta.data$Replicate == "R1" & FAP@meta.data$Diet == "LFD", ])])))
  Averages[i, 9] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "4" & FAP@meta.data$Replicate == "R1" & FAP@meta.data$Diet == "HFD", ])])))
  Averages[i, 10] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "1" & FAP@meta.data$Replicate == "R2" & FAP@meta.data$Diet == "LFD", ])])))
  Averages[i, 11] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "1" & FAP@meta.data$Replicate == "R2" & FAP@meta.data$Diet == "HFD", ])])))
  Averages[i, 12] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "2" & FAP@meta.data$Replicate == "R2" & FAP@meta.data$Diet == "LFD", ])])))
  Averages[i, 13] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "2" & FAP@meta.data$Replicate == "R2" & FAP@meta.data$Diet == "HFD", ])])))
  Averages[i, 14] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "3" & FAP@meta.data$Replicate == "R2" & FAP@meta.data$Diet == "LFD", ])])))
  Averages[i, 15] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "3" & FAP@meta.data$Replicate == "R2" & FAP@meta.data$Diet == "HFD", ])])))
  Averages[i, 16] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "4" & FAP@meta.data$Replicate == "R2" & FAP@meta.data$Diet == "LFD", ])])))
  Averages[i, 17] <- log1p(mean(expm1(FAP@assays$RNA@data[rownames(FAP@assays$RNA@data) %in% Averages[i, 1], colnames(FAP@assays$RNA@data) %in% rownames(FAP@meta.data[ FAP@meta.data$Label == "4" & FAP@meta.data$Replicate == "R2" & FAP@meta.data$Diet == "HFD", ])])))
}
rownames(Averages) <- Averages[, 1]

# Scale each experiment separately
Averages[, 2:9] <- t(scale(t(Averages[, 2:9])))
Averages[, 10:17] <- t(scale(t(Averages[, 10:17])))

# Cluster the data and highlight some genes (set seed to make it reproducible)
ToHighlight <- c("Ccl8", "Ccl9", "Ctss", "Ctsd", "Ar", "Pdgfra", "Id3", "Sox9", "Acaca", "Fasn", "Srebf2", "Insig1", "Smurf2", "Tgfbr2", "Tgfbr3", "Ltpb1", "Rps5", "Rpl4", "Eif1", "Eif5a", "Rps6", "Rpl36", "Rps1", "Rpl27", "Col1a1", "Col3a1", "Col4a2", "Col8a1")
set.seed(100)
HTM <- ComplexHeatmap::Heatmap(Averages[, c(2,10,3,11,4,12,5,13,6,14,7,15,8,16,9,17)], cluster_columns=F, km = 7, show_row_names = F, clustering_distance_rows = "manhattan", row_km_repeats = 100, right_annotation = rowAnnotation(foo = anno_mark(at = which(rownames(Averages) %in% ToHighlight), labels = rownames(Averages)[which(rownames(Averages) %in% ToHighlight)])))

# Extract the cluster order and save it
set.seed(100)
Clusters <- row_order(HTM)
saveRDS(Clusters, "eWAT_FAP_Diet_Clusters.Rds")

# Plot the heatmap
set.seed(100)
HTM
```
[Back to start](../README.md)<br>