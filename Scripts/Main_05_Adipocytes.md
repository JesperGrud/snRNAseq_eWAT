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
Adipocytes <- subset(eWAT, subset = Annotation == "Adipocyte")

# Embedding and clustering of adipocytes into subpopulations  
# NOTE: THIS STEP IS NON-DETERMINISTIC (Harmony) - RESULTS VARY FROM RUN TO RUN. 
# Thus, in order to reproduce downstream analyses using the same Harmony results, please download eWAT_Adipocytes.Rds
Adipocytes <- FindVariableFeatures(Adipocytes, nfeatures=5000)
VariableFeatures(Adipocytes) <- VariableFeatures(Adipocytes)[!(VariableFeatures(Adipocytes) %in% Amb)]
Adipocytes <- ScaleData(Adipocytes)
Adipocytes <- RunPCA(Adipocytes)
Adipocytes <- RunHarmony(Adipocytes, group.by.vars="Dataset")
Adipocytes <- RunUMAP(Adipocytes, dims=1:10, reduction="harmony")
Adipocytes <- FindNeighbors(object = Adipocytes, dims = 1:10, reduction = "harmony")
Adipocytes <- FindClusters(Adipocytes, resolution = 0.2, algorithm = 1)
Adipocytes@project.name <- "Adipocytes"
# Adipocytes <- readRDS("eWAT_Adipocytes.Rds") 

## Adjust cluster labels from 0-indexed to 1-indexed
Adipocytes$Label <- 0
Adipocytes@meta.data[ Adipocytes@meta.data$seurat_clusters == 0,"Label"] <- 1
Adipocytes@meta.data[ Adipocytes@meta.data$seurat_clusters == 1,"Label"] <- 2
Adipocytes@meta.data[ Adipocytes@meta.data$seurat_clusters == 2,"Label"] <- 3
Adipocytes <- SetIdent(Adipocytes, value=Adipocytes$Label)

## Detection of marker genes and dietary effects
# Setup
Data <- Adipocytes
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
saveRDS(Result, "eWAT_Adipocyte_Markers.Rds")

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
Adipocytes$Subtype <- "NA"
Adipocytes@meta.data[Adipocytes@meta.data$Label %in% 1, "Subtype"] <- "LSA"
Adipocytes@meta.data[Adipocytes@meta.data$Label %in% 2, "Subtype"] <- "LGA"
Adipocytes@meta.data[Adipocytes@meta.data$Label %in% 3, "Subtype"] <- "SLSA"
Adipocytes <- SetIdent(Adipocytes, value = Adipocytes$Subtype)

# Transfer labels to the eWAT object
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Subtype == "LSA",]),"Subtype"] <- "LSA"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Subtype == "LGA",]),"Subtype"] <- "LGA"
eWAT@meta.data[ rownames(eWAT@meta.data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Subtype == "SLSA",]),"Subtype"] <- "SLSA"

# Save results
saveRDS(Adipocytes, "eWAT_Adipocytes.Rds") ## This object is downloadable from Open Science Framework
saveRDS(eWAT, "eWAT_Annotated.Rds") ## This object is downloadable from Open Science Framework with additional annotations added by subclustering each major cell type

## DimPlots
DimPlot(Adipocytes, label = T)
DimPlot(Adipocytes, label = T, split.by="Diet")
DimPlot(Adipocytes, label = T, split.by="Replicate")

## Adipocyte marker heatmap
Averages <- data.frame(gene = c("Pparg","Lipe","Plin4","Plin1","Lep","Cd36", "Adrb3", "Acsl1","Adipor2","Dgat1","Dgat2"))
for (i in 1:nrow(Averages)) {
  Averages[i,2] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation != "Adipocyte" & eWAT@meta.data$Replicate == "R1",])])))
  Averages[i,3] <- log1p(mean(expm1(eWAT@assays$RNA@data[ rownames(eWAT@assays$RNA@data) %in% Averages[i,1], colnames(eWAT@assays$RNA@data) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Annotation != "Adipocyte" & eWAT@meta.data$Replicate == "R2",])])))
  Averages[i,4] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "1" & Adipocytes@meta.data$Replicate == "R1",])])))
  Averages[i,5] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "1" & Adipocytes@meta.data$Replicate == "R2",])])))
  Averages[i,6] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "2" & Adipocytes@meta.data$Replicate == "R1",])])))
  Averages[i,7] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "2" & Adipocytes@meta.data$Replicate == "R2",])])))
  Averages[i,8] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "3" & Adipocytes@meta.data$Replicate == "R1",])])))
  Averages[i,9] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "3" & Adipocytes@meta.data$Replicate == "R2",])])))
}
rownames(Averages) <- Averages[,1]

# Scale each experiment
Exp1 <- Averages[,c(3,5,7,9)]
Exp1 <- t(scale(t(Exp1)))
Exp2 <- Averages[,c(2,4,6,8)]
Exp2 <- t(scale(t(Exp2)))
Averages[,c(3,5,7,9)] <- Exp1
Averages[,c(2,4,6,8)] <- Exp2

# Plot a heatmap
Heatmap(Averages[,c(2,3,6,7,4,5,8,9)], cluster_columns=F)

## Marker gene heatmap
# Select marker genes for each cell type
LSA <- data.frame(gene = c("Nnat", "Lrp3", "Car3", "Abcd2"), type = "LSA")
LGA <- data.frame(gene = c("Pparg", "Acaca", "Elovl6","Acly","Igf2r"), type = "LGA")
Shared <- data.frame(gene = c("Apoe", "Abcg1", "Trem2","Cd36"), type = "Shared")
SLSA <- data.frame(gene = c("Hif1a", "Nedd9", "Gadd45g","Rab7","Lep"), type = "SLSA")

# Combine results
Averages <- rbind(LSA, LGA, Shared, SLSA)

# Calculate average expression for selected marker genes
for (i in 1:nrow(Averages)) {
  Averages[i,3] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "1" & Adipocytes@meta.data$Replicate == "R1",])])))
  Averages[i,4] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "1" & Adipocytes@meta.data$Replicate == "R2",])])))
  Averages[i,5] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "2" & Adipocytes@meta.data$Replicate == "R1",])])))
  Averages[i,6] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "2" & Adipocytes@meta.data$Replicate == "R2",])])))
  Averages[i,7] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "3" & Adipocytes@meta.data$Replicate == "R1",])])))
  Averages[i,8] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[ rownames(Adipocytes@assays$RNA@data) %in% Averages[i,1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[ Adipocytes@meta.data$Label == "3" & Adipocytes@meta.data$Replicate == "R2",])])))
}
rownames(Averages) <- Averages[,1]

# Scale each experiment
Exp1 <- Averages[,c(3,5,7)]
Exp1 <- t(scale(t(Exp1)))
Exp2 <- Averages[,c(4,6,8)]
Exp2 <- t(scale(t(Exp2)))
Averages[,c(3,5,7)] <- Exp1
Averages[,c(4,6,8)] <- Exp2

# Plot the heatmap
Heatmap(Averages[c(5:9,1:4,10:18),c(5,6,3,4,7,8)], cluster_columns=F, cluster_rows=F)

## Gene module scoring
# Extract marker gene results
Clust1 <- Result[[1]] 
Clust2 <- Result[[2]] 
Clust3 <- Result[[3]] 

# Subset to genes that are enriched or exclusive (!= NS)
Clust1 <- Clust1[Clust1$Marker != "NS", ]
Clust2 <- Clust2[Clust2$Marker != "NS", ]
Clust3 <- Clust3[Clust3$Marker != "NS", ]

# Extracting the minimum logFC across replicates in pairwise cluster comparisons
# RATIONALE: The genes with the highest minimum fold change are the most specific ones for any given cluster
Clust1$logFC_OvO <- apply(Clust1[, grep("logFC_Cluster", colnames(Clust1))], 1, FUN = "min")
Clust2$logFC_OvO <- apply(Clust2[, grep("logFC_Cluster", colnames(Clust2))], 1, FUN = "min")
Clust3$logFC_OvO <- apply(Clust3[, grep("logFC_Cluster", colnames(Clust3))], 1, FUN = "min")

# Computing gene module scores using the top 50 most specific marker genes
Adipocytes <- AddModuleScore(Adipocytes, features = list(Cluster1 = Clust1[order(-Clust1$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster2 = Clust2[order(-Clust2$logFC_OvO), ][1:50, ]$Symbol,
                                             Cluster3 = Clust3[order(-Clust3$logFC_OvO), ][1:50, ]$Symbol))

# Scale the module scores (across all nuclei)
Adipocytes$Cluster1 <- scale(Adipocytes$Cluster1)
Adipocytes$Cluster2 <- scale(Adipocytes$Cluster2)
Adipocytes$Cluster3 <- scale(Adipocytes$Cluster3)

# Computing mean scaled module scores in each cluster in each replicate
Averages <- as.data.frame(matrix(ncol=length(Result), nrow=length(Result)*2))

# Set column and row names
Labels <- Adipocytes@meta.data[, c("Subtype","Label")]
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
      Averages[Counter,module] <- mean(Adipocytes@meta.data[ Adipocytes@meta.data$Label == label & Adipocytes$Replicate == rep,paste("Cluster",module,sep="")])
      Counter <- Counter + 1
    }
  }
}

## Plot mean scaled module scores
col <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
Heatmap(as.matrix(Averages), cluster_columns = F, cluster_rows = F, col = col)

## Subtype proportions
# Calculate proportions in each dataset
Composition <- data.frame(matrix(ncol=length(unique(Adipocytes$Subtype)), nrow=4))
colnames(Composition) <- unique(Adipocytes$Subtype)
for (i in 1:length(unique(Adipocytes$Subtype))) {
  Composition[1,i] <- nrow(Adipocytes@meta.data[ Adipocytes@meta.data$Dataset == "LFD_R1" & Adipocytes@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Adipocytes@meta.data[ Adipocytes@meta.data$Dataset == "LFD_R1",])
  Composition[2,i] <- nrow(Adipocytes@meta.data[ Adipocytes@meta.data$Dataset == "LFD_R2" & Adipocytes@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Adipocytes@meta.data[ Adipocytes@meta.data$Dataset == "LFD_R2",])
  Composition[3,i] <- nrow(Adipocytes@meta.data[ Adipocytes@meta.data$Dataset == "HFD_R1" & Adipocytes@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Adipocytes@meta.data[ Adipocytes@meta.data$Dataset == "HFD_R1",])
  Composition[4,i] <- nrow(Adipocytes@meta.data[ Adipocytes@meta.data$Dataset == "HFD_R2" & Adipocytes@meta.data$Subtype == colnames(Composition)[i],]) / nrow(Adipocytes@meta.data[ Adipocytes@meta.data$Dataset == "HFD_R2",])
}

# Calculate stats (mean, SEM, P)
Stats <- data.frame(matrix(ncol=5, nrow=length(unique(Adipocytes$Subtype))))
colnames(Stats) <- c("LFD_mean","LFD_SEM","HFD_mean","HFD_SEM","Pvalue")
rownames(Stats) <- colnames(Composition)
for (i in 1:length(unique(Adipocytes$Subtype))) {
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

### Deconvolution
## Mouse data
# Process microarray data
gset <- getGEO("GSE39549", GSEMatrix =TRUE, getGPL=TRUE, AnnotGPL=TRUE)
Data <- gset$GSE39549_series_matrix.txt.gz@assayData$exprs
Features <- pData(gset$GSE39549_series_matrix.txt.gz@featureData)
Samples <- pData(gset$GSE39549_series_matrix.txt.gz@phenoData)
Data <- as.data.frame(Data[ , colnames(Data) %in% rownames(Samples[ Samples[,43] == "epididymal adipose tissue",])])
Data$gene <- Features[,3]
Data <- Data[ Data$gene != "",]
Data$Mean <- abs(apply(Data[, colnames(Data) %in% rownames(Samples[ Samples[,43] == "epididymal adipose tissue" & Samples[,44] == "High fat diet",])],1,FUN="mean")-apply(Data[, colnames(Data) %in% rownames(Samples[ Samples[,43] == "epididymal adipose tissue" & Samples[,44] == "Normal diet",])],1,FUN="mean"))
Data <- Data[ order(Data$gene, -Data$Mean),]
Data <- Data[ duplicated(Data$gene) == F,]
Features <- as.data.frame(Features)
Features <- Features[ rownames(Features) %in% rownames(Data),]
rownames(Features) <- Features[,3]
Features <- Features[ match(Data$gene, rownames(Features)),]
Samples <- as.data.frame(Samples)
Samples <- Samples[ rownames(Samples) %in% colnames(Data),]
rownames(Data) <- Data$gene
Data <- Data[,1:30]
Microarray <- Data

### Create pseudo-bulk from each subtype in each dataset
# Drop cell types not intrinsic to the tissue (i.e. epididymal and spermatozoa)
eWAT <- subset(eWAT, subset = Annotation %in% c("Immune","Adipocyte","FAP","Endothelial","Mesothelial"))

# For each subtype in eWAT, sum the raw counts
Subtypes <- c("NPVM","PVM","P-LAM","RM","LAM","CEM","DCs","B-cells","T-cells","FAP2","FAP3","FAP1","FAP4","LSA","SLSA","LGA","MC","IMC","LEC","VEC","EPC")
Counts <- as.data.frame(matrix(ncol=4, nrow=nrow(eWAT)))
Counter <- 1
for (i in Subtypes) {
  Counts[,Counter] <- rowSums(eWAT@assays$RNA@counts[, colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Subtype == i & eWAT@meta.data$Dataset == "LFD_R1",])])
  colnames(Counts)[Counter] <- paste(i, "LFD_R1",sep="_")
  Counter <- Counter + 1
  Counts[,Counter] <- rowSums(eWAT@assays$RNA@counts[, colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Subtype == i & eWAT@meta.data$Dataset == "LFD_R2",])])
  colnames(Counts)[Counter] <- paste(i, "LFD_R2",sep="_")
  Counter <- Counter + 1
  Counts[,Counter] <- rowSums(eWAT@assays$RNA@counts[, colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Subtype == i & eWAT@meta.data$Dataset == "HFD_R1",])])
  colnames(Counts)[Counter] <- paste(i, "HFD_R1",sep="_")
  Counter <- Counter + 1
  Counts[,Counter] <- rowSums(eWAT@assays$RNA@counts[, colnames(eWAT@assays$RNA@counts) %in% rownames(eWAT@meta.data[ eWAT@meta.data$Subtype == i & eWAT@meta.data$Dataset == "HFD_R2",])])
  colnames(Counts)[Counter] <- paste(i, "HFD_R2",sep="_")
  Counter <- Counter + 1
}

# Normalize the data (with edgeR)
Normalized <- cpm(calcNormFactors(DGEList(Counts)), log = T, normalized.lib.size=T)

# Average across replicates in each diet
LFD <- as.data.frame(matrix(ncol=length(Subtypes), nrow=nrow(eWAT)))
HFD <- as.data.frame(matrix(ncol=length(Subtypes), nrow=nrow(eWAT)))
colnames(LFD) <- Subtypes
colnames(HFD) <- Subtypes
for (i in 1:length(Subtypes)) { LFD[,i] <- log2(rowMeans(2^Normalized[,(4 * (i-1) + 1):(4 * (i-1) + 2)])) }
for (i in 1:length(Subtypes)) { HFD[,i] <- log2(rowMeans(2^Normalized[,(4 * i - 1):(4 * i)])) }
rownames(LFD) <- rownames(eWAT)
rownames(HFD) <- rownames(eWAT)

# Subset to genes detected in both snRNA-seq and the microarray study
LFD <- LFD[ rownames(LFD) %in% rownames(Microarray),]
HFD <- HFD[ rownames(HFD) %in% rownames(Microarray),]
Microarray <- Microarray[ rownames(Microarray) %in% rownames(LFD),]
Microarray <- Microarray[ match(rownames(LFD), rownames(Microarray)),]

# Average across diets
Average <- as.data.frame(matrix(ncol=length(Subtypes), nrow = nrow(LFD)))
rownames(Average) <- rownames(LFD)
colnames(Average) <- colnames(LFD)
for (i in 1:length(Subtypes)) { Average[,i] <- log2(rowMeans(2^cbind(HFD[,i], LFD[,i]))) }

# Dispersion across diets
Var <- as.data.frame(matrix(ncol=length(Subtypes), nrow = nrow(LFD)))
rownames(Var) <- rownames(LFD)
colnames(Var) <- colnames(LFD)
for (i in 1:length(Subtypes)) { Var[,i] <- apply(cbind(LFD[,i], HFD[,i]),1,FUN="sd")/apply(cbind(LFD[,i], HFD[,i]),1,FUN="mean") }

# Subset genes based on having low index of dispersion across diets (i.e. consistently expressed genes)
Reference <- Average[ rownames(Average) %in% names(which(apply(Var,1,FUN="max") < 0.5)),]
Query <- Microarray[ rownames(Microarray) %in% names(which(apply(Var,1,FUN="max") < 0.5)),]

# Estimate regression coefficients for each subtype in the query set (microarray) in each sample based on the reference (snRNA-seq) using limits to avoid extreme solutions using cross-validated ridge regression
Inferred <- as.data.frame(matrix(ncol=30, nrow=length(unique(eWAT$Subtype))+1))
for (i in 1:30) {
  Sample <- cv.glmnet(y = Query[,i], x = as.matrix(Reference), alpha = 0, intercept = TRUE, standardize = TRUE, lower.limits = -0.15, upper.limits = 0.15)
  Inferred[,i] <- coef(Sample, "lambda.min")[,1]
  colnames(Inferred)[i] <- colnames(Microarray)[i]
  rownames(Inferred) <- rownames(coef(Sample))
}

# Scale inferred coefficients between 0 and 1 to estimate proportions, after getting rid of the intercept
Inferred <- t(Inferred)
Inferred <- Inferred[,2:(length(unique(eWAT$Subtype))+1)]
for (i in 1:nrow(Inferred)) { Inferred[i,] <- (Inferred[i,]-min(Inferred[i,]))/sum(Inferred[i,]-min(Inferred[i,])) }

# Calculate the relative distribution of SLSA, LSA and LGA among the proportion of adipocytes
Sum <- rowSums(Inferred[, which(colnames(Inferred) %in% c("SLSA","LGA","LSA"))])
for (i in 1:nrow(Inferred)) { 
  Inferred[i,colnames(Inferred) == "SLSA"] <- Inferred[i,colnames(Inferred) == "SLSA"] / Sum[i]
  Inferred[i,colnames(Inferred) == "LSA"] <- Inferred[i,colnames(Inferred) == "LSA"] / Sum[i]
  Inferred[i,colnames(Inferred) == "LGA"] <- Inferred[i,colnames(Inferred) == "LGA"] / Sum[i]
  }

## Collect statistics (mean and SEM)
# Normal diet
Stats_ND <- as.data.frame(matrix(ncol=5, nrow=5))
TypeCounter <- 1
for (m in c("SLSA","LSA","LGA")) {
  WeekCounter <- 1
  for (i in c("2weeks","4weeks","8weeks","20weeks","24weeks")) {
    Stats_ND[WeekCounter,TypeCounter] <- mean(Inferred[ rownames(Inferred) %in% rownames(Samples[Samples[,42] %in%  i & Samples[,44] == "Normal diet",]),which(colnames(Inferred) == m)])
    Stats_ND[WeekCounter,TypeCounter+1] <- sd(Inferred[ rownames(Inferred) %in% rownames(Samples[Samples[,42] %in%  i & Samples[,44] == "Normal diet",]),which(colnames(Inferred) == m)])/sqrt(3)
    WeekCounter <- WeekCounter + 1
  }
  TypeCounter <- TypeCounter + 2
}

# High fat diet
Stats_HFD <- as.data.frame(matrix(ncol=5, nrow=5))
TypeCounter <- 1
for (m in c("SLSA","LSA","LGA")) {
  WeekCounter <- 1
  for (i in c("2weeks","4weeks","8weeks","20weeks","24weeks")) {
    Stats_HFD[WeekCounter,TypeCounter] <- mean(Inferred[ rownames(Inferred) %in% rownames(Samples[Samples[,42] %in%  i & Samples[,44] == "High fat diet",]),which(colnames(Inferred) == m)])
    Stats_HFD[WeekCounter,TypeCounter+1] <- sd(Inferred[ rownames(Inferred) %in% rownames(Samples[Samples[,42] %in%  i & Samples[,44] == "High fat diet",]),which(colnames(Inferred) == m)])/sqrt(3)
    WeekCounter <- WeekCounter + 1
  }
  TypeCounter <- TypeCounter + 2
}

# Plot the deconvolution results
par(mfcol=c(1,3))
# SLSA
plot(Stats_ND[,1], ylim=c(0,0.4), col="black", las=1, ylab = "Inferred fraction", type="o", cex = 1, pch = 16, main = "SLSA")
arrows(seq(1,5,by=1), Stats_ND[,1] + Stats_ND[,2], seq(1,5,by=1), Stats_ND[,1] - Stats_ND[,2], code = 3, length = 0.05, angle = 90)
lines(Stats_HFD[,1], type = "o", col = "blue", pch = 16)
arrows(seq(1,5,by=1), Stats_HFD[,1] + Stats_HFD[,2], seq(1,5,by=1), Stats_HFD[,1] - Stats_HFD[,2], code = 3, length = 0.05, angle = 90, col = "blue")

# LSA
plot(Stats_ND[,3], ylim=c(0,0.8), col="black", las=1, ylab = "Inferred fraction", type="o", cex = 1, pch = 16, main = "LSA")
arrows(seq(1,5,by=1), Stats_ND[,3] + Stats_ND[,4], seq(1,5,by=1), Stats_ND[,3] - Stats_ND[,4], code = 3, length = 0.05, angle = 90)
lines(Stats_HFD[,3], type = "o", col = "blue", pch = 16)
arrows(seq(1,5,by=1), Stats_HFD[,3] + Stats_HFD[,4], seq(1,5,by=1), Stats_HFD[,3] - Stats_HFD[,4], code = 3, length = 0.05, angle = 90, col = "blue")

# LGA
plot(Stats_ND[,5], ylim=c(0,0.4), col="black", las=1, ylab = "Inferred fraction", type="o", cex = 1, pch = 16, main = "LGA")
arrows(seq(1,5,by=1), Stats_ND[,5] + Stats_ND[,6], seq(1,5,by=1), Stats_ND[,5] - Stats_ND[,6], code = 3, length = 0.05, angle = 90)
lines(Stats_HFD[,5], type = "o", col = "blue", pch = 16)
arrows(seq(1,5,by=1), Stats_HFD[,5] + Stats_HFD[,6], seq(1,5,by=1), Stats_HFD[,5] - Stats_HFD[,6], code = 3, length = 0.05, angle = 90, col = "blue")

## Human data
# Process RNA-seq data
RNA <- read.delim("~/Human_Adipose.counts", header=T)
rownames(RNA) <- RNA$Symbol
Human <- as.data.frame(cpm(calcNormFactors(DGEList(RNA[,9:ncol(RNA)])), log=T, normalized.lib.sizes=T))
Human$Gene.name <- rownames(Human)

# Convert to mouse gene symbols using homology information (from Ensembl)
# For duplicated genes, keep the one with the largest fold change between lean and obese donors
Convert <- read.delim("~/Human_Mouse.txt", header=T)
Human <- merge(Human, Convert, by="Gene.name")
Human <- Human[ Human$Mouse.gene.name %in% rownames(Tmp),]
Human$Mean <- abs(rowMeans(Human[,grep("Lean", colnames(Human))]) - rowMeans(Human[,grep("Obese", colnames(Human))]))
Human <- Human[ order(Human$Mouse.gene.name, -Human$Mean),]
Human <- Human[ duplicated(Human$Mouse.gene.name) == F,]
rownames(Human) <- Human$Mouse.gene.name

### Process pseudo-bulk from each subtype in each dataset (generated above under mouse data)
# Average across replicates in each diet
LFD <- as.data.frame(matrix(ncol=length(Subtypes), nrow=nrow(eWAT)))
HFD <- as.data.frame(matrix(ncol=length(Subtypes), nrow=nrow(eWAT)))
colnames(LFD) <- Subtypes
colnames(HFD) <- Subtypes
for (i in 1:length(Subtypes)) { LFD[,i] <- log2(rowMeans(2^Normalized[,(4 * (i-1) + 1):(4 * (i-1) + 2)])) }
for (i in 1:length(Subtypes)) { HFD[,i] <- log2(rowMeans(2^Normalized[,(4 * i - 1):(4 * i)])) }
rownames(LFD) <- rownames(eWAT)
rownames(HFD) <- rownames(eWAT)

# Subset to genes detected in both snRNA-seq and the microarray study
LFD <- LFD[ rownames(LFD) %in% rownames(Human),]
HFD <- HFD[ rownames(HFD) %in% rownames(Human),]
Human <- Human[ rownames(Human) %in% rownames(LFD),]
Human <- Human[ match(rownames(LFD), rownames(Human)),]

# Average across diets
Average <- as.data.frame(matrix(ncol=length(Subtypes), nrow = nrow(LFD)))
rownames(Average) <- rownames(LFD)
colnames(Average) <- colnames(LFD)
for (i in 1:length(Subtypes)) { Average[,i] <- log2(rowMeans(2^cbind(HFD[,i], LFD[,i]))) }

# Dispersion across diets
Var <- as.data.frame(matrix(ncol=length(Subtypes), nrow = nrow(LFD)))
rownames(Var) <- rownames(LFD)
colnames(Var) <- colnames(LFD)
for (i in 1:length(Subtypes)) { Var[,i] <- apply(cbind(LFD[,i], HFD[,i]),1,FUN="sd")/apply(cbind(LFD[,i], HFD[,i]),1,FUN="mean") }

# Subset genes based on having low index of dispersion across diets (i.e. consistently expressed genes)
Reference <- Average[ rownames(Average) %in% names(which(apply(Var,1,FUN="max") < 0.5)),]
Query <- Human[ rownames(Human) %in% names(which(apply(Var,1,FUN="max") < 0.5)),2:13]

# Estimate regression coefficients for each subtype in the query set (rna-seq) in each sample based on the reference (snRNA-seq) using limits to avoid extreme solutions using cross-validated ridge regression
# RATIONALE: Limits set slightly broader for human dataset, as the expected consistently with mouse is lower. Therefore, the analysis should allow more extreme values.
Inferred <- as.data.frame(matrix(ncol=ncol(Query), nrow=length(Subtypes)+1))
for (i in 1:ncol(Query)) {
  Sample <- cv.glmnet(y = 2^Query[,i], x = 2^as.matrix(Reference), alpha = 0, intercept = TRUE, standardize = TRUE, lower.limits = -0.25, upper.limits = 0.25)
  Inferred[,i] <- coef(Sample, "lambda.min")[,1]
  colnames(Inferred)[i] <- colnames(Normalized)[i]
  rownames(Inferred) <- rownames(coef(Sample))
}

# Scale inferred coefficients between 0 and 1 to estimate proportions, after getting rid of the intercept
Inferred <- t(Inferred)
Inferred <- Inferred[,2:(length(unique(eWAT$Subtype))+1)]
for (i in 1:nrow(Inferred)) { Inferred[i,] <- (Inferred[i,]-min(Inferred[i,]))/sum(Inferred[i,]-min(Inferred[i,])) }

# Calculate the relative distribution of SLSA, LSA and LGA among the proportion of adipocytes
Sum <- rowSums(Inferred[, which(colnames(Inferred) %in% c("SLSA","LGA","LSA"))])
for (i in 1:nrow(Inferred)) { 
  Inferred[i,which(colnames(Inferred) == "SLSA")] <- Inferred[i,which(colnames(Inferred) == "SLSA")] / Sum[i]
  Inferred[i,which(colnames(Inferred) == "LSA")] <- Inferred[i,which(colnames(Inferred) == "LSA")] / Sum[i]
  Inferred[i,which(colnames(Inferred) == "LGA")] <- Inferred[i,which(colnames(Inferred) == "LGA")] / Sum[i]
}

# Calculate statistic (mean and SEM)
Stats <- as.data.frame(matrix(ncol=2, nrow=6))
Stats[1,1] <- mean(Inferred[1:6,which(colnames(Inferred) == "LGA")])
Stats[2,1] <- mean(Inferred[7:12,which(colnames(Inferred) == "LGA")])
Stats[3,1] <- mean(Inferred[1:6,which(colnames(Inferred) == "LSA")])
Stats[4,1] <- mean(Inferred[7:12,which(colnames(Inferred) == "LSA")])
Stats[5,1] <- mean(Inferred[1:6,which(colnames(Inferred) == "SLSA")])
Stats[6,1] <- mean(Inferred[7:12,which(colnames(Inferred) == "SLSA")])
Stats[1,2] <- sd(Inferred[1:6,which(colnames(Inferred) == "LGA")])/sqrt(6)
Stats[2,2] <- sd(Inferred[7:12,which(colnames(Inferred) == "LGA")])/sqrt(6)
Stats[3,2] <- sd(Inferred[1:6,which(colnames(Inferred) == "LSA")])/sqrt(6)
Stats[4,2] <- sd(Inferred[7:12,which(colnames(Inferred) == "LSA")])/sqrt(6)
Stats[5,2] <- sd(Inferred[1:6,which(colnames(Inferred) == "SLSA")])/sqrt(6)
Stats[6,2] <- sd(Inferred[7:12,which(colnames(Inferred) == "SLSA")])/sqrt(6)

# Plot the deconvolution results
Bplt <- barplot(Stats[,1], las=1, ylim=c(0,0.8))
arrows(Bplt, Stats[,1] + Stats[,2], Bplt, Stats[,1] - Stats[,2], length=0.05, code=3, angle=90) 

### Diet-dependent regulation of gene expression in FAP subtypes
# Extract diet-regulated genes from all SLSA and LSA
LSA_Diet <- Result[[1]]
LSA_Diet <- LSA_Diet[LSA_Diet$Diet != "NS", ]
SLSA_Diet <- Result[[3]]
SLSA_Diet <- SLSA_Diet[SLSA_Diet$Diet != "NS", ]

# Calculate average expression for diet-regulated genes
Averages <- as.data.frame(matrix(ncol = 1, nrow = length(unique(sort(c(LSA_Diet[, 1], SLSA_Diet[, 1]))))))
Averages[, 1] <- unique(sort(c(LSA_Diet[, 1], SLSA_Diet[, 1])))
for (i in 1:nrow(Averages)) {
  Averages[i, 2] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Label == "1" & Adipocytes@meta.data$Diet == "LFD" & Adipocytes@meta.data$Replicate == "R1", ])])))
  Averages[i, 3] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Label == "1" & Adipocytes@meta.data$Diet == "LFD" & Adipocytes@meta.data$Replicate == "R2", ])])))
  Averages[i, 4] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Label == "1" & Adipocytes@meta.data$Diet == "HFD" & Adipocytes@meta.data$Replicate == "R1", ])])))
  Averages[i, 5] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Label == "1" & Adipocytes@meta.data$Diet == "HFD" & Adipocytes@meta.data$Replicate == "R2", ])])))
  Averages[i, 6] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Label == "3" & Adipocytes@meta.data$Diet == "LFD" & Adipocytes@meta.data$Replicate == "R1", ])])))
  Averages[i, 7] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Label == "3" & Adipocytes@meta.data$Diet == "LFD" & Adipocytes@meta.data$Replicate == "R2", ])])))
  Averages[i, 8] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Label == "3" & Adipocytes@meta.data$Diet == "HFD" & Adipocytes@meta.data$Replicate == "R1", ])])))
  Averages[i, 9] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Label == "3" & Adipocytes@meta.data$Diet == "HFD" & Adipocytes@meta.data$Replicate == "R2", ])])))
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
ToHighlight <- c("Pparg","Lpl","Lipe","Car3","Adipor2","Fgf10","Nr3c1","Plin1","Irs1","Insr","Cd36","Lipa","Mif","Ccl6","Ccl8","Ccl9","Col1a1","Col1a2","C1qa","C1qb","Ccl2","Cxcl12","Cd300lg","Col5a2","Col6a1","Col6a2","Bmp3","Bmp2k","Lbp","Deptor")
set.seed(100)
HTM <- Heatmap(Averages[, 2:9], cluster_columns = F, km = 3, show_row_names = F, right_annotation = rowAnnotation(foo = anno_mark(at = which(rownames(Averages) %in% ToHighlight), labels = rownames(Averages)[which(rownames(Averages) %in% ToHighlight)])), clustering_distance_rows = "manhattan", row_km_repeats = 100)

# Extract the cluster order and save it
set.seed(100)
Clusters <- row_order(HTM)
saveRDS(Clusters, "eWAT_Adipocyte_Diet_Clusters.Rds")

# Plot the heatmap
set.seed(100)
HTM

### Adipokines
# Sources of information for lists of adipokines:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3881510/
# https://journals.physiology.org/doi/full/10.1152/ajpendo.00297.2015
Adipokines <- c("Lep", "Adipoq", "Retn", "Tnf", "Il6", "Il1b", "Rbp4", "Ccl2", "Serpine1", "Nampt", "Rarres2", "Cfd", "Agt", "Apln")

# Extracting mean expression values across replicates and diets
Averages <- data.frame(gene = Adipokines)
for (i in 1:nrow(Averages)) {
  Averages[i, 2] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Replicate == "R1" & Adipocytes@meta.data$Diet == "LFD", ])])))
  Averages[i, 3] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Replicate == "R1" & Adipocytes@meta.data$Diet == "HFD", ])])))
  Averages[i, 4] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Replicate == "R2" & Adipocytes@meta.data$Diet == "LFD", ])])))
  Averages[i, 5] <- log1p(mean(expm1(Adipocytes@assays$RNA@data[rownames(Adipocytes@assays$RNA@data) %in% Averages[i, 1], colnames(Adipocytes@assays$RNA@data) %in% rownames(Adipocytes@meta.data[Adipocytes@meta.data$Replicate == "R2" & Adipocytes@meta.data$Diet == "HFD", ])])))
}
rownames(Averages) <- Averages[, 1]
names(Averages)[2:5] <- c("LFD_R1", "HFD_R1", "LFD_R2", "HFD_R2")

# Remove genes without signal
Averages <- Averages[!is.na(Averages[, 2]), ]

# Calculate logFC in each replicate between diets
Averages[, "logFC_R1"] <- Averages[, 3] - Averages[, 2]
Averages[, "logFC_R2"] <- Averages[, 5] - Averages[, 4]

# Add mean logFC and SEM
Averages[, "Mean_logFC"] <- apply(Averages[, c(6, 7)], 1, FUN = "mean")
Averages[, "sem"] <- apply(Averages[, c(6, 7)], 1, FUN = "sd") / sqrt(2)

# Add maximum expression values across diets
Averages[, "Max"] <- apply(cbind(rowMeans(Averages[, c(2, 4)]), rowMeans(Averages[, c(3, 5)])), 1, FUN = "max")

# Plot logFC against maximum expression
plot(Averages[, 10], Averages[, 8], xlab = "Maximum log expression", ylab = "logFC", pch = 16, las = 1, main = "Adipokines")
abline(h = 0)
```
[Back to start](../README.md)<br>