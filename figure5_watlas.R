### Setup
## Options
options(future.globals.maxSize= 20000*1024^2)

## Dependencies
source("dependencies.R")

## Libraries
library(Seurat)
library(Matrix)
library(JOINTLY)
library(BiocParallel)
library(HGC)
library(presto)
library(SeuratDisk)
library(writexl)
library(gamlss)
library(emmeans)
library(BisqueRNA)
library(recount3)
library(compositions)
library(ComplexHeatmap)
library(edgeR)

## Functions
source("functions.R")

##### PROCESS DATASETS
#### Tabula Sapiens data (https://pubmed.ncbi.nlm.nih.gov/35549404/)
#### Download from here: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5
### Prepare the dataset 
## Import the data
data <- readRDS("data/Adipose/raw/TS/local.rds")
data <- subset(data, tissue == "adipose tissue")
data$donor_id <- as.character(data$donor_id)

## Quality control
counts <- data@assays$RNA@counts
logCounts <- log2(colSums(counts))
logFeatures <- log2(colSums(counts > 0))
res <- resid(lm(logFeatures ~ logCounts))
plot(logCounts, logFeatures, col = ifelse(abs(res) >= 2, "red","blue"))
keep <- names(which(abs(res) <= 2))
counts <- counts[,colnames(counts) %in% keep]
logCounts <- log2(colSums(counts))
logFeatures <- log2(colSums(counts > 0))
res <- resid(lm(logFeatures ~ logCounts))
plot(logCounts, logFeatures, col = ifelse(abs(res) >= 0.8, "red","blue"))
keep <- names(which(abs(res) <= 0.8))
counts <- counts[,colnames(counts) %in% keep]
logFeatures <- colSums(counts > 0)
percent.mt <- colSums(counts[ rownames(counts) %in% human[ human$Chr == "MT","Ensembl"],]) / colSums(counts)
plot(y = percent.mt, x = logFeatures)
keep <- names(which(percent.mt <= 0.2))
counts <- counts[,colnames(counts) %in% keep]
ts <- subset(data, cells = colnames(counts))

### Labelling
## Initial embedding using JOINTLY
ts <- NormalizeData(ts, verbose = FALSE)
preprocessed <- preprocess(ts, batch.var = "donor_id")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(ts), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
ts[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")

## Initial clusters and labels
ts <- RunUMAP(ts, reduction = "jointly", dims = 1:15)
snn <- SNN.Construction(ts@reductions$jointly@cell.embeddings, k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 12)
ts$clusters <- cl
ts$labels <- NA
ts@meta.data[ ts@meta.data$clusters == 1, "labels"] <- "FAP"
ts@meta.data[ ts@meta.data$clusters == 2, "labels"] <- "FAP"
ts@meta.data[ ts@meta.data$clusters == 3, "labels"] <- "Endothelial"
ts@meta.data[ ts@meta.data$clusters == 4, "labels"] <- "NKT"
ts@meta.data[ ts@meta.data$clusters == 5, "labels"] <- "FAP"
ts@meta.data[ ts@meta.data$clusters == 6, "labels"] <- "MPS"
ts@meta.data[ ts@meta.data$clusters == 7, "labels"] <- "FAP"
ts@meta.data[ ts@meta.data$clusters == 8, "labels"] <- "Endothelial"
ts@meta.data[ ts@meta.data$clusters == 9, "labels"] <- "MPS"
ts@meta.data[ ts@meta.data$clusters == 10, "labels"] <- "Endothelial"
ts@meta.data[ ts@meta.data$clusters == 11, "labels"] <- "Mural"
DimPlot(ts, label = TRUE, group.by="labels") & NoLegend()

## Recluster MPS cells
ts.mps <- subset(ts, labels == "MPS")
ts.mps <- NormalizeData(ts.mps, verbose = FALSE)
preprocessed <- preprocess(ts.mps, batch.var = "donor_id")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(ts.mps), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
ts.mps[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")
snn <- SNN.Construction(ts.mps@reductions$jointly@cell.embeddings[,1:15])
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 10)
ts.mps$clusters <- cl
ts@meta.data[ rownames(ts@meta.data) %in% rownames(ts.mps@meta.data[ ts.mps@meta.data$clusters %in% c(8),]), "labels"] <- "Doublets"

# Remove doublets and cluster again
ts.mps <- subset(ts, labels == "MPS")
ts.mps <- NormalizeData(ts.mps, verbose = FALSE)
preprocessed <- preprocess(ts.mps, batch.var = "donor_id")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(ts.mps), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
ts.mps[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")
snn <- SNN.Construction(ts.mps@reductions$jointly@cell.embeddings[,1:15])
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 10)
ts.mps$clusters <- cl
ts@meta.data[ rownames(ts@meta.data) %in% rownames(ts.mps@meta.data[ ts.mps@meta.data$clusters %in% c(10),]), "labels"] <- "Doublets"

# Remove doublets and cluster again
ts.mps <- subset(ts, labels == "MPS")
ts.mps <- NormalizeData(ts.mps, verbose = FALSE)
preprocessed <- preprocess(ts.mps, batch.var = "donor_id")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(ts.mps), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
ts.mps[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")
snn <- SNN.Construction(ts.mps@reductions$jointly@cell.embeddings[,1:15])
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 8)
ts.mps$clusters <- cl
ts@meta.data[ rownames(ts@meta.data) %in% rownames(ts.mps@meta.data[ ts.mps@meta.data$clusters %in% c(1,6,7,2,5),]), "labels"] <- "MonoMac"
ts@meta.data[ rownames(ts@meta.data) %in% rownames(ts.mps@meta.data[ ts.mps@meta.data$clusters %in% c(8),]), "labels"] <- "Doublets"
ts@meta.data[ rownames(ts@meta.data) %in% rownames(ts.mps@meta.data[ ts.mps@meta.data$clusters %in% c(3),]), "labels"] <- "MonoMac"
ts@meta.data[ rownames(ts@meta.data) %in% rownames(ts.mps@meta.data[ ts.mps@meta.data$clusters %in% c(4),]), "labels"] <- "DC"
ts <- subset(ts, labels != "Doublets")

## Recluster FAP cells
ts.fap <- subset(ts, labels == "FAP")
ts.fap <- NormalizeData(ts.fap, verbose = FALSE)
preprocessed <- preprocess(ts.fap, batch.var = "donor_id")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(ts.fap), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
ts.fap[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")
snn <- SNN.Construction(ts.fap@reductions$jointly@cell.embeddings[,1:15])
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 5)
ts.fap$clusters <- cl
DimPlot(ts.fap, label = TRUE, group.by="clusters")
ts@meta.data[ rownames(ts@meta.data) %in% rownames(ts.fap@meta.data[ ts.fap@meta.data$clusters %in% c(5),]), "labels"] <- "Endothelial"

## Save the final object
saveRDS(ts, "data/Adipose/processed/tabula_sapiens.rds")

#### Jaitin data (https://pubmed.ncbi.nlm.nih.gov/31257031/)
#### Download from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128518
### Prepare the dataset
## Import genes
load("data/Adipose/genes.rda")

## Quality control sample GSM3876308
counts <- readMM("data/Adipose/raw/Jaitin/GSM3876308_DAH170_matrix.mtx.gz")
barcodes <- read.delim("data/Adipose/raw/Jaitin/GSM3876308_DAH170_barcodes.tsv.gz", header = FALSE)
features <- read.delim("data/Adipose/raw/Jaitin/GSM3876308_DAH170_features.tsv.gz", header = FALSE)
colnames(counts) <- barcodes[,1]
rownames(counts) <- features[,1]
metadata <- read.delim("data/Adipose/raw/Jaitin/10xHumanImmuneMatrix.csv", sep=",")
colnames(counts) <- paste("DAH170_22", colnames(counts), sep="_")
counts <- counts[ , colnames(counts) %in% metadata$CellID,]
logCounts <- log2(colSums(counts))
logFeatures <- log2(colSums(counts > 0))
res <- resid(lm(logFeatures ~ logCounts))
plot(logCounts, logFeatures, col = ifelse(abs(res) >= 0.8, "red","blue"))
keep <- names(which(abs(res) <= 0.8))
counts <- counts[,colnames(counts) %in% keep]
logFeatures <- colSums(counts > 0)
percent.mt <- colSums(counts[ rownames(counts) %in% human[ human$Chr == "MT","Ensembl"],]) / colSums(counts)
plot(y = percent.mt, x = logFeatures)
keep <- names(which(percent.mt <= 0.15))
counts <- counts[,colnames(counts) %in% keep]
md <- metadata[ metadata$CellID %in% colnames(counts),]
md <- md[ match(colnames(counts), md$CellID),]
rownames(md) <- md$CellID
md <- md[,"cell.type",drop=FALSE]
seuA <- CreateSeuratObject(counts, project = "Jaitin", meta.data = md, min.cells = 0, min.features = 0)
seuA$sample <- "DAH170_22"
seuA$age <- 52
seuA$bmi <- 46
seuA$gender <- "female"
seuA$tissue <- "VAT"

## Quality control sample GSM3876309
counts <- readMM("data/Adipose/raw/Jaitin/GSM3876309_DAH170_55_matrix.mtx.gz")
barcodes <- read.delim("data/Adipose/raw/Jaitin/GSM3876309_DAH170_55_barcodes.tsv.gz", header = FALSE)
features <- read.delim("data/Adipose/raw/Jaitin/GSM3876309_DAH170_55_features.tsv.gz", header = FALSE)
colnames(counts) <- barcodes[,1]
rownames(counts) <- features[,1]
metadata <- read.delim("data/Adipose/raw/Jaitin/10xHumanImmuneMatrix.csv", sep=",")
colnames(counts) <- paste("DAH170_55", colnames(counts), sep="_")
counts <- counts[ , colnames(counts) %in% metadata$CellID,]
logCounts <- log2(colSums(counts))
logFeatures <- log2(colSums(counts > 0))
res <- resid(lm(logFeatures ~ logCounts))
plot(logCounts, logFeatures, col = ifelse(abs(res) >= 0.7, "red","blue"))
keep <- names(which(abs(res) <= 0.7))
counts <- counts[,colnames(counts) %in% keep]
logFeatures <- colSums(counts > 0)
percent.mt <- colSums(counts[ rownames(counts) %in% human[ human$Chr == "MT","Ensembl"],]) / colSums(counts)
plot(y = percent.mt, x = logFeatures)
keep <- names(which(percent.mt <= 0.15))
counts <- counts[,colnames(counts) %in% keep]
md <- metadata[ metadata$CellID %in% colnames(counts),]
md <- md[ match(colnames(counts), md$CellID),]
rownames(md) <- md$CellID
md <- md[,"cell.type",drop=FALSE]
seuB <- CreateSeuratObject(counts, project = "Jaitin", meta.data = md, min.cells = 0, min.features = 0)
seuB$sample <- "DAH170_55"
seuB$age <- 52
seuB$bmi <- 46
seuB$gender <- "female"
seuB$tissue <- "VAT"

## Merge datasets
jaitin <- merge(seuA, seuB)

### Labelling
## Initial embedding using JOINTLY
jaitin <- NormalizeData(jaitin, verbose = FALSE)
preprocessed <- preprocess(jaitin, batch.var = "sample")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(jaitin), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
jaitin[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")

## Clusters and labels
jaitin <- RunUMAP(jaitin, reduction = "jointly", dims = 1:15)
snn <- SNN.Construction(jaitin@reductions$jointly@cell.embeddings[,1:15], k = 10, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 12)
jaitin$clusters <- cl
jaitin$labels <- NA
jaitin@meta.data[ jaitin@meta.data$clusters == 1, "labels"] <- "NKT"
jaitin@meta.data[ jaitin@meta.data$clusters == 2, "labels"] <- "NKT"
jaitin@meta.data[ jaitin@meta.data$clusters == 3, "labels"] <- "MonoMac"
jaitin@meta.data[ jaitin@meta.data$clusters == 4, "labels"] <- "NKT"
jaitin@meta.data[ jaitin@meta.data$clusters == 5, "labels"] <- "NKT"
jaitin@meta.data[ jaitin@meta.data$clusters == 6, "labels"] <- "Bcell"
jaitin@meta.data[ jaitin@meta.data$clusters == 7, "labels"] <- "MonoMac"
jaitin@meta.data[ jaitin@meta.data$clusters == 8, "labels"] <- "NKT"
jaitin@meta.data[ jaitin@meta.data$clusters == 9, "labels"] <- "NKT"
jaitin@meta.data[ jaitin@meta.data$clusters == 10, "labels"] <- "DC"
jaitin@meta.data[ jaitin@meta.data$clusters == 11, "labels"] <- "MonoMac"
jaitin@meta.data[ jaitin@meta.data$clusters == 12, "labels"] <- "MonoMac"

## Save final object
saveRDS(jaitin, "data/Adipose/processed/jaitin.rds")

#### Vijay data (https://pubmed.ncbi.nlm.nih.gov/32066997/)
#### Download from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129363
### Prepare the dataset
## Import the data
data <- readMM("data/Adipose/raw/Vijay/GSE129363_Discovery_Cohort_matrix.mtx.gz")
barcodes <- read.delim("data/Adipose/raw/Vijay/GSE129363_Discovery_Cohort_barcodes.tsv.gz", header = FALSE)
genes <- read.delim("data/Adipose/raw/Vijay/GSE129363_Discovery_Cohort_genes.tsv.gz", header = FALSE)
rownames(data) <- genes[,1]
colnames(data) <- barcodes[,1]

## Quality control
logCounts <- log2(colSums(data))
logFeatures <- log2(colSums(data > 0))
plot(y = logCounts, x = logFeatures, pch=16)
keep <- names(which(logFeatures >= 6))
keep <- keep[ keep %in% names(which(logCounts >= 6))]
data <- data[,colnames(data) %in% keep]
logCounts <- log2(colSums(data))
logFeatures <- log2(colSums(data > 0))
res <- resid(lm(logFeatures ~ logCounts))
plot(logCounts, logFeatures, col = ifelse(abs(res) >= 1, "red","blue"))
keep <- names(which(abs(res) <= 1.0))
data <- data[,colnames(data) %in% keep]
logFeatures <- colSums(data > 0)
percent.mt <- colSums(data[ rownames(data) %in% human[ human$Chr == "MT","Ensembl"],]) / colSums(data)
plot(y = percent.mt, x = logFeatures)
keep <- names(which(percent.mt <= 0.2))
data <- data[,colnames(data) %in% keep]

## Metadata
metadata <- read.delim("data/Adipose/raw/Vijay/GSE129363_Discovery_Cohort_CellAnnotation.txt.gz", header = TRUE)
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]
colnames(metadata) <- c("sample", "condition", "tissue")
data <- data[ ,colnames(data) %in% rownames(metadata)]
metadata <- metadata[ rownames(metadata) %in% colnames(data),]
metadata <- metadata[ match(colnames(data), rownames(metadata)),]
vijay <- CreateSeuratObject(data, project = "Vijay", meta.data = metadata, min.cells = 0, min.features = 0)

## Remove type 2 diabetic individuals
vijay <- subset(vijay, condition == "NonDiabetic")

### Labelling
## Initial embedding using JOINTLY
vijay <- NormalizeData(vijay, verbose = FALSE)
preprocessed <- preprocess(vijay, batch.var = "sample")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(vijay), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
vijay[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")

## Clusters and labels
vijay <- RunUMAP(vijay, reduction = "jointly", dims = 1:15)
snn <- SNN.Construction(vijay@reductions$jointly@cell.embeddings, k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 5)
vijay$clusters <- cl
vijay$labels <- NA
vijay@meta.data[ vijay@meta.data$clusters == 1, "labels"] <- "Immune"
vijay@meta.data[ vijay@meta.data$clusters == 2, "labels"] <- "Endothelial"
vijay@meta.data[ vijay@meta.data$clusters == 3, "labels"] <- "FAP"
vijay@meta.data[ vijay@meta.data$clusters == 4, "labels"] <- "NKT"
vijay@meta.data[ vijay@meta.data$clusters == 5, "labels"] <- "Mesothelial"

## Recluster immune cells
vijay.immune <- subset(vijay, labels == "Immune")
snn <- SNN.Construction(vijay.immune@reductions$jointly@cell.embeddings[,1:15], k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 20)
vijay.immune$clusters <- cl
vijay@meta.data[ rownames(vijay@meta.data) %in% rownames(vijay.immune@meta.data[ vijay.immune@meta.data$clusters %in% c(12),]), "labels"] <- "Bcell"
vijay@meta.data[ rownames(vijay@meta.data) %in% rownames(vijay.immune@meta.data[ vijay.immune@meta.data$clusters %in% c(2,7,8,11,17,18, 20),]), "labels"] <- "DC"
vijay@meta.data[ rownames(vijay@meta.data) %in% rownames(vijay.immune@meta.data[ vijay.immune@meta.data$clusters %in% c(1, 4, 5, 10, 13, 15,16,19),]), "labels"] <- "MonoMac"
vijay@meta.data[ rownames(vijay@meta.data) %in% rownames(vijay.immune@meta.data[ vijay.immune@meta.data$clusters %in% c(3, 6, 9, 14, 17),]), "labels"] <- "MonoMac"

## Save final object
saveRDS(vijay, "data/Adipose/processed/vijay.rds")

#### Hildreth data (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8102391/)
#### Download from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155960
### Prepare the dataset
## Setup files and metadata
files <- dir("data/Adipose/raw/Hildreth/")
md <- data.frame(id = c("L1","L2","L3","O1","O2","O3"), gender = "female", bmi = c(20.5,34.8,16.9,33.4,30,31.6), age = c(39,32,43,39,36,35), tissue = "SAT")
seu.list <- list()

## Import files in a loop and create a list of Seurat objects with proper metadata
for (file in files) {
  counts <- read.csv(paste("/work/NMF_project/adipose/raw/", file, sep=""))
  rownames(counts) <- counts[,1]
  counts <- counts[,-1]
  counts <- as.matrix(counts)
  counts <- as(counts, "dgCMatrix")
  seuA <- CreateSeuratObject(counts, project = "Hildreth", min.cells = 0, min.features = 0)
  seuA$sample <- md[ md$id %in% substr(file, regexpr("-", file)+1,regexpr("-", file)+2),1]
  seuA$gender <- md[ md$id %in% substr(file, regexpr("-", file)+1,regexpr("-", file)+2),2]
  seuA$bmi <- md[ md$id %in% substr(file, regexpr("-", file)+1,regexpr("-", file)+2),3]
  seuA$age <- md[ md$id %in% substr(file, regexpr("-", file)+1,regexpr("-", file)+2),4]
  seuA$tissue <- md[ md$id %in% substr(file, regexpr("-", file)+1,regexpr("-", file)+2),5]
  seuA$enrichment <- substr(file, regexpr("_",file)+1, regexpr("_",file)+5)
  seu.list[[length(seu.list)+1]] <- seuA
}

## Combine objects
hildreth <- merge(seu.list[[1]], seu.list[2:length(seu.list)])

## Quality control
counts <- hildreth@assays$RNA@counts
logCounts <- log2(colSums(counts))
logFeatures <- log2(colSums(counts > 0))
res <- resid(lm(logFeatures ~ logCounts))
plot(sort(abs(res)))
plot(logCounts, logFeatures, col = ifelse(abs(res) >= 0.8, "red","blue"))
keep <- names(which(abs(res) <= 0.8))
counts <- counts[,colnames(counts) %in% keep]
logFeatures <- colSums(counts > 0)
percent.mt <- colSums(counts[ rownames(counts) %in% human[ human$Chr == "MT","Ensembl"],]) / colSums(counts)
plot(y = percent.mt, x = logFeatures)
keep <- names(which(percent.mt <= 0.15))
counts <- counts[,colnames(counts) %in% keep]
hildreth <- subset(hildreth, cells = colnames(counts))
hildreth$percent.mt <- colSums(hildreth@assays$RNA@counts[ rownames(hildreth@assays$RNA@counts) %in% human[ human$Chr == "MT","Ensembl"],]) / colSums(hildreth@assays$RNA@counts)
hildreth <- subset(hildreth, percent.mt <= 0.20)
hildreth <- subset(hildreth, nCount_RNA >= 500)

### Labelling
## Initial embedding using JOINTLY
hildreth <- NormalizeData(hildreth, verbose = FALSE)
preprocessed <- preprocess(hildreth, batch.var = "sample")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(hildreth), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
hildreth[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")

## Clusters and labels
hildreth <- RunUMAP(hildreth, reduction = "jointly", dims = 1:15)
snn <- SNN.Construction(hildreth@reductions$jointly@cell.embeddings[,1:15], k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 12)
hildreth$clusters <- cl
hildreth$labels <- NA
hildreth@meta.data[ hildreth@meta.data$clusters == 1, "labels"] <- "FAP"
hildreth@meta.data[ hildreth@meta.data$clusters == 2, "labels"] <- "FAP"
hildreth@meta.data[ hildreth@meta.data$clusters == 3, "labels"] <- "MPS"
hildreth@meta.data[ hildreth@meta.data$clusters == 4, "labels"] <- "MPS"
hildreth@meta.data[ hildreth@meta.data$clusters == 5, "labels"] <- "Endothelial"
hildreth@meta.data[ hildreth@meta.data$clusters == 6, "labels"] <- "Bcell"
hildreth@meta.data[ hildreth@meta.data$clusters == 7, "labels"] <- "NKT"
hildreth@meta.data[ hildreth@meta.data$clusters == 8, "labels"] <- "Mural"
hildreth@meta.data[ hildreth@meta.data$clusters == 9, "labels"] <- "NKT"
hildreth@meta.data[ hildreth@meta.data$clusters == 10, "labels"] <- "NKT"
hildreth@meta.data[ hildreth@meta.data$clusters == 11, "labels"] <- "NKT"
hildreth@meta.data[ hildreth@meta.data$clusters == 12, "labels"] <- "NKT"

## Recluster the MPS
hildreth.mps <- subset(hildreth, labels == "MPS")
snn <- SNN.Construction(hildreth.mps@reductions$jointly@cell.embeddings[,1:15], k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 20)
hildreth.mps$clusters <- cl
hildreth@meta.data[ rownames(hildreth@meta.data) %in% rownames(hildreth.mps@meta.data[ hildreth.mps@meta.data$clusters %in% c(13,15,12,17,10,14,3,5),]), "labels"] <- "DC"
hildreth@meta.data[ rownames(hildreth@meta.data) %in% rownames(hildreth.mps@meta.data[ hildreth.mps@meta.data$clusters %in% c(2,16),]), "labels"] <- "MonoMac"
hildreth@meta.data[ rownames(hildreth@meta.data) %in% rownames(hildreth.mps@meta.data[ hildreth.mps@meta.data$clusters %in% c(8,20,19,16,18,9,4,1,7,11,6),]), "labels"] <- "MonoMac"

## Recluster the B cells
hildreth.b <- subset(hildreth, labels == "Bcell")
snn <- SNN.Construction(hildreth.b@reductions$jointly@cell.embeddings[,1:15], k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 3)
hildreth.b$clusters <- cl
hildreth@meta.data[ rownames(hildreth@meta.data) %in% rownames(hildreth.b@meta.data[ hildreth.b@meta.data$clusters %in% c(1),]), "labels"] <- "NKT"
hildreth@meta.data[ rownames(hildreth@meta.data) %in% rownames(hildreth.b@meta.data[ hildreth.b@meta.data$clusters %in% c(2,3),]), "labels"] <- "Bcell"

## Save the final object
saveRDS(hildreth, "data/Adipose/processed/hildreth.rds")

#### Emont data (https://pubmed.ncbi.nlm.nih.gov/35296864/)
#### Download data from here: https://singlecell.broadinstitute.org/single_cell/study/SCP1376/a-single-cell-atlas-of-human-and-mouse-white-adipose-tissue
#### Dataset already QC'd
### Import data
emont <- readRDS("data/Adipose/raw/Emont/emont.rds")

### Labelling the SVF
## Initial embedding using JOINTLY
emont.SVF <- subset(emont, tissue == "SVF")
emont.SVF <- NormalizeData(emont.SVF, verbose = FALSE)
preprocessed <- preprocess(emont.SVF, batch.var = "individual")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(emont.SVF), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
emont.SVF[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")

## Clusters and labels
emont.SVF <- RunUMAP(emont.SVF, reduction = "jointly", dims = 1:15)
snn <- SNN.Construction(emont.SVF@reductions$jointly@cell.embeddings[,1:15], k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 9)
emont.SVF$clusters <- cl
emont.SVF$labels <- NA
emont.SVF@meta.data[ emont.SVF@meta.data$clusters == 1, "labels"] <- "Immune"
emont.SVF@meta.data[ emont.SVF@meta.data$clusters == 2, "labels"] <- "FAP"
emont.SVF@meta.data[ emont.SVF@meta.data$clusters == 3, "labels"] <- "NKT"
emont.SVF@meta.data[ emont.SVF@meta.data$clusters == 4, "labels"] <- "FAP"
emont.SVF@meta.data[ emont.SVF@meta.data$clusters == 5, "labels"] <- "FAP"
emont.SVF@meta.data[ emont.SVF@meta.data$clusters == 6, "labels"] <- "Endothelial"
emont.SVF@meta.data[ emont.SVF@meta.data$clusters == 7, "labels"] <- "Proliferating"
emont.SVF@meta.data[ emont.SVF@meta.data$clusters == 8, "labels"] <- "FAP"
emont.SVF@meta.data[ emont.SVF@meta.data$clusters == 9, "labels"] <- "FAP"

## Reclustering immune cells
emont.immune <- subset(emont.SVF, labels == "Immune")
snn <- SNN.Construction(emont.immune@reductions$jointly@cell.embeddings[,1:15], k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 20)
emont.immune$clusters <- cl
emont.SVF@meta.data[ rownames(emont.SVF@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(14),]), "labels"] <- "Mast"
emont.SVF@meta.data[ rownames(emont.SVF@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(1,6,7,8,15,16, 20),]), "labels"] <- "DC"
emont.SVF@meta.data[ rownames(emont.SVF@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(4, 10,13,17,18,19),]), "labels"] <- "MonoMac"
emont.SVF@meta.data[ rownames(emont.SVF@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(2,3,5,9,11,12),]), "labels"] <- "MonoMac"

## Save the final object
saveRDS(emont.SVF, "data/Adipose/processed/emont_svf.rds")

### Labelling the SVF
## Initial embedding using JOINTLY
emont.AT <- subset(emont, tissue != "SVF")
emont.AT <- NormalizeData(emont.AT, verbose = FALSE)
preprocessed <- preprocess(emont.AT, batch.var = "individual")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(emont.AT), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
emont.AT[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")

## Clusters and labels
emont.AT <- RunUMAP(emont.AT, reduction = "jointly", dims = 1:15)
snn <- SNN.Construction(emont.AT@reductions$jointly@cell.embeddings[,1:15], k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 9)
emont.AT$clusters <- cl
emont.AT$labels <- NA
emont.AT@meta.data[ emont.AT@meta.data$clusters == 1, "labels"] <- "FAP"
emont.AT@meta.data[ emont.AT@meta.data$clusters == 2, "labels"] <- "Mesothelial"
emont.AT@meta.data[ emont.AT@meta.data$clusters == 3, "labels"] <- "Mural"
emont.AT@meta.data[ emont.AT@meta.data$clusters == 4, "labels"] <- "Adipocyte"
emont.AT@meta.data[ emont.AT@meta.data$clusters == 5, "labels"] <- "Immune"
emont.AT@meta.data[ emont.AT@meta.data$clusters == 6, "labels"] <- "Endothelial"
emont.AT@meta.data[ emont.AT@meta.data$clusters == 7, "labels"] <- "Immune"
emont.AT@meta.data[ emont.AT@meta.data$clusters == 8, "labels"] <- "Immune"
emont.AT@meta.data[ emont.AT@meta.data$clusters == 9, "labels"] <- "Endothelial"
DimPlot(emont.AT, label = TRUE, group.by="labels") & NoLegend()

## Reclustering mural cells
emont.mural <- subset(emont.AT, labels == "Mural")
snn <- SNN.Construction(emont.mural@reductions$jointly@cell.embeddings, k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 2)
emont.mural$clusters <- cl
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.mural@meta.data[ emont.mural@meta.data$clusters %in% c(1),]), "labels"] <- "Mural"
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.mural@meta.data[ emont.mural@meta.data$clusters %in% c(2),]), "labels"] <- "Endothelial"

## Reclustering immune cells
emont.immune <- subset(emont.AT, labels == "Immune")
snn <- SNN.Construction(emont.immune@reductions$jointly@cell.embeddings, k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 3)
emont.immune$clusters <- cl
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(1),]), "labels"] <- "Immune"
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(2),]), "labels"] <- "NKT"
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(3),]), "labels"] <- "Endothelial"

## Reclustering of immune cells
emont.immune <- subset(emont.AT, labels == "Immune")
emont.immune <- NormalizeData(emont.immune, verbose = FALSE)
preprocessed <- preprocess(emont.immune, batch.var = "individual")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(emont.immune), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
emont.immune[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")
emont.immune <- RunUMAP(emont.immune, reduction = "jointly", dims=1:15)
snn <- SNN.Construction(emont.immune@reductions$jointly@cell.embeddings, k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 20)
emont.immune$clusters <- cl
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(12),]), "labels"] <- "Mast"
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(16),]), "labels"] <- "Bcell"
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(1,2,4,5,6,7,8,10,11,13,14,15,17,18,19),]), "labels"] <- "MonoMac"
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(3,20),]), "labels"] <- "MonoMac"
emont.AT@meta.data[ rownames(emont.AT@meta.data) %in% rownames(emont.immune@meta.data[ emont.immune@meta.data$clusters %in% c(9),]), "labels"] <- "DC"

## Save the final object
saveRDS(emont.AT, "data/Adipose/processed/emont_at.rds")

#### Barboza (https://www.biorxiv.org/content/10.1101/2022.06.29.496888v1.full)
#### The dataset used was shared by Carey Lumeng directly. The data is now downloadable from SCP here: https://singlecell.broadinstitute.org/single_cell/study/SCP1903/human-adipose-tissue-single-nuclear-rna-sequencing
#### Dataset already QC'd
### Import data
load("data/Adipose/raw/Barboza/lumeng_at_022023.seurat")

### Labelling the SVF
## Initial embedding using JOINTLY
DefaultAssay(vat_sat_lo_int) <- "RNA"
vat_sat_lo_int[["integrated"]] <- NULL
preprocessed <- preprocess(vat_sat_lo_int, batch.var = "patient")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(vat_sat_lo_int), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
vat_sat_lo_int[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")

## Remove low quality cells using clusters
vat_sat_lo_int <- RunUMAP(vat_sat_lo_int, reduction = "jointly", dims=1:15)
snn <- SNN.Construction(vat_sat_lo_int@reductions$jointly@cell.embeddings, k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 30)
vat_sat_lo_int$clusters <- cl
vat_sat_lo_int <- subset(vat_sat_lo_int, cells = rownames(vat_sat_lo_int@meta.data[ !(vat_sat_lo_int@meta.data$clusters %in% c(1,4,19, 15)),]))

## Reintegrate with JOINTLY
preprocessed <- preprocess(vat_sat_lo_int, batch.var = "patient")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(vat_sat_lo_int), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
vat_sat_lo_int[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")

## Clusters and labels
vat_sat_lo_int <- RunUMAP(vat_sat_lo_int, reduction = "jointly", dims=1:15)
snn <- SNN.Construction(vat_sat_lo_int@reductions$jointly@cell.embeddings, k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 20)
vat_sat_lo_int$clusters <- cl
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(vat_sat_lo_int@meta.data[ vat_sat_lo_int@meta.data$clusters %in% c(4,9),]), "labels"] <- "Doublets"
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(vat_sat_lo_int@meta.data[ vat_sat_lo_int@meta.data$clusters %in% c(5),]), "labels"] <- "Mural"
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(vat_sat_lo_int@meta.data[ vat_sat_lo_int@meta.data$clusters %in% c(14,6),]), "labels"] <- "FAP"
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(vat_sat_lo_int@meta.data[ vat_sat_lo_int@meta.data$clusters %in% c(1,2,7),]), "labels"] <- "Mesothelial"
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(vat_sat_lo_int@meta.data[ vat_sat_lo_int@meta.data$clusters %in% c(10,13),]), "labels"] <- "Endothelial"
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(vat_sat_lo_int@meta.data[ vat_sat_lo_int@meta.data$clusters %in% c(12,11,3,16,15),]), "labels"] <- "Adipocyte"
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(vat_sat_lo_int@meta.data[ vat_sat_lo_int@meta.data$clusters %in% c(19,20,8,18,17),]), "labels"] <- "Immune"

## Recluster immune cells
baboza.subset <- subset(vat_sat_lo_int, labels == "Immune")
preprocessed <- preprocess(baboza.subset, batch.var = "patient")
cpca <- cpca(preprocessed, bpparam = MulticoreParam())
prepared <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(kernel.list = prepared$kernels, snn.list = prepared$snn, rare.list = prepared$rareity, cpca.result = cpca, share.objects = FALSE, bpparam = MulticoreParam())
embed <- solved$Hmat
embed <- t(do.call("cbind", embed))
embed <- embed[ match(colnames(baboza.subset), rownames(embed)),]
embed <- t(scale(t(scale(embed))))
colnames(embed) <- paste("jointly_", 1:15, sep="")
baboza.subset[["jointly"]] <- CreateDimReducObject(embed, assay = "RNA")
baboza.subset <- RunUMAP(baboza.subset, reduction = "jointly", dims=1:15)
snn <- SNN.Construction(baboza.subset@reductions$jointly@cell.embeddings, k = 20, threshold = 1/15)
tree <- HGC.dendrogram(G = snn)
cl <- cutree(tree, k = 50)
baboza.subset$clusters <- cl
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(baboza.subset@meta.data[ baboza.subset@meta.data$clusters %in% c(21,26),]), "labels"] <- "Doublets"
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(baboza.subset@meta.data[ baboza.subset@meta.data$clusters %in% c(22,1),]), "labels"] <- "DC"
vat_sat_lo_int@meta.data[ rownames(vat_sat_lo_int@meta.data) %in% rownames(baboza.subset@meta.data[ baboza.subset@meta.data$clusters %in% c(42,43,40,30,34,38,18,49,23,45,16,17,3,20),]), "labels"] <- "NKT"
vat_sat_lo_int@meta.data[ vat_sat_lo_int@meta.data$labels == "Immune","labels"] <- "MonoMac"
vat_sat_lo_int <- subset(vat_sat_lo_int, cells = rownames(vat_sat_lo_int@meta.data[ !(vat_sat_lo_int@meta.data$labels %in% "Doublets"),]))

## Save final object
saveRDS(vat_sat_lo_int, "data/Adipose/processed/barboza.rds")

##### Harmonize DATASETS
## Combine objects and write files
# Read datasets
emont_svf <- readRDS("data/Adipose/processed/emont_svf.rds")
emont_at <- readRDS("data/Adipose/processed/emont_at.rds")
hildreth <- readRDS("data/Adipose/processed/hildreth.rds")
vijay <- readRDS("data/Adipose/processed/vijay.rds")
jaitin <- readRDS("data/Adipose/processed/jaitin.rds")
ts <- readRDS("data/Adipose/processed/tabula_sapiens.rds")
barboza <- readRDS("data/Adipose/processed/barboza.rds")

# Streamline metadata
metadata <- emont_svf@meta.data
metadata$enrichment <- "None"
metadata$sample <- metadata$individual
metadata <- metadata[,c("orig.ident","sample","tissue","depot", "enrichment","labels")]
metadata$orig.ident <- "Emont"
emont_svf@meta.data <- metadata

metadata <- emont_at@meta.data
metadata$enrichment <- "None"
metadata$sample <- paste(metadata$individual, metadata$depot, sep="_")
metadata <- metadata[,c("orig.ident","sample","tissue","depot", "enrichment","labels")]
metadata$tissue <- "AT"
metadata$orig.ident <- "Emont"
emont_at@meta.data <- metadata

metadata <- hildreth@meta.data
metadata$depot <- metadata$tissue
metadata$tissue <- "SVF"
metadata <- metadata[,c("orig.ident","sample","tissue","depot", "enrichment","labels")]
hildreth@meta.data <- metadata

metadata <- vijay@meta.data
metadata$enrichment <- "None"
metadata$depot <- metadata$tissue
metadata$tissue <- "SVF"
metadata <- metadata[,c("orig.ident","sample","tissue","depot", "enrichment","labels")]
vijay@meta.data <- metadata

metadata <- jaitin@meta.data
metadata$enrichment <- "None"
metadata$depot <- metadata$tissue
metadata$tissue <- "SVF"
metadata <- metadata[,c("orig.ident","sample","tissue","depot", "enrichment","labels")]
metadata$orig.ident <- "Jaitin"
jaitin@meta.data <- metadata

metadata <- ts@meta.data
metadata$enrichment <- "None"
metadata$depot <- "SAT"
metadata$tissue <- "SVF"
metadata$sample <- metadata$donor_id
metadata$orig.ident <- "Tabula_Sapiens"
metadata <- metadata[,c("orig.ident","sample","tissue","depot", "enrichment","labels")]
ts@meta.data <- metadata

metadata <- barboza@meta.data
metadata$enrichment <- "None"
metadata$depot <- toupper(metadata$depot)
metadata$tissue <- "AT"
metadata$sample <- metadata$patient.depot
metadata$orig.ident <- "Baboza"
metadata <- metadata[,c("orig.ident","sample","tissue","depot", "enrichment","labels")]
barboza@meta.data <- metadata

## Streamline gene names
# Tabula sapiens
load("data/Adipose/genes.rda")
counts <- ts@assays$RNA@counts
genes <- human[ human$Ensembl %in% rownames(counts),]
genes <- genes[ duplicated(genes$Symbol)==FALSE,]
genes <- genes[ genes$Symbol != "",]
genes <- genes[ !is.na(genes$Symbol),]
counts <- counts[ order(rowSums(counts), decreasing = TRUE),]
counts <- counts[ rownames(counts) %in% genes$Ensembl,]
genes <- genes[ duplicated(genes$Ensembl)==FALSE,]
genes <- genes[ match(rownames(counts), genes$Ensembl),]
rownames(counts) <- genes$Symbol
ts[["RNA"]] <- CreateAssayObject(counts = counts)

# Hildreth
load("data/Adipose/genes.rda")
counts <- hildreth@assays$RNA@counts
genes <- human[ human$Ensembl %in% rownames(counts),]
genes <- genes[ duplicated(genes$Symbol)==FALSE,]
genes <- genes[ genes$Symbol != "",]
genes <- genes[ !is.na(genes$Symbol),]
counts <- counts[ order(rowSums(counts), decreasing = TRUE),]
counts <- counts[ rownames(counts) %in% genes$Ensembl,]
genes <- genes[ duplicated(genes$Ensembl)==FALSE,]
genes <- genes[ match(rownames(counts), genes$Ensembl),]
rownames(counts) <- genes$Symbol
hildreth[["RNA"]] <- CreateAssayObject(counts = counts)

# Vijay
load("data/Adipose/genes.rda")
counts <- vijay@assays$RNA@counts
genes <- human[ human$Ensembl %in% rownames(counts),]
genes <- genes[ duplicated(genes$Symbol)==FALSE,]
genes <- genes[ genes$Symbol != "",]
genes <- genes[ !is.na(genes$Symbol),]
counts <- counts[ order(rowSums(counts), decreasing = TRUE),]
counts <- counts[ rownames(counts) %in% genes$Ensembl,]
genes <- genes[ duplicated(genes$Ensembl)==FALSE,]
genes <- genes[ match(rownames(counts), genes$Ensembl),]
rownames(counts) <- genes$Symbol
vijay[["RNA"]] <- CreateAssayObject(counts = counts)

# Jaitin
load("data/Adipose/genes.rda")
counts <- jaitin@assays$RNA@counts
genes <- human[ human$Ensembl %in% rownames(counts),]
genes <- genes[ duplicated(genes$Symbol)==FALSE,]
genes <- genes[ genes$Symbol != "",]
genes <- genes[ !is.na(genes$Symbol),]
counts <- counts[ order(rowSums(counts), decreasing = TRUE),]
counts <- counts[ rownames(counts) %in% genes$Ensembl,]
genes <- genes[ duplicated(genes$Ensembl)==FALSE,]
genes <- genes[ match(rownames(counts), genes$Ensembl),]
rownames(counts) <- genes$Symbol
jaitin[["RNA"]] <- CreateAssayObject(counts = counts)

## Superset
all.genes <- c(rownames(jaitin), rownames(vijay), rownames(barboza), rownames(ts), rownames(hildreth), rownames(emont_svf))
all.genes <- names(which(table(all.genes) >= 3))

## Set default assays
DefaultAssay(emont_at) <- "RNA"
DefaultAssay(emont_svf) <- "RNA"
DefaultAssay(hildreth) <- "RNA"
DefaultAssay(vijay) <- "RNA"
DefaultAssay(ts) <- "RNA"
DefaultAssay(jaitin) <- "RNA"
DefaultAssay(barboza) <- "RNA"

## Get rid of additional assays
emont_at[["SCT"]] <- NULL
emont_at[["integrated"]] <- NULL
emont_svf[["SCT"]] <- NULL
emont_svf[["integrated"]] <- NULL

## Merge all datasets and fix label granularity
WATLAS <- merge(emont_at, list(emont_svf, hildreth, vijay, ts, jaitin, barboza))
WATLAS <- subset(WATLAS, features = all.genes)
gene.counts <- rowSums(WATLAS@assays$RNA@counts)
WATLAS <- subset(WATLAS, features = names(which(gene.counts >= 10)) )
cell.counts <- rowSums(WATLAS@assays$RNA@counts > 0)
WATLAS <- subset(WATLAS, features = names(which(cell.counts >= 5)) )
WATLAS@meta.data[ WATLAS@meta.data$labels == "Macrophage", "labels"] <- "MonoMac"
WATLAS@meta.data[ WATLAS@meta.data$labels == "Monocyte", "labels"] <- "MonoMac"

## Save objects as .rds and .h5ad
saveRDS(WATLAS, "data/Adipose/watlas.rds")
SaveH5Seurat(WATLAS, "data/Adipose/watlas.h5Seurat", overwrite = TRUE)
Convert("data/Adipose/watlas.h5ad", dest = "h5ad", overwrite = TRUE)

##### Supervised integration using scANVI (runInitialIntegration.py)

##### Relabel datasets at high granularity
## Import latent dimensions
latent <- read.csv("data/Adipose/latent/latent_initial.csv", header= FALSE)
rownames(latent) <- colnames(WATLAS)
colnames(latent) <- paste("ANVI_", seq(1,15, 1), sep="")
WATLAS[["anvi"]] <- CreateDimReducObject(as.matrix(latent), assay = "RNA")
WATLAS <- RunUMAP(WATLAS, reduction = "anvi", dims = 1:15)

## Cluster
WATLAS <- FindNeighbors(WATLAS, reduction = "anvi", dims = 1:15)
WATLAS <- FindClusters(WATLAS, res = 0.01)
DimPlot(WATLAS, label = TRUE, repel = TRUE) & NoLegend()

## Define fine and coarse labels
WATLAS$labels.l1 <- "NA"
WATLAS$labels.l2 <- "NA"
WATLAS@meta.data[ WATLAS@meta.data$seurat_clusters %in% c(0),"labels.l1"] <- "FAP"
WATLAS@meta.data[ WATLAS@meta.data$seurat_clusters %in% c(1),"labels.l1"] <- "NKT"
WATLAS@meta.data[ WATLAS@meta.data$seurat_clusters %in% c(2),"labels.l1"] <- "Mesothelial"
WATLAS@meta.data[ WATLAS@meta.data$seurat_clusters %in% c(3),"labels.l1"] <- "MPS"
WATLAS@meta.data[ WATLAS@meta.data$seurat_clusters %in% c(4),"labels.l1"] <- "Endothelial"
WATLAS@meta.data[ WATLAS@meta.data$seurat_clusters %in% c(5),"labels.l1"] <- "Adipocyte"

## Subclustering of adipocytes
subset <- subset(WATLAS, labels.l1 == "Adipocyte")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 0.1)
DimPlot(subset, label = TRUE)
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0),]),"labels.l1"] <- "Adipocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0),]),"labels.l2"] <- "Adipocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(1,2),]),"labels.l1"] <- "Low_quality"

## Re-subclustering of adipocytes
subset <- subset(WATLAS, labels.l1 == "Adipocyte")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 2)
DimPlot(subset, label = TRUE)

subset$labels.l2 <- "NOT"
subset@meta.data[ subset@meta.data$seurat_clusters %in% c(), "labels.l2"] <- "DGAT2"
subset@meta.data[ subset@meta.data$seurat_clusters %in% c(), "labels.l2"] <- "PPM1L"
subset@meta.data[ subset@meta.data$seurat_clusters %in% c(13), "labels.l2"] <- "DCN"
subset@meta.data[ subset@meta.data$seurat_clusters %in% c(), "labels.l2"] <- ""

WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(28),]),"labels.l1"] <- "Low_quality"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(1,20,19,15,10,9,4,8),]),"labels.l2"] <- "DGAT2+ Adipocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,5,6,11,12,17,18,22,23,24,26,16),]),"labels.l2"] <- "CLSTN2+ Adipocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(13),]),"labels.l2"] <- "DCN+ Adipocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(21,7,4,3,2,1,14,25,27),]),"labels.l2"] <- "PRSS23+ Adipocyte"

WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0),]),"labels.l2"] <- "Adipocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(1),]),"labels.l2"] <- "SRPX2+ Adipocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(2),]),"labels.l2"] <- "Adipocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(3),]),"labels.l2"] <- "DGAT2+ Adipocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(4),]),"labels.l1"] <- "Low_quality"

## Subclustering of mesothelial
subset <- subset(WATLAS, labels.l1 == "Mesothelial")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 0.1)
DimPlot(subset, label = TRUE)
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0),]),"labels.l1"] <- "Mesothelial"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0),]),"labels.l2"] <- "Mesothelial"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(1),]),"labels.l1"] <- "Low_quality"

## Subclustering of FAP
subset <- subset(WATLAS, labels.l1 == "FAP")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 2)
DimPlot(subset, label = TRUE)
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data),"labels.l2"] <- "FAP"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(29),]),"labels.l1"] <- "Endometrium"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(29),]),"labels.l2"] <- "Endometrium"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(26,32,35),]),"labels.l1"] <- "Low_quality"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(30),]),"labels.l1"] <- "Mesothelial"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(30),]),"labels.l2"] <- "Mesothelial"

## Subclustering of endothelial + mural
subset <- subset(WATLAS, labels.l1 == "Endothelial")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 1)
DimPlot(subset, label = TRUE)

WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(3),]),"labels.l1"] <- "LEC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(6,13,14),]),"labels.l1"] <- "SMC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(7,17),]),"labels.l1"] <- "Pericyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(11,12,15,16,18,21,22),]),"labels.l1"] <- "Low_quality"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,1,2,4,5,8,9,10,12,19),]),"labels.l1"] <- "VEC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(20),]),"labels.l1"] <- "Schwann"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(3),]),"labels.l2"] <- "LEC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(1,19,10,9),]),"labels.l2"] <- "ACKR1+ VEC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(4,2),]),"labels.l2"] <- "SOX5+ VEC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,8,5,12),]),"labels.l2"] <- "BTNL9+ VEC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(20),]),"labels.l2"] <- "Schwann"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(6,13,14),]),"labels.l2"] <- "SMC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(7,17),]),"labels.l2"] <- "Pericyte"

## Subclustering of Schwann cells
subset <- subset(WATLAS, labels.l1 == "Schwann")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 2)
DimPlot(subset, label = TRUE)
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(4),]),"labels.l1"] <- "Schwann"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(4),]),"labels.l2"] <- "Schwann"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,1,2,3,5),]),"labels.l1"] <- "Low_quality"

## Subclustering of NKT cells
subset <- subset(WATLAS, labels.l1 == "NKT")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 1)
DimPlot(subset, label = TRUE)
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(17),]),"labels.l1"] <- "Mast"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(17),]),"labels.l2"] <- "CPA3+ Mast"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(16),]),"labels.l1"] <- "Low_quality"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(20),]),"labels.l1"] <- "ILC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(20),]),"labels.l2"] <- "KIT+ ILC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(8,9,12,13),]),"labels.l1"] <- "NK cell"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(13),]),"labels.l2"] <- "XCL1+ NK cell"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(8,9,12),]),"labels.l2"] <- "CD247+ NK cell"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,2,3,4,5,14,15,18),]),"labels.l1"] <- "CD4+ T cell"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(18),]),"labels.l2"] <- "SELL+ CD4+ T cell" 
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(15),]),"labels.l2"] <- "CTLA4+ CD4+ T cell"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,2,4,5),]),"labels.l2"] <- "CD40LG+ CD4+ T cell"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(3,14),]),"labels.l2"] <- "SKAP1+ CD4+ T cell"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(1,6,7,10,11,19),]),"labels.l1"] <- "CD8+ T cell" 
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(11),]),"labels.l2"] <- "SLC4A10+ CD8+ T cell" 
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(7,10),]),"labels.l2"] <- "GNLY+ CD8+ T cell" 
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(1,6),]),"labels.l2"] <- "GZMK+ CD8+ T cell" 
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(19),]),"labels.l2"] <- "MT1E+ CD8+ T cell" 

## Subclustering of MPS cells
subset <- subset(WATLAS, labels.l1 == "MPS")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 1)
DimPlot(subset, label = TRUE)
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(11),]),"labels.l1"] <- "B cell"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(21),]),"labels.l1"] <- "pDC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(24),]),"labels.l1"] <- "Plasmablast"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(4,6,18,23),]),"labels.l1"] <- "DC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(3,12,14,15,16,17,19),]),"labels.l1"] <- "Monocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,1,2,5,7,8,9,10,13,22),]),"labels.l1"] <- "Macrophage"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(20),]),"labels.l2"] <- "Low_quality"

WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(11),]),"labels.l2"] <- "MS4A1+ B cell"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(21),]),"labels.l2"] <- "PLD4+ pDC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(24),]),"labels.l2"] <- "JCHAIN+ Plasmablast"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,1,7,8,9,22),]),"labels.l2"] <- "LYVE1+ Macrophage"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(2,13),]),"labels.l2"] <- "LPL+ Macrophage"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(5,10),]),"labels.l1"] <- "Macrophage"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(5,10),]),"labels.l2"] <- "CYP27A1+ Macrophage"

## Subclustering of DC
subset <- subset(WATLAS, labels.l1 == "DC")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 1)
DimPlot(subset, label = TRUE)
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(8,9,12),]),"labels.l2"] <- "CLEC9A+ DC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(11),]),"labels.l2"] <- "CCR7+ DC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(1,3,4,5,6,7,10),]),"labels.l2"] <- "CD1C+ DC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,2),]),"labels.l1"] <- "Monocyte"

## Subclustering of monocytes
subset <- subset(WATLAS, labels.l1 == "Monocyte")
subset <- FindNeighbors(subset, reduction = "anvi", dims = 1:15)
subset <- FindClusters(subset, res = 2)
DimPlot(subset, label = TRUE)
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data),"labels.l2"] <- "FCN1+ Monocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(17,4,24,21,5,22),]),"labels.l2"] <- "HES4+ Monocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(13,10),]),"labels.l2"] <- "CSF3R+ Monocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(13,10),]),"labels.l1"] <- "Monocyte"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,20),]),"labels.l2"] <- "CD1C+ DC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(0,20),]),"labels.l1"] <- "DC"
WATLAS@meta.data[ rownames(WATLAS@meta.data) %in% rownames(subset@meta.data[ subset@meta.data$seurat_clusters %in% c(3,15,16),]),"labels.l2"] <- "LYVE1+ Macrophage"

## Filter cells
WATLAS <- subset(WATLAS, labels.l1 != "Low_quality") # 8097 cells (2.6%)
WATLAS <- subset(WATLAS, labels.l2 != "Low_quality")

## Save final objects
saveRDS(WATLAS, "data/Adipose/watlas_filtered.Rds")
SaveH5Seurat(WATLAS, "data/Adipose/watlas_filtered.h5Seurat", overwrite = TRUE)
Convert("data/Adipose/watlas_filtered.h5Seurat", dest = "h5ad", overwrite = TRUE)

##### Supervised integration using scANVI (runFinalIntegration.py)

##### Final setup
### Import latent dimensions and mebed
latent <- read.csv("data/Adipose/latent/latent_final.csv", header= FALSE)
rownames(latent) <- colnames(WATLAS)
colnames(latent) <- paste("ANVI_", seq(1,15, 1), sep="")
WATLAS[["anvi"]] <- CreateDimReducObject(as.matrix(latent), assay = "RNA")
WATLAS <- RunUMAP(WATLAS, reduction = "anvi", dims = 1:15)
WATLAS <- NormalizeData(WATLAS)

### Regroup into major groups
WATLAS$group <- WATLAS$labels.l1
WATLAS@meta.data[ WATLAS@meta.data$labels.l1 %in% c("Plasmablast","ILC","pDC","NK cell","CD8+ T cell", "CD4+ T cell", "B cell"),"group"] <- "Lymphoid immune cells"
WATLAS@meta.data[ WATLAS@meta.data$labels.l1 %in% c("Monocyte","Macrophage","DC","Mast"),"group"] <- "Myeloid immune cells"
WATLAS@meta.data[ WATLAS@meta.data$labels.l1 %in% c("Pericyte","SMC"),"group"] <- "Mural"

### Fix metadata
# Setup
WATLAS$bmi <- "NA"
WATLAS$age <- "NA"
WATLAS$gender <- "NA"
WATLAS$wtstatus <- "NA"
WATLAS$dm <- "NA"

# Emont et al.
emont <- readRDS("data/Adipose/raw/Emont/emont.rds")
emont <- emont@meta.data
emont_Ad <- emont[ emont$tissue == "adipose",]
emont_Ad$sample <- paste(emont_Ad$individual, emont_Ad$depot, sep="_")
emont_Ad <- emont_Ad[ duplicated(emont_Ad$sample) == FALSE,]
for (donor in emont_Ad$sample) {
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"bmi"] <- emont_Ad[ emont_Ad$sample == donor,"bmi"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"age"] <- emont_Ad[ emont_Ad$sample == donor,"age"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"gender"] <- emont_Ad[ emont_Ad$sample == donor,"sex"]
}

emont_SVF <- emont[ emont$tissue == "SVF",]
emont_SVF <- emont_SVF[ duplicated(emont_SVF$individual) == FALSE,]
for (donor in emont_SVF$individual) {
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"bmi"] <- emont_SVF[ emont_SVF$individual == donor,"bmi"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"age"] <- emont_SVF[ emont_SVF$individual == donor,"age"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"gender"] <- emont_SVF[ emont_SVF$individual == donor,"sex"]
}
WATLAS@meta.data[ WATLAS@meta.data$age == "nan", "age"] <- "NA"

# Barboza et al.
barboza <- read.delim("data/Adipose/raw/Barboza/metadata.tsv")
barboza <- barboza[-1,]
barboza <- barboza[ duplicated(barboza$biosample_id) == FALSE,]
barboza$bmi <- "NA"
barboza$age <- "NA"
barboza$dm <- "NA"
barboza[ barboza$biosample_id == "1_vat","bmi"] <- 43
barboza[ barboza$biosample_id == "1_vat","age"] <- 52
barboza[ barboza$biosample_id == "1_vat","dm"] <- "DM"
barboza[ barboza$biosample_id == "1_sat","bmi"] <- 43
barboza[ barboza$biosample_id == "1_sat","age"] <- 52
barboza[ barboza$biosample_id == "1_sat","dm"] <- "DM"
barboza[ barboza$biosample_id == "2_vat","bmi"] <- 28
barboza[ barboza$biosample_id == "2_vat","age"] <- 46
barboza[ barboza$biosample_id == "2_vat","dm"] <- "NDM"
barboza[ barboza$biosample_id == "2_sat","bmi"] <- 28
barboza[ barboza$biosample_id == "2_sat","age"] <- 46
barboza[ barboza$biosample_id == "2_sat","dm"] <- "NDM"
barboza[ barboza$biosample_id == "3_vat","bmi"] <- 63
barboza[ barboza$biosample_id == "3_vat","age"] <- 34
barboza[ barboza$biosample_id == "3_vat","dm"] <- "NDM"
barboza[ barboza$biosample_id == "3_sat","bmi"] <- 63
barboza[ barboza$biosample_id == "3_sat","age"] <- 34
barboza[ barboza$biosample_id == "3_sat","dm"] <- "NDM"
barboza[ barboza$biosample_id == "4_vat","bmi"] <- 26
barboza[ barboza$biosample_id == "4_vat","age"] <- 64
barboza[ barboza$biosample_id == "4_vat","dm"] <- "NDM"
barboza[ barboza$biosample_id == "4_sat","bmi"] <- 26
barboza[ barboza$biosample_id == "4_sat","age"] <- 64
barboza[ barboza$biosample_id == "4_sat","dm"] <- "NDM"
barboza[ barboza$biosample_id == "5_vat","bmi"] <- 24
barboza[ barboza$biosample_id == "5_vat","age"] <- 47
barboza[ barboza$biosample_id == "5_vat","dm"] <- "NDM"
barboza[ barboza$biosample_id == "6_vat","bmi"] <- 43
barboza[ barboza$biosample_id == "6_vat","age"] <- 36
barboza[ barboza$biosample_id == "6_vat","dm"] <- "NDM"
barboza[ barboza$biosample_id == "6_sat","bmi"] <- 43
barboza[ barboza$biosample_id == "6_sat","age"] <- 36
barboza[ barboza$biosample_id == "6_sat","dm"] <- "NDM"
for (donor in barboza$biosample_id) {
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"wtstatus"] <- barboza[ barboza$biosample_id == donor,"wtstatus"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"gender"] <- barboza[ barboza$biosample_id == donor,"sex"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"age"] <- barboza[ barboza$biosample_id == donor,"age"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"bmi"] <- barboza[ barboza$biosample_id == donor,"bmi"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"dm"] <- barboza[ barboza$biosample_id == donor,"dm"]
}

# Hildreth et al.
hildreth <- readRDS("data/Adipose/processed/hildreth.rds")
hildreth <- hildreth@meta.data
hildreth <- hildreth[ duplicated(hildreth$sample) == FALSE,]
for (donor in hildreth$sample) {
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"bmi"] <- hildreth[ hildreth$sample == donor,"bmi"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"age"] <- hildreth[ hildreth$sample == donor,"age"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"gender"] <- hildreth[ hildreth$sample == donor,"gender"]
}

# Jaitin et al.
jaitin <- readRDS("data/Adipose/processed/jaitin.rds")
jaitin <- jaitin@meta.data
jaitin <- jaitin[ duplicated(jaitin$sample) == FALSE,]
for (donor in jaitin$sample) {
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"bmi"] <- jaitin[ jaitin$sample == donor,"bmi"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"age"] <- jaitin[ jaitin$sample == donor,"age"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"gender"] <- jaitin[ jaitin$sample == donor,"gender"]
}

# Tabula Sapiens
tabula_sapiens <- readRDS("data/Adipose/processed/tabula_sapiens.rds")
tabula_sapiens <- tabula_sapiens@meta.data
tabula_sapiens <- tabula_sapiens[ duplicated(tabula_sapiens$donor_id) == FALSE,]
tabula_sapiens$sex <- as.character(tabula_sapiens$sex)
for (donor in tabula_sapiens$donor_id) {
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"gender"] <- tabula_sapiens[ tabula_sapiens$donor_id == donor,"sex"]
}

# Vijay
vijay <- readRDS("data/Adipose/processed/vijay.rds")
vijay <- vijay@meta.data
vijay <- vijay[ duplicated(vijay$sample) == FALSE,]
vijay$condition <- as.character(factor(vijay$condition, labels = "NDM"))
for (donor in vijay$sample) {
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"dm"] <- vijay[ vijay$sample == donor,"condition"]
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"wtstatus"] <- "obese"
  WATLAS@meta.data[ WATLAS@meta.data$sample %in% donor,"gender"] <- "female"
}
WATLAS_vijay <- subset(WATLAS, orig.ident == "Vijay")
WATLAS_vijay <- SetIdent(WATLAS_vijay, value = WATLAS_vijay$sample)
WATLAS_vijay <- AverageExpression(WATLAS_vijay)$RNA
barplot(WATLAS_vijay[ rownames(WATLAS_vijay) == "DDX3Y",]) # Male specific
barplot(WATLAS_vijay[ rownames(WATLAS_vijay) == "XIST",]) # Female specific
WATLAS@meta.data[ WATLAS@meta.data$sample %in% c("Ate-12-VAT", "Ate-12-SAT", "Ate-14-VAT", "Ate-14-SAT"),"gender"] <- "male"
WATLAS@meta.data[ WATLAS@meta.data$bmi >= 35, "wtstatus"] <- "obese"
WATLAS@meta.data[ WATLAS@meta.data$bmi < 35, "wtstatus"] <- "lean"
WATLAS@meta.data[ WATLAS@meta.data$bmi == "NA", "wtstatus"] <- "NA"

## Save the final object
saveRDS(WATLAS, "data/Adipose/WATLAS_final_metadata.Rds")

#### Analysis of marker genes
### DotPlot to visualize
## Coarse labels
DotPlot(WATLAS, features = c("ACACB","LPL","ADIPOQ","DCN","PDGFRA","GSN", "MSLN","UPK3B","BNC1","PRLR","PGR","FAP", "XKR4","SOX10","NCAM1","PROX1","MMRN1","FLT4","VWF","CDH5","PECAM1","NOTCH3","KCNAB1","RGS6","CD14","C1QA","CD68","CD8A","IL7R","CD40LG"), group.by = "group")

## Fine labels
# Mural
subs <- subset(WATLAS, group == "Mural")
DotPlot(subs, features = c("MYOCD","MYH11","RCAN2","STEAP4","POSTN","COL25A1"), group.by = "labels.l2")

# VEC
subs <- subset(WATLAS, group == "VEC")
DotPlot(subs, features = c("ACKR1","TSHZ2","TLL1","BTNL9","CD36","ADGRF5","SOX5","GJA5","ARL15"), group.by = "labels.l2")

# FAP
subs <- subset(WATLAS, group == "FAP")
DotPlot(subs, features = c("DPP4","PI16","CD55","ICAM1","NAMPT","MT2A","PPARG","ACACA","RORA","CXCL14","APOD","GSN"), group.by = "labels.l2")

# Adipocytes
subs <- subset(WATLAS, group == "Adipocyte")
DotPlot(subs, features = c("PRSS23", "PDE5A", "ADIPOQ", "DGAT2","GPAM", "LPL", "DCN", "VIM", "RORA", "CLSTN2","CRIM1", "ITGA1"), group.by = "labels.l2")

# Myeloid cells
subs <- subset(WATLAS, group == "Myeloid immune cells")
DotPlot(subs, features = c("CPA3","HDC","KIT","LAMP3","CCR7", "EBI3","CD1C","CLEC10A","FCER1A","CLEC9A","CADM1","BATF3","HES4","FCGR3A","MBD2","FCN1","LYZ","VCAN","CSF3R","S100A8","CXCR2","LPL","TREM2","CD9","LYVE1","MRC1","F13A1","CYP27A1","DST","MARCO"), group.by = "labels.l2")

# NK cells
subs <- subset(WATLAS, labels.l2 %in% c("XCL1+ NK cell", "CD247+ NK cell"))
DotPlot(subs, features = c("CD247","FCGR3A","SPON2", "XCL1","CD160","CCL3"), group.by = "labels.l2")

# B cells
subs <- subset(WATLAS, labels.l2 %in% c("MS4A1+ B cell", "JCHAIN+ Plasmablast"))
DotPlot(subs, features = c("MS4A1","CD79A","CD19","JCHAIN","SDC1","MZB1"), group.by = "labels.l2")

# Other
subs <- subset(WATLAS, labels.l2 %in% c("KIT+ ILC", "PLD4+ pDC"))
DotPlot(subs, features = c("PLD4","SCT","LILRA4", "KIT","TLE1","LST1"), group.by = "labels.l2")

# CD4+ T cells
subs <- subset(WATLAS, labels.l2 %in% c("SKAP1+ CD4+ T cell", "SELL+ CD4+ T cell", "CTLA4+ CD4+ T cell","CD40LG+ CD4+ T cell"))
DotPlot(subs, features = c("FOXP3","CTLA4","TIGIT","SELL","CCR7","TSHZ2","CD40LG","GZMA","GZMK","SKAP1","PLCB1","ZBTB16"), group.by = "labels.l2")

# CD8+ T cells
subs <- subset(WATLAS, labels.l2 %in% c("SLC4A10+ CD8+ T cell","MT1E+ CD8+ T cell","GZMK+ CD8+ T cell","GNLY+ CD8+ T cell"))
DotPlot(subs, features = c("SLC4A10","CEBPD","NCR3","GZMK","SH2D1A","YBX3","GNLY","NKG7","KLRD1","MT1E","MT1X","MT2A"), group.by = "labels.l2")

## Tables
# Marker genes
markers <- wilcoxauc(WATLAS, group_by = "labels.l2")
groups <- unique(markers$group)
marker.list <- list()
for (group in groups) {
  tmp <- markers[ markers$group == group,]
  tmp <- tmp[,c(1,3,4,7,8,6,9,10)]
  marker.list[[length(marker.list)+1]] <- tmp
  names(marker.list)[length(marker.list)] <- group
}
writexl::write_xlsx(marker.list, col_names = TRUE, format_headers = TRUE, path = "/work/NMF_project/Figures/V1/Supplementary_Table_S1.xlsx")

# Pseudo-bulk
pseudoBulk <- AverageExpression(WATLAS, group.by = "labels.l2")$RNA
pseudoBulk <- as.data.frame(pseudoBulk)
writexl::write_xlsx(pseudoBulk, col_names = TRUE, format_headers = TRUE, path = "/work/NMF_project/Figures/V1/Supplementary_Table_S2.xlsx")

#### UMAPs
DimPlot(WATLAS, group.by = "labels.l2", label = TRUE, raster = FALSE) & NoLegend()
DimPlot(WATLAS, group.by = "orig.ident", raster = FALSE)
DimPlot(WATLAS, group.by = "depot", raster = FALSE)

#### Analysis of compositions
### Label granularity: L1
## Setup the metadata
md <- WATLAS@meta.data
md <- md[ duplicated(md$sample) == FALSE,]
md <- md[ md$wtstatus != "NA",]
mat_l1 <- as.matrix(table(WATLAS$sample, WATLAS$labels.l1))
mat_l1 <- mat_l1 / rowSums(mat_l1)
mat_l1 <- mat_l1[ rownames(mat_l1) %in% md$sample,]
comp <- mat_l1
md <- md[ match(rownames(comp), md$sample),]

## Calculate results using regression and marginal means
results <- data.frame(Celltype = "NA", Variable = "NA", Difference = 0, ci.min = 0, ci.max = 0)
counter <- 1
for (i in 1:20) {
  model <- lm(comp[,i] ~ md$wtstatus + md$depot + md$gender)
  out <- eff_size(emmeans(model, specs = "wtstatus"), sigma = sigma(model), edf = df.residual(model))
  results[counter,1] <- colnames(comp)[i]
  results[counter,2] <- "Weight"
  results[counter,3] <- summary(out)$effect.size
  results[counter,4] <- summary(out)$lower.CL
  results[counter,5] <- summary(out)$upper.CL
  counter <- counter + 1
  out <- eff_size(emmeans(model, specs = "gender"), sigma = sigma(model), edf = df.residual(model))
  results[counter,1] <- colnames(comp)[i]
  results[counter,2] <- "Gender"
  results[counter,3] <- summary(out)$effect.size
  results[counter,4] <- summary(out)$lower.CL
  results[counter,5] <- summary(out)$upper.CL
  counter <- counter + 1
  out <- eff_size(emmeans(model, specs = "depot"), sigma = sigma(model), edf = df.residual(model))
  results[counter,1] <- colnames(comp)[i]
  results[counter,2] <- "Depot"
  results[counter,3] <- summary(out)$effect.size
  results[counter,4] <- summary(out)$lower.CL
  results[counter,5] <- summary(out)$upper.CL
  counter <- counter + 1
}

# Plots
ggplot(data = results[ results$Variable == "Weight",], aes(x = Difference, y = Celltype, xmin = ci.min, xmax = ci.max)) + geom_point() + geom_pointrange() + theme_minimal()
ggplot(data = results[ results$Variable == "Depot",], aes(x = Difference, y = Celltype, xmin = ci.min, xmax = ci.max)) + geom_point() + geom_pointrange() + theme_minimal()
ggplot(data = results[ results$Variable == "Gender",], aes(x = Difference, y = Celltype, xmin = ci.min, xmax = ci.max)) + geom_point() + geom_pointrange() + theme_minimal()

# Tables
writexl::write_xlsx(results, "/work/NMF_project/Figures/V3/Supplementary_Table_Composition_L1.xlsx")

### Label granularity: L2
## Setup the metadata
mat_l2 <- as.matrix(table(WATLAS$sample, WATLAS$labels.l2))
mat_l2 <- mat_l2 / rowSums(mat_l2)
mat_l2 <- mat_l2[ rownames(mat_l2) %in% md$sample,]
comp <- mat_l2
md <- md[ match(rownames(comp), md$sample),]

## Calculate results using regression and marginal means
results_l2 <- data.frame(Celltype = "NA", State = "NA", Variable = "NA", Difference = 0, ci.min = 0, ci.max = 0)
counter <- 1
for (i in 1:43) {
    model <- lm(comp[,i] ~ md$wtstatus + md$depot + md$gender)
    out <- eff_size(emmeans(model, specs = "wtstatus"), sigma = sigma(model), edf = df.residual(model))
    results_l2[counter,1] <- names(which.max(table(WATLAS@meta.data[WATLAS@meta.data$labels.l2 %in% colnames(comp)[i],"labels.l1"])))
    results_l2[counter,2] <- colnames(comp)[i]
    results_l2[counter,3] <- "Weight"
    results_l2[counter,4] <- summary(out)$effect.size
    results_l2[counter,5] <- summary(out)$lower.CL
    results_l2[counter,6] <- summary(out)$upper.CL
    counter <- counter + 1
    out <- eff_size(emmeans(model, specs = "gender"), sigma = sigma(model), edf = df.residual(model))
    results_l2[counter,1] <- names(which.max(table(WATLAS@meta.data[WATLAS@meta.data$labels.l2 %in% colnames(comp)[i],"labels.l1"])))
    results_l2[counter,2] <- colnames(comp)[i]
    results_l2[counter,3] <- "Gender"
    results_l2[counter,4] <- summary(out)$effect.size
    results_l2[counter,5] <- summary(out)$lower.CL
    results_l2[counter,6] <- summary(out)$upper.CL
    counter <- counter + 1
    out <- eff_size(emmeans(model, specs = "depot"), sigma = sigma(model), edf = df.residual(model))
    results_l2[counter,1] <- names(which.max(table(WATLAS@meta.data[WATLAS@meta.data$labels.l2 %in% colnames(comp)[i],"labels.l1"])))
    results_l2[counter,2] <- colnames(comp)[i]
    results_l2[counter,3] <- "Depot"
    results_l2[counter,4] <- summary(out)$effect.size
    results_l2[counter,5] <- summary(out)$lower.CL
    results_l2[counter,6] <- summary(out)$upper.CL
    counter <- counter + 1
}

## Plots
ggplot(data = results_l2[ results_l2$Variable == "Weight",], aes(x = Difference, y = State, xmin = ci.min, xmax = ci.max)) + geom_point() + geom_pointrange() + theme_minimal()
ggplot(data = results_l2[ results_l2$Variable == "Depot",], aes(x = Difference, y = State, xmin = ci.min, xmax = ci.max)) + geom_point() + geom_pointrange() + theme_minimal()
ggplot(data = results_l2[ results_l2$Variable == "Gender",], aes(x = Difference, y = State, xmin = ci.min, xmax = ci.max)) + geom_point() + geom_pointrange() + theme_minimal()

## Tables
writexl::write_xlsx(results_l2, "/work/NMF_project/Figures/V3/Supplementary_Table_Composition_L2.xlsx")


#### Heatmaps
### Adipocyte marker genes
WATLAS$labels.new <- WATLAS$labels.l2
WATLAS@meta.data[ WATLAS@meta.data$labels.l1 != "Adipocyte","labels.new"] <- "Other"
WATLAS$id <- paste(WATLAS$labels.new, WATLAS$depot, sep="_")
ids <- unique(WATLAS$id)
pb <- as.data.frame(matrix(ncol=1, nrow=nrow(WATLAS)))
for (i in 1:length(ids)) {
  pb[,i] <- rowMeans(WATLAS@assays$RNA@data[,colnames(WATLAS@assays$RNA@data) %in% names(which(WATLAS$id == ids[i]))])
}
rownames(pb) <- rownames(WATLAS)
colnames(pb) <- ids
Heatmap(t(scale(t(pb[ rownames(pb) %in% c("LPL","ACSL1","PPARG","LIPE","PLIN4","DGAT2","DGAT1","PLIN1","LEP","CD36"),]))), cluster_columns = FALSE, cluster_rows = FALSE)   

### Adipocyte modules
## Subset to only adipocytes from Emont and Baboza et al.
ad <- subset(WATLAS, labels.l1 == "Adipocyte")
ad <- subset(ad, orig.ident %in% c("Emont","Baboza"))
ad$id <- paste(ad$labels.l2, ad$depot, ad$sample, sep="_")
ids <- unique(ad$id)

## Calculate pseudo-bulk and save metadata
pb <- as.data.frame(matrix(ncol=1, nrow=nrow(ad)))
md <- as.data.frame(matrix(ncol=5, nrow=1))
counter <- 1
for (i in 1:length(ids)) {
  if (table(ad$id)[names(table(ad$id)) == ids[i]] >= 5) {
    pb[,counter] <- rowSums(ad@assays$RNA@counts[,colnames(ad@assays$RNA@counts) %in% names(which(ad$id == ids[i]))])
    colnames(pb)[counter] <- ids[i]
    md[counter,1] <- ids[i]
    md[counter,2] <- unique(ad@meta.data[ ad@meta.data$id == ids[i],"orig.ident"])
    md[counter,3] <- unique(ad@meta.data[ ad@meta.data$id == ids[i],"sample"])
    md[counter,4] <- unique(ad@meta.data[ ad@meta.data$id == ids[i],"depot"])
    md[counter,5] <- unique(ad@meta.data[ ad@meta.data$id == ids[i],"labels.l2"])
    counter <- counter + 1
  }
}
rownames(pb) <- rownames(ad)

## Filter genes
keep <- intersect(names(which(rowSums(ad@assays$RNA@counts[,colnames(ad@assays$RNA@counts) %in% rownames(ad@meta.data[ ad@meta.data$orig.ident == "Emont",])] > 0) >= 20)), names(which(rowSums(ad@assays$RNA@counts[,colnames(ad@assays$RNA@counts) %in% rownames(ad@meta.data[ ad@meta.data$orig.ident == "Baboza",])] > 0) >= 20)))
pb <- pb[ rownames(pb) %in% keep,]
keep <- filterByExpr(pb[ ,colnames(pb) %in% md[ md$V2 == "Emont", "V1"]], group = md[ md$V2 == "Emont", "V5"])
pb.subset <- pb[ keep,]
keep <- filterByExpr(pb[ ,colnames(pb) %in% md[ md$V2 != "Emont", "V1"]], group = md[ md$V2 != "Emont", "V5"])
pb.subset <- pb.subset[ rownames(pb.subset) %in% rownames(pb)[keep],]
mm <- model.matrix(~ factor(md$V5) + factor(md$V4) + factor(md$V2) + 0)

## Estimate dispersions and fit a glm model
DGE <- DGEList(counts = pb.subset, group = md$V5)
DGE <- calcNormFactors(DGE)
DGE <- estimateDisp(DGE, design = mm)
fit <- glmFit(DGE,mm)

## CLSTN2+ adipocytes
# Pairwise tests
testA <- glmLRT(fit, contrast = c(1,-1,0,0,0,0))
testB <- glmLRT(fit, contrast = c(1,0,-1,0,0,0))
testC <- glmLRT(fit, contrast = c(1,0,0,-1,0,0))
testD <- glmLRT(fit, contrast = c(1,-1/3,-1/3,-1/3,0,0))
testA <- as.data.frame(topTags(testA, n = nrow(pb.subset), sort.by = "none"))
testB <- as.data.frame(topTags(testB, n = nrow(pb.subset), sort.by = "none"))
testC <- as.data.frame(topTags(testC, n = nrow(pb.subset), sort.by = "none"))
testD <- as.data.frame(topTags(testD, n = nrow(pb.subset), sort.by = "none"))
testD <- testD[ testD$FDR <= 0.05 & testD$logFC >= log2(1.5),]
testD <- testD[ rownames(testD) %in% rownames(testA)[ testA$logFC >= log2(1.3)],]
testD <- testD[ rownames(testD) %in% rownames(testB)[ testB$logFC >= log2(1.3)],]
testD <- testD[ rownames(testD) %in% rownames(testC)[ testC$logFC >= log2(1.3)],]
testD <- rownames(testD)

# Normalized expression values
norm <- as.data.frame(cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))
norm <- norm[ , order(md$V5, md$V2, md$V4)]
norm <- norm[ rownames(norm) %in% testD,]
norm <- removeBatchEffect(norm, batch = md$V4, batch2 = md$V2, design = model.matrix(~ factor(md$V5)))
norm <- t(scale(t(norm)))

# Setup and plot a heatmap
mod1 <- rownames(norm)
highlights <- c("CLSTN2","PGAP1", "ACER2", "LEPR", "ENPP1","CBLB","KSR1","CFAP69") #hAd5 from Emont et al.
htm1 <- Heatmap(t(as.matrix(norm)), cluster_columns = TRUE, cluster_rows = FALSE, left_annotation = rowAnnotation(Study = md[ order(md$V5, md$V2, md$V4), "V2"], Depot = md[ order(md$V5, md$V2, md$V4), "V4"]), top_annotation = HeatmapAnnotation(foo = anno_mark(at = which(rownames(norm) %in% highlights), labels = rownames(norm)[which(rownames(norm) %in% highlights)])), show_column_names = FALSE, show_column_dend = FALSE)

## DCN+ adipocytes
# Pairwise tests
testA <- glmLRT(fit, contrast = c(-1,1,0,0,0,0))
testB <- glmLRT(fit, contrast = c(0,1,-1,0,0,0))
testC <- glmLRT(fit, contrast = c(0,1,0,-1,0,0))
testD <- glmLRT(fit, contrast = c(-1/3,1,-1/3,-1/3,0,0))
testA <- as.data.frame(topTags(testA, n = nrow(pb.subset), sort.by = "none"))
testB <- as.data.frame(topTags(testB, n = nrow(pb.subset), sort.by = "none"))
testC <- as.data.frame(topTags(testC, n = nrow(pb.subset), sort.by = "none"))
testD <- as.data.frame(topTags(testD, n = nrow(pb.subset), sort.by = "none"))
testD <- testD[ testD$FDR <= 0.05 & testD$logFC >= log2(1.5),]
testD <- testD[ rownames(testD) %in% rownames(testA)[ testA$logFC >= log2(1.3)],]
testD <- testD[ rownames(testD) %in% rownames(testB)[ testB$logFC >= log2(1.3)],]
testD <- testD[ rownames(testD) %in% rownames(testC)[ testC$logFC >= log2(1.3)],]
testD <- rownames(testD)

# Normalized expression values
norm <- as.data.frame(cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))
norm <- norm[ , order(md$V5, md$V2, md$V4)]
norm <- norm[ rownames(norm) %in% testD,]
norm <- removeBatchEffect(norm, batch = md$V4, batch2 = md$V2, design = model.matrix(~ factor(md$V5)))
norm <- t(scale(t(norm)))

# Setup and plot a heatmap
mod2 <- rownames(norm)
highlights <- c("DCN","CXCL14","COL1A2","APOD","CFD","COL6A3","LUM")
htm2 <- Heatmap(t(as.matrix(norm)), cluster_columns = TRUE, cluster_rows = FALSE, left_annotation = rowAnnotation(Study = md[ order(md$V5, md$V2, md$V4), "V2"], Depot = md[ order(md$V5, md$V2, md$V4), "V4"]), top_annotation = HeatmapAnnotation(foo = anno_mark(at = which(rownames(norm) %in% highlights), labels = rownames(norm)[which(rownames(norm) %in% highlights)])), show_column_names = FALSE, show_column_dend = FALSE)

## DGAT2+ adipocytes
# Pairwise tests
testA <- glmLRT(fit, contrast = c(-1,0,1,0,0,0))
testB <- glmLRT(fit, contrast = c(0,-1,1,0,0,0))
testC <- glmLRT(fit, contrast = c(0,0,1,-1,0,0))
testD <- glmLRT(fit, contrast = c(-1/3,-1/3,1,-1/3,0,0))
testA <- as.data.frame(topTags(testA, n = nrow(pb.subset), sort.by = "none"))
testB <- as.data.frame(topTags(testB, n = nrow(pb.subset), sort.by = "none"))
testC <- as.data.frame(topTags(testC, n = nrow(pb.subset), sort.by = "none"))
testD <- as.data.frame(topTags(testD, n = nrow(pb.subset), sort.by = "none"))
testD <- testD[ testD$FDR <= 0.05 & testD$logFC >= log2(1.5),]
testD <- testD[ rownames(testD) %in% rownames(testA)[ testA$logFC >= log2(1.3)],]
testD <- testD[ rownames(testD) %in% rownames(testB)[ testB$logFC >= log2(1.3)],]
testD <- testD[ rownames(testD) %in% rownames(testC)[ testC$logFC >= log2(1.3)],]
testD <- rownames(testD)

# Normalized expression values
norm <- as.data.frame(cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))
norm <- norm[ , order(md$V5, md$V2, md$V4)]
norm <- norm[ rownames(norm) %in% testD,]
norm <- removeBatchEffect(norm, batch = md$V4, batch2 = md$V2, design = model.matrix(~ factor(md$V5)))
norm <- t(scale(t(norm)))

# Setup and plot a heatmap
mod3 <- rownames(norm)
highlights <- c("DGAT2","GPAM","MOGAT1","ENPP3","GYS2","FGFR2","FGF13","PRKCD","CDH20","CADM1") #hAd4 from Emont et al.
htm3 <- Heatmap(t(as.matrix(norm)), cluster_columns = TRUE, cluster_rows = FALSE, left_annotation = rowAnnotation(Study = md[ order(md$V5, md$V2, md$V4), "V2"], Depot = md[ order(md$V5, md$V2, md$V4), "V4"]), top_annotation = HeatmapAnnotation(foo = anno_mark(at = which(rownames(norm) %in% highlights), labels = rownames(norm)[which(rownames(norm) %in% highlights)])), show_column_names = FALSE, show_column_dend = FALSE)

## PRSS23+ adipocytes
# Pairwise tests
testA <- glmLRT(fit, contrast = c(-1,0,0,1,0,0))
testB <- glmLRT(fit, contrast = c(0,-1,0,1,0,0))
testC <- glmLRT(fit, contrast = c(0,0,-1,1,0,0))
testD <- glmLRT(fit, contrast = c(-1/3,-1/3,-1/3,1,0,0))
testA <- as.data.frame(topTags(testA, n = nrow(pb.subset), sort.by = "none"))
testB <- as.data.frame(topTags(testB, n = nrow(pb.subset), sort.by = "none"))
testC <- as.data.frame(topTags(testC, n = nrow(pb.subset), sort.by = "none"))
testD <- as.data.frame(topTags(testD, n = nrow(pb.subset), sort.by = "none"))
testD <- testD[ testD$FDR <= 0.05 & testD$logFC >= log2(1.5),]
testD <- testD[ rownames(testD) %in% rownames(testA)[ testA$logFC >= log2(1.3)],]
testD <- testD[ rownames(testD) %in% rownames(testB)[ testB$logFC >= log2(1.3)],]
testD <- testD[ rownames(testD) %in% rownames(testC)[ testC$logFC >= log2(1.3)],]
testD <- rownames(testD)

# Normalized expression values
norm <- as.data.frame(cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))
norm <- norm[ , order(md$V5, md$V2, md$V4)]
norm <- norm[ rownames(norm) %in% testD,]
norm <- removeBatchEffect(norm, batch = md$V4, batch2 = md$V2, design = model.matrix(~ factor(md$V5)))
norm <- t(scale(t(norm)))

# Setup and plot a heatmap
mod4 <- rownames(norm)
highlights <- c("PRSS23", "PDE5A", "ADIPOQ-AS1","SULT1C3","SULT2B1","SULT1C2","SULT1C4","GRIN2A","GRIN3A","GRIN2C") #hAd 3 from Emont
htm4 <- Heatmap(t(as.matrix(norm)), cluster_columns = TRUE, cluster_rows = FALSE, left_annotation = rowAnnotation(Study = md[ order(md$V5, md$V2, md$V4), "V2"], Depot = md[ order(md$V5, md$V2, md$V4), "V4"]), top_annotation = HeatmapAnnotation(foo = anno_mark(at = which(rownames(norm) %in% highlights), labels = rownames(norm)[which(rownames(norm) %in% highlights)])), show_column_names = FALSE, show_column_dend = FALSE)

## Plot the combined heatmap
htm1 + htm2 + htm3 + htm4

### Calculate module scores, setup for heatmap and plot it
ad <- UCell::AddModuleScore_UCell(ad, features = list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4))
md.full <- ad@meta.data
md.full$id <- paste(md.full$labels.l2, md.full$depot, md.full$orig.ident, sep="_")
md.full$mod1_UCell <- scale(md.full$mod1_UCell)
md.full$mod2_UCell <- scale(md.full$mod2_UCell)
md.full$mod3_UCell <- scale(md.full$mod3_UCell)
md.full$mod4_UCell <- scale(md.full$mod4_UCell)
modules <- aggregate(md.full[,c("mod1_UCell", "mod2_UCell", "mod3_UCell","mod4_UCell")], by = list(md.full$id), FUN="median")
rownames(modules) <- modules[,1]
modules <- modules[,-1]
Heatmap(as.matrix(modules), cluster_rows = FALSE, cluster_columns = FALSE)

### Adipocyte depots
# CLSTN2+ adipocytes
depot.counts <- pb.subset[,grep("CLSTN", colnames(pb.subset))]
depot.md <- md[match(colnames(depot.counts), md[,1]),]
mm <- model.matrix(~ factor(depot.md$V2) + factor(depot.md$V4))
DGE <- DGEList(counts = depot.counts, group = depot.md$V4)
DGE <- calcNormFactors(DGE)
DGE <- estimateDisp(DGE, design = mm)
fit <- glmFit(DGE,mm)
test <- glmLRT(fit, coef = "factor(depot.md$V4)VAT")
test <- as.data.frame(topTags(test, n = nrow(depot.counts), sort.by = "none"))
test <- test[ test$FDR <= 0.05 & abs(test$logFC) >= log2(1.5),]
norm_clstn <- as.data.frame(cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))
norm_clstn <- removeBatchEffect(norm_clstn, batch = depot.md$V2, design = model.matrix(~ factor(depot.md$V4)))
norm_clstn <- t(scale(t(norm_clstn)))
genes_clstn <- rownames(test)

# DGAT2+ adipocytes
depot.counts <- pb.subset[,grep("DGAT", colnames(pb.subset))]
depot.md <- md[match(colnames(depot.counts), md[,1]),]
mm <- model.matrix(~ factor(depot.md$V2) + factor(depot.md$V4))
DGE <- DGEList(counts = depot.counts, group = depot.md$V4)
DGE <- calcNormFactors(DGE)
DGE <- estimateDisp(DGE, design = mm)
fit <- glmFit(DGE,mm)
test <- glmLRT(fit, coef = "factor(depot.md$V4)VAT")
test <- as.data.frame(topTags(test, n = nrow(depot.counts), sort.by = "none"))
test <- test[ test$FDR <= 0.05 & abs(test$logFC) >= log2(1.5),]
norm_dgat <- as.data.frame(cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))
norm_dgat <- removeBatchEffect(norm_dgat, batch = depot.md$V2, design = model.matrix(~ factor(depot.md$V4)))
norm_dgat <- t(scale(t(norm_dgat)))
genes_dgat <- rownames(test)

# PRSS23+ adipocytes
depot.counts <- pb.subset[,grep("PRSS23", colnames(pb.subset))]
depot.md <- md[match(colnames(depot.counts), md[,1]),]
mm <- model.matrix(~ factor(depot.md$V2) + factor(depot.md$V4))
DGE <- DGEList(counts = depot.counts, group = depot.md$V4)
DGE <- calcNormFactors(DGE)
DGE <- estimateDisp(DGE, design = mm)
fit <- glmFit(DGE,mm)
test <- glmLRT(fit, coef = "factor(depot.md$V4)VAT")
test <- as.data.frame(topTags(test, n = nrow(depot.counts), sort.by = "none"))
test <- test[ test$FDR <= 0.05 & abs(test$logFC) >= log2(1.5),]
norm_prss <- as.data.frame(cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))
norm_prss <- removeBatchEffect(norm_prss, batch = depot.md$V2, design = model.matrix(~ factor(depot.md$V4)))
norm_prss <- t(scale(t(norm_prss)))
genes_prss <- rownames(test)

# DCN+ adipocytes
depot.counts <- pb.subset[,grep("DCN", colnames(pb.subset))]
depot.md <- md[match(colnames(depot.counts), md[,1]),]
mm <- model.matrix(~ factor(depot.md$V2) + factor(depot.md$V4))
DGE <- DGEList(counts = depot.counts, group = depot.md$V4)
DGE <- calcNormFactors(DGE)
DGE <- estimateDisp(DGE, design = mm)
fit <- glmFit(DGE,mm)
test <- glmLRT(fit, coef = "factor(depot.md$V4)VAT")
test <- as.data.frame(topTags(test, n = nrow(depot.counts), sort.by = "none"))
test <- test[ test$FDR <= 0.05 & abs(test$logFC) >= log2(1.5),]
norm_dcn <- as.data.frame(cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))
norm_dcn <- removeBatchEffect(norm_dcn, batch = depot.md$V2, design = model.matrix(~ factor(depot.md$V4)))
norm_dcn <- t(scale(t(norm_dcn)))
genes_dcn <- rownames(test)

## Combine genes and plot heatmaps
depot_genes <- unique(c(genes_dcn, genes_clstn, genes_dgat, genes_prss))

## Combine results and cluster
norm_counts <- cbind(norm_dcn[ rownames(norm_dcn) %in% depot_genes,],norm_dgat[ rownames(norm_dgat) %in% depot_genes,],norm_clstn[ rownames(norm_clstn) %in% depot_genes,],norm_prss[ rownames(norm_prss) %in% depot_genes,])
md_match <- md[ match(colnames(norm_counts), md[,1]),]
norm_counts <- norm_counts[ , order(md_match$V5, md_match$V4, md_match$V2)]
md_match <- md_match[order(md_match$V5, md_match$V4, md_match$V2),]
annot <- HeatmapAnnotation(Celltype = md_match$V5, Depot = md_match$V4, Study = md_match$V2)
Heatmap(norm_counts, cluster_columns = FALSE, km = 2, bottom_annotation = annot)

#### Deconvolution using Bisque
### SAT
## Setup the data
WATLAS.sat <- subset(WATLAS, depot == "SAT")
WATLAS.sat <- subset(WATLAS, labels.l1 != "Mesothelial")

## Setup GTEx data
human_projects <- available_projects()
adipose <- create_rse(human_projects[ grep("ADIPOSE", human_projects$project),])
assay(adipose, "counts") <- compute_read_counts(adipose)
counts <- assay(adipose, "counts")
genes <- as.data.frame(rowRanges(adipose))
genes <- genes[ duplicated(genes$gene_name)==FALSE,]
genes <- genes[ genes$gene_name != "",]
genes <- genes[ !is.na(genes$gene_name),]
counts <- counts[ order(rowSums(counts), decreasing = TRUE),]
counts <- counts[ rownames(counts) %in% genes$gene_id,]
genes <- genes[ duplicated(genes$gene_id)==FALSE,]
genes <- genes[ genes$gene_id %in% rownames(counts),]
genes <- genes[ match(rownames(counts), genes$gene_id),]
rownames(counts) <- genes$gene_name
counts <- counts[ ,grep("Subcut", adipose$gtex.smtsd)]

## Define states and celltypes
celltype.list = list(
  ILC = c("KIT+ ILC"),
  CD4_T = c("SKAP1+ CD4+ T cell","CD40LG+ CD4+ T cell","CTLA4+ CD4+ T cell","SELL+ CD4+ T cell"),
  CD8_T = c("GZMK+ CD8+ T cell","GNLY+ CD8+ T cell", "MT1E+ CD8+ T cell","SLC4A10+ CD8+ T cell"),
  NK = c("CD247+ NK cell","XCL1+ NK cell"),
  Mural = c("SMC","Pericyte"), 
  Adipose = c("PRSS23+ Adipocyte","CLSTN2+ Adipocyte","DCN+ Adipocyte","DGAT2+ Adipocyte"),
  FAP = c("CXCL14+ FAP","PPARG+ FAP","DPP4+ FAP","ICAM1+ FAP","Endometrium"),
  Vascular = c("ACKR1+ VEC","BTNL9+ VEC","SOX5+ VEC"),
  LEC = c("LEC"),
  Monocyte = c("FCN1+ Monocyte","CSF3R+ Monocyte","HES4+ Monocyte"),
  Mast = c("CPA3+ Mast"),
  Macrophage = c("CYP27A1+ Macrophage","LYVE1+ Macrophage","LPL+ Macrophage"),
  cDC = c("CLEC9A+ DC","CD1C+ DC","CCR7+ DC"),
  pDC = c("PLD4+ pDC"),
  Bcells = c("MS4A1+ B cell","JCHAIN+ Plasmablast"),
  Schwann = c("Schwann")
)

## Set deconGroup
WATLAS.sat$deconGroup <- "NA"
for (i in 1:length(celltype.list)) {
  WATLAS.sat@meta.data[ WATLAS.sat@meta.data$labels.l2 %in% celltype.list[[i]], "deconGroup"] <- names(celltype.list)[i]
}

## Define marker genes
markers <- wilcoxauc(WATLAS.sat, group_by = "deconGroup")
markers <- names(which(table(markers[ markers$auc >= 0.6,"feature"]) == 1))

## Setup scRNA-seq data
sc.pheno <- base::data.frame(check.names = F, check.rows = F,stringsAsFactors = F, row.names = rownames(WATLAS.sat@meta.data), SubjectName = WATLAS.sat$sample, cellType = WATLAS.sat$deconGroup)
sc.meta <- base::data.frame(labelDescription = base::c("SubjectName", "cellType"), row.names = base::c("SubjectName", "cellType"))
sc.pdata <- methods::new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta)
sc.data <- GetAssayData(WATLAS.sat, "count")
sc.data <- sc.data[ rownames(sc.data) %in% rownames(counts),]
sc.eset <- Biobase::ExpressionSet(assayData = as.matrix(sc.data), phenoData = sc.pdata)

## Setup bulk RNA-seq data
bulk.eset <- ExpressionSet(assayData = counts)

## Deconvolute
decon.sat <- ReferenceBasedDecomposition(bulk.eset = bulk.eset, sc.eset = sc.eset, markers = markers, use.overlap = FALSE)
saveRDS(decon.sat, "data/Adipose/decon/decon_SAT.Rds")

## Plots and stats
plot(rowMeans(decon.sat$bulk.props),rowMeans(decon.sat$sc.props), las = 1, xlab="Average deconvoluted fraction", ylab="Average scRNA-seq fraction", pch = 16, main = "SAT")
cor(rowMeans(decon.sat$bulk.props),rowMeans(decon.sat$sc.props)) # 0.998428
sqrt(mean((rowMeans(decon.sat$bulk.props)-rowMeans(decon.sat$sc.props))^2)) # 0.005592918

### VAT
## Setup the data
WATLAS.vat <- subset(WATLAS, depot == "VAT")

## Setup GTEx data
counts <- assay(adipose, "counts")
genes <- as.data.frame(rowRanges(adipose))
genes <- genes[ duplicated(genes$gene_name)==FALSE,]
genes <- genes[ genes$gene_name != "",]
genes <- genes[ !is.na(genes$gene_name),]
counts <- counts[ order(rowSums(counts), decreasing = TRUE),]
counts <- counts[ rownames(counts) %in% genes$gene_id,]
genes <- genes[ duplicated(genes$gene_id)==FALSE,]
genes <- genes[ genes$gene_id %in% rownames(counts),]
genes <- genes[ match(rownames(counts), genes$gene_id),]
rownames(counts) <- genes$gene_name
counts <- counts[ ,grep("Subcut", adipose$gtex.smtsd, invert = TRUE)]

## Define states and celltypes
celltype.list = list(
  ILC = c("KIT+ ILC"),
  CD4_T = c("SKAP1+ CD4+ T cell","CD40LG+ CD4+ T cell","CTLA4+ CD4+ T cell","SELL+ CD4+ T cell"),
  CD8_T = c("GZMK+ CD8+ T cell","GNLY+ CD8+ T cell", "MT1E+ CD8+ T cell","SLC4A10+ CD8+ T cell"),
  NK = c("CD247+ NK cell","XCL1+ NK cell"),
  Mural = c("SMC","Pericyte"), 
  Adipose = c("PRSS23+ Adipocyte","CLSTN2+ Adipocyte","DCN+ Adipocyte","DGAT2+ Adipocyte"),
  FAP = c("CXCL14+ FAP","PPARG+ FAP","DPP4+ FAP","ICAM1+ FAP","Endometrium"),
  Vascular = c("ACKR1+ VEC","BTNL9+ VEC","SOX5+ VEC"),
  LEC = c("LEC"),
  Monocyte = c("FCN1+ Monocyte","CSF3R+ Monocyte","HES4+ Monocyte"),
  Mast = c("CPA3+ Mast"),
  Macrophage = c("CYP27A1+ Macrophage","LYVE1+ Macrophage","LPL+ Macrophage"),
  cDC = c("CLEC9A+ DC","CD1C+ DC","CCR7+ DC"),
  pDC = c("PLD4+ pDC"),
  Bcells = c("MS4A1+ B cell","JCHAIN+ Plasmablast"),
  Schwann = c("Schwann"),
  Mesothelial = c("RACK1+ Mesothelial","PAMR1+ Mesothelial","LRRTM3+ Mesothelial")
)

## Set deconGroup
WATLAS.vat$deconGroup <- "NA"
for (i in 1:length(celltype.list)) {
  WATLAS.vat@meta.data[ WATLAS.vat@meta.data$labels.l2 %in% celltype.list[[i]], "deconGroup"] <- names(celltype.list)[i]
}

## Define marker genes
markers <- wilcoxauc(WATLAS.vat, group_by = "deconGroup")
markers <- names(which(table(markers[ markers$auc >= 0.6,"feature"]) == 1))

## Setup scRNA-seq data
sc.pheno <- base::data.frame(check.names = F, check.rows = F,stringsAsFactors = F, row.names = rownames(WATLAS.vat@meta.data), SubjectName = WATLAS.vat$sample, cellType = WATLAS.vat$deconGroup)
sc.meta <- base::data.frame(labelDescription = base::c("SubjectName", "cellType"), row.names = base::c("SubjectName", "cellType"))
sc.pdata <- methods::new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta)
sc.data <- GetAssayData(WATLAS.vat, "count")
sc.data <- sc.data[ rownames(sc.data) %in% rownames(counts),]
sc.eset <- Biobase::ExpressionSet(assayData = as.matrix(sc.data), phenoData = sc.pdata)

## Setup bulk RNA-seq data
bulk.eset <- ExpressionSet(assayData = counts)

## Deconvolute
decon.vat <- ReferenceBasedDecomposition(bulk.eset = bulk.eset, sc.eset = sc.eset, markers = markers, use.overlap = FALSE)
saveRDS(decon.vat, "data/Adipose/decon/decon_VAT.Rds")

## Plots and stats
plot(rowMeans(decon.vat$bulk.props),rowMeans(decon.vat$sc.props), las = 1, xlab="Average deconvoluted fraction", ylab="Average scRNA-seq fraction", pch = 16, main = "VAT")
cor(rowMeans(decon.vat$bulk.props),rowMeans(decon.vat$sc.props)) # 0.9984365
sqrt(mean((rowMeans(decon.vat$bulk.props)-rowMeans(decon.vat$sc.props))^2)) # 0.006805062

### Analyze the combined deconvolution results
## Combine the datasets from VAT and SAT
decon.sat.bulk <- t(decon.sat$bulk.props)
decon.vat.bulk <- t(decon.vat$bulk.props)
decon.vat.bulk <- decon.vat.bulk[, colnames(decon.vat.bulk) %in% colnames(decon.sat.bulk)]
decon.vat.bulk <- decon.vat.bulk[,match(colnames(decon.sat.bulk), colnames(decon.vat.bulk))]
decon.bulk <- rbind(decon.sat.bulk, decon.vat.bulk)
decon.bulk <- as.data.frame(decon.bulk)
md <- as.data.frame(colData(adipose))
md <- md[ match(rownames(decon.bulk), rownames(md)),]
decon.bulk$Sex <- factor(md$gtex.sex)
decon.bulk$Age <- factor(md$gtex.age)
decon.bulk$Depot <- factor(ifelse(grepl("Subcuta", md$gtex.smtsd),"SAT","VAT"))
saveRDS(decon.bulk, "data/Adipose/decon/decon_bulk.Rds")

## Scale each row to one
decon.bulk[,1:16] <- decon.bulk[,1:16] / rowSums(decon.bulk[,1:16])

## Differential abundance
results_bulk <- data.frame(Celltype = "NA", Variable = "NA", Difference = 0, ci.min = 0, ci.max = 0)
counter <- 1
for (i in 1:16) {
  model <- lm(decon.bulk[,i] ~ decon.bulk$Depot + decon.bulk$Sex)
  out <- eff_size(emmeans(model, specs = "Sex"), sigma = sigma(model), edf = df.residual(model))
  results_bulk[counter,1] <- colnames(decon.bulk)[i]
  results_bulk[counter,2] <- "Gender"
  results_bulk[counter,3] <- summary(out)$effect.size
  results_bulk[counter,4] <- summary(out)$lower.CL
  results_bulk[counter,5] <- summary(out)$upper.CL
  counter <- counter + 1
  out <- eff_size(emmeans(model, specs = "Depot"), sigma = sigma(model), edf = df.residual(model))
  results_bulk[counter,1] <- colnames(decon.bulk)[i]
  results_bulk[counter,2] <- "Depot"
  results_bulk[counter,3] <- summary(out)$effect.size
  results_bulk[counter,4] <- summary(out)$lower.CL
  results_bulk[counter,5] <- summary(out)$upper.CL
  counter <- counter + 1
}

## Plots
ggplot(data = results_bulk[ results_bulk$Variable == "Depot",], aes(x = Difference, y = Celltype, xmin = ci.min, xmax = ci.max)) + geom_point() + geom_pointrange() + theme_minimal()
ggplot(data = results_bulk[ results_bulk$Variable == "Gender",], aes(x = Difference, y = Celltype, xmin = ci.min, xmax = ci.max)) + geom_point() + geom_pointrange() + theme_minimal()
