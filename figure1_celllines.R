### Setup
## Options
options(future.globals.maxSize= 20000*1024^2)

## Dependencies
source("dependencies.R")

## Libraries
library(Seurat)
library(JOINTLY)
library(inflection)
library(UCell)
library(HGC)

## Functions
source("functions.R")

## Generate batches and noise for cell line data
celllines <- readRDS("data/CellLine/CellLines.rds")
noise <- simulateNoise(celllines, label_var = "author_label", mean = 0.5)
saveRDS(noise$no, "data/CellLines_no_noise.rds")
saveRDS(noise$noise, "data/CellLine/CellLines_noise.rds")
rm(list = c("celllines","noise"))

## Process dataset without noise
# Import
dataset <- readRDS("data/CellLine/CellLines_no_noise.rds")

# Run JOINTLY
if (file.exists("data/CellLine/CellLines_no_noise_JOINTLY.rds")) {
  solved <- readRDS("data/CellLine/CellLines_no_noise_JOINTLY.rds")
} else {
  proc <- preprocess(dataset, "batch")
  cpca <- cpca(proc)
  inputs <- prepareData(cpca$cpca)
  solved.list <- list()
  for (rep in 1:5) {
    solved.list[[rep]] <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca, bpparam = BiocParallel::MulticoreParam(), share.objects = FALSE)
  }

  # Choose best clustering
  results <- as.data.frame(matrix(ncol=4, nrow=1))
  for (rep in 1:5) {
    # Insert embedding
    H <- t(do.call("cbind",solved.list[[rep]]$Hmat))
    H <- H[ match(colnames(dataset), rownames(H)),]
    H <- scale(H)
    H <- t(scale(t(H)))
    colnames(H) <- paste("JOINTLY", 1:ncol(H), sep="_")
    dataset[["jointly"]] <- CreateDimReducObject(as.matrix(H), assay = "RNA")
    
    # Clustering
    dend <- HGC.dendrogram(G = HGC::SNN.Construction(H))
    clusters <- cutree(dend, k = 8)
    names(clusters) <- colnames(dataset)
    dataset$clusters <- clusters
    
    # Metrics
    results[rep,1] <- rep
    results[rep,2] <- aricode::ARI(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.999819
    results[rep,3] <- aricode::NMI(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.9996423
    results[rep,4] <- clevr::v_measure(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.9996565
  }
  solved <- solved.list[[5]]
  saveRDS(solved, "data/CellLine/CellLines_no_noise_JOINTLY.rds")
}

# Insert best embedding
H <- t(do.call("cbind",solved$Hmat))
H <- H[ match(colnames(dataset), rownames(H)),]
H <- scale(H)
H <- t(scale(t(H)))
colnames(H) <- paste("JOINTLY", 1:ncol(H), sep="_")
dataset[["jointly"]] <- CreateDimReducObject(as.matrix(H), assay = "RNA")

# Visualize
dataset <- NormalizeData(dataset)
dataset <- RunUMAP(dataset, reduction = "jointly", dims = 1:15)

# Clustering
dend <- HGC.dendrogram(G = HGC::SNN.Construction(H))
clusters <- cutree(dend, k = 8)
names(clusters) <- colnames(dataset)
dataset$clusters <- clusters

# Metrics
aricode::ARI(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.999819
aricode::NMI(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.9996423
clevr::v_measure(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.9996565

# Calculate module scores
modules <- getModules(solved)
dataset <- UCell::AddModuleScore_UCell(dataset, features = modules)

# Default PCA using the same features
VariableFeatures(dataset) <- cpca$features
dataset <- NormalizeData(dataset)
dataset <- ScaleData(dataset)
dataset <- RunPCA(dataset)
dataset <- RunUMAP(dataset, reduction = "pca", dims = 1:15, reduction.name = "Default_UMAP")

# Clustering on default PCA
dend <- HGC.dendrogram(G = HGC::SNN.Construction(dataset@reductions$pca@cell.embeddings[,1:15]))
clusters <- cutree(dend, k = 8)
names(clusters) <- colnames(dataset)
dataset$default_clusters <- clusters

# Metrics
aricode::ARI(dataset@meta.data[ ,"default_clusters"], dataset@meta.data[ ,"author_label"]) # 0.9998022
aricode::NMI(dataset@meta.data[ ,"default_clusters"], dataset@meta.data[ ,"author_label"]) # 0.9996399
clevr::v_measure(dataset@meta.data[ ,"default_clusters"], dataset@meta.data[ ,"author_label"]) # 0.9996522

# Save the object
saveRDS(dataset, "data/CellLine/CellLines_no_noise_result.Rds")

# Plots
DimPlot(dataset, group.by = "author_label", label = TRUE, reduction = "Default_UMAP") + DimPlot(dataset, group.by = "default_clusters", label = TRUE, reduction = "Default_UMAP") + DimPlot(dataset, group.by = "author_label", label = TRUE) + DimPlot(dataset, group.by = "clusters", label = TRUE) & NoLegend()
DimPlot(dataset, group.by = "batch", label = TRUE, reduction = "Default_UMAP") + DimPlot(dataset, group.by = "batch", label = TRUE) & NoLegend()
FeaturePlot(dataset, c("factor_1_UCell", "factor_2_UCell","factor_3_UCell","factor_4_UCell","factor_5_UCell","factor_7_UCell","factor_9_UCell","factor_12_UCell"), min.cutoff = "q80")
FeaturePlot(dataset, c("factor_1_UCell", "factor_2_UCell","factor_3_UCell","factor_4_UCell","factor_5_UCell","factor_7_UCell","factor_9_UCell","factor_12_UCell"))

## Process dataset with noise
# Import
dataset <- readRDS("data/CellLine/CellLines_noise.rds")

# Run JOINTLY
if (file.exists("data/CellLine/CellLines_noise_JOINTLY.rds")) {
  solved <- readRDS("data/CellLine/CellLines_noise_JOINTLY.rds")
} else {
  proc <- preprocess(dataset, "batch")
  cpca <- cpca(proc)
  inputs <- prepareData(cpca$cpca)
  solved.list <- list()
  for (rep in 1:5) {
    solved.list[[rep]] <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca, bpparam = BiocParallel::MulticoreParam(), share.objects = FALSE)
  }
  
  # Choose best clustering
  results <- as.data.frame(matrix(ncol=4, nrow=1))
  for (rep in 1:5) {
    # Insert embedding
    H <- t(do.call("cbind",solved.list[[rep]]$Hmat))
    H <- H[ match(colnames(dataset), rownames(H)),]
    H <- scale(H)
    H <- t(scale(t(H)))
    colnames(H) <- paste("JOINTLY", 1:ncol(H), sep="_")
    dataset[["jointly"]] <- CreateDimReducObject(as.matrix(H), assay = "RNA")
    
    # Clustering
    dend <- HGC.dendrogram(G = HGC::SNN.Construction(H))
    clusters <- cutree(dend, k = 8)
    names(clusters) <- colnames(dataset)
    dataset$clusters <- clusters
    
    # Metrics
    results[rep,1] <- rep
    results[rep,2] <- aricode::ARI(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.999819
    results[rep,3] <- aricode::NMI(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.9996423
    results[rep,4] <- clevr::v_measure(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.9996565
  }
  solved <- solved.list[[2]]
  saveRDS(solved, "data/CellLine/CellLines_noise_JOINTLY.rds")
}

# Insert best embedding
H <- t(do.call("cbind",solved$Hmat))
H <- H[ match(colnames(dataset), rownames(H)),]
H <- scale(H)
H <- t(scale(t(H)))
colnames(H) <- paste("JOINTLY", 1:ncol(H), sep="_")
dataset[["jointly"]] <- CreateDimReducObject(as.matrix(H), assay = "RNA")

# Visualize
dataset <- NormalizeData(dataset)
dataset <- RunUMAP(dataset, reduction = "jointly", dims = 1:15)

# Clustering
dend <- HGC.dendrogram(G = HGC::SNN.Construction(H))
clusters <- cutree(dend, k = 8)
names(clusters) <- colnames(dataset)
dataset$clusters <- clusters

# Metrics
aricode::ARI(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.998979
aricode::NMI(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.9982883
clevr::v_measure(dataset@meta.data[ ,"clusters"], dataset@meta.data[ ,"author_label"]) # 0.998317

# Calculate scores
modules <- getModules(solved)
dataset <- UCell::AddModuleScore_UCell(dataset, features = modules)

## Default PCA
VariableFeatures(dataset) <- cpca$features
dataset <- NormalizeData(dataset)
dataset <- ScaleData(dataset)
dataset <- RunPCA(dataset)
dataset <- RunUMAP(dataset, reduction = "pca", dims = 1:15, reduction.name = "Default_UMAP")

# Clustering on default PCA
dend <- HGC.dendrogram(G = HGC::SNN.Construction(dataset@reductions$pca@cell.embeddings[,1:15]))
clusters <- cutree(dend, k = 8)
names(clusters) <- colnames(dataset)
dataset$default_clusters <- clusters

# Metrics
aricode::ARI(dataset@meta.data[ ,"default_clusters"], dataset@meta.data[ ,"author_label"]) # 0.8853705
aricode::NMI(dataset@meta.data[ ,"default_clusters"], dataset@meta.data[ ,"author_label"]) # 0.9442867
clevr::v_measure(dataset@meta.data[ ,"default_clusters"], dataset@meta.data[ ,"author_label"]) # 0.95545161

# Save the object
saveRDS(dataset, "data/CellLine/CellLines_noise_result.Rds")

# Plots
DimPlot(dataset, group.by = "author_label", label = TRUE, reduction = "Default_UMAP") + DimPlot(dataset, group.by = "default_clusters", label = TRUE, reduction = "Default_UMAP") + DimPlot(dataset, group.by = "author_label", label = TRUE) + DimPlot(dataset, group.by = "clusters", label = TRUE) & NoLegend()
DimPlot(dataset, group.by = "batch", label = TRUE, reduction = "Default_UMAP") + DimPlot(dataset, group.by = "batch", label = TRUE) & NoLegend()
FeaturePlot(dataset, c("factor_1_UCell", "factor_2_UCell", "factor_6_UCell","factor_15_UCell","factor_7_UCell", "factor_10_UCell","factor_12_UCell","factor_14_UCell"), min.cutoff = "q80")
FeaturePlot(dataset, c("factor_1_UCell", "factor_2_UCell", "factor_6_UCell","factor_15_UCell","factor_7_UCell", "factor_10_UCell","factor_12_UCell","factor_14_UCell"))

### Analysis of W matrices
# Import solution from noisy datasets
solved <- readRDS("data/CellLine/CellLines_noise_JOINTLY.rds")

# Intepretation
W <- solved$Wmat
for (i in 1:length(W)) { 
  W.tmp <- W[[i]]
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W[[i]] <- W.tmp
}

# Correlation between batches
batch1_self <- cor(W[[1]], W[[1]])
rownames(batch1_self) <- colnames(batch1_self) <- paste("Batch1_Factor", 1:15, sep="")
batch1_self <- batch1_self[hclust(dist(batch1_self))$order,hclust(dist(batch1_self))$order]
batch1_htm <- Heatmap(batch1_self, cluster_columns = FALSE, cluster_rows = FALSE, name = "Batch1")

batch2_self <- cor(W[[2]], W[[2]])
rownames(batch2_self) <- colnames(batch2_self) <- paste("Batch2_Factor", 1:15, sep="")
batch2_self <- batch2_self[hclust(dist(batch2_self))$order,hclust(dist(batch2_self))$order]
batch2_htm <- Heatmap(batch2_self, cluster_columns = FALSE, cluster_rows = FALSE, name = "Batch2")

batch_cross <- cor(W[[1]], W[[2]])
rownames(batch_cross) <- paste("Batch1_Factor", 1:15, sep="")
colnames(batch_cross) <- paste("Batch2_Factor", 1:15, sep="")
batch_cross <- batch_cross[hclust(dist(batch_cross))$order,hclust(dist(batch_cross))$order]
batch_cross_htm <- Heatmap(batch_cross, cluster_columns = FALSE, cluster_rows = FALSE, name = "Cross-batch")

# Heatmap
batch1_htm + batch_cross_htm + batch2_htm

## Comparison of modules with and without batch effects
## Get modules with noise
noise_modules <- getModules(solved)

## Get modules without noise
# Import solution from noisy datasets
solved <- readRDS("data/CellLine/CellLines_no_noise_JOINTLY.rds")
clean_modules <- getModules(solved)

## Calculate Jaccard index
jacc <- do.call("rbind", lapply(clean_modules, FUN = function(y) { unlist(lapply(noise_modules, y = y, FUN = function(x,y) { sum(x %in% y) / ( length(x) + length(y) - sum(x %in% y)) })) }))
jacc <- jacc[apply(jacc, 2, FUN="which.max"),]
rownames(jacc) <- paste("Clean_", rownames(jacc), sep="")
colnames(jacc) <- paste("Noise_", colnames(jacc), sep="")

# Heatmap
Heatmap(as.matrix(jacc), cluster_rows = FALSE, cluster_columns = FALSE, col = circlize::colorRamp2(c(0, 0.33, 0.66, 1), c("white", RColorBrewer::brewer.pal(3, "Blues"))))