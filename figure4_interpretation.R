### Setup
## Options
options(future.globals.maxSize= 20000*1024^2)

## Dependencies
source("dependencies.R")

## Libraries
library(Seurat)
library(JOINTLY)
library(BiocParallel)
library(inflection)
library(UCell)
library(scCustomize)
library(enrichR)
library(ComplexHeatmap)
library(HGC)
library(presto)

## Functions
source("functions.R")

## Import pancreas data
dataset <- readRDS("data/Adipose/raw/Emont/emont.rds")
dataset[["SCT"]] <- NULL
dataset[["integrated"]] <- NULL
dataset <- subset(dataset, tissue == "adipose")
dataset <- subset(dataset, depot == "VAT")
dataset$individual <- as.character(dataset$individual)

## Run JOINTLY
proc <- preprocess(dataset, "individual")
cpca <- cpca(proc, bpparam = MulticoreParam())
inputs <- prepareData(cpca$cpca)
solved <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca, bpparam = BiocParallel::MulticoreParam())
saveRDS(solved, "data/Interpretation/jointly.Rds")


## Embed using the H matrix
H <- solved$Hmat.scaled
colnames(H) <- paste("JOINTLY", 1:ncol(H), sep="_")
dataset[["jointly"]] <- CreateDimReducObject(as.matrix(H), assay = "RNA")
dataset <- RunUMAP(dataset, reduction = "jointly", dims = 1:15)

## Interpret W matrix
# Scale the factors
W <- solved$Wmat
for (i in 1:length(W)) { 
  W.tmp <- W[[i]]
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W[[i]] <- W.tmp
}

# Sum the W matrix
W.sum <- W[[1]]
for (i in 2:length(W)) { W.sum <- W.sum + W[[i]]}
W.sum <- W.sum / length(solved$Wmat)
rownames(W.sum) <- rownames(cpca$normalized[[1]])

# Scale the sum matrix
W.tmp <- W.sum
W.tmp <- scale(W.tmp)
W.tmp <- t(W.tmp)
W.tmp <- scale(W.tmp)
W.tmp <- t(W.tmp)
W.sum <- W.tmp

# Get modules
modules <- list()
for (i in 1:15) {
  modules[[length(modules)+1]] <- names(sort(W.sum[,i], decreasing = TRUE))[1:inflection::uik(y = sort(W.sum[,i], decreasing = TRUE), x = seq(1, nrow(W.sum),1))]
  names(modules)[length(modules)] <- paste("factor_", i, sep="")
}
saveRDS(modules, "data/Interpretation/modules.Rds")

# Calculate scores
dataset <- NormalizeData(dataset)
dataset <- UCell::AddModuleScore_UCell(dataset, features = modules)
saveRDS(dataset, "/work/NMF_project/Final_datasets/Interpretation/dataset.Rds")

## Reimport data
dataset <- readRDS("data/Interpretation/dataset.Rds")
modules <- readRDS("data/Interpretation/modules.Rds")

# Cluster the data
dend <- HGC.dendrogram(G = HGC::SNN.Construction(dataset@reductions$jointly@cell.embeddings))
clusters <- cutree(dend, k = 7)
names(clusters) <- colnames(dataset)
dataset$clusters <- clusters

dataset$label <- "NULL"
dataset@meta.data[ dataset@meta.data$clusters == 1, "label"] <- "FAP"
dataset@meta.data[ dataset@meta.data$clusters == 2, "label"] <- "Mesothelial"
dataset@meta.data[ dataset@meta.data$clusters == 3, "label"] <- "Mural"
dataset@meta.data[ dataset@meta.data$clusters == 4, "label"] <- "Adipocyte"
dataset@meta.data[ dataset@meta.data$clusters == 5, "label"] <- "Immune"
dataset@meta.data[ dataset@meta.data$clusters == 6, "label"] <- "Endothelial"
dataset@meta.data[ dataset@meta.data$clusters == 7, "label"] <- "LEC"

# Save the object
saveRDS(dataset, "data/Interpretation/adipose.rds")

# Subcluster immune cells
dataset.subset <- subset(dataset, label == "Immune")

# Run JOINTLY
proc <- preprocess(dataset.subset, "individual")
cpca <- cpca(proc, bpparam = MulticoreParam())
inputs <- prepareData(cpca$cpca)
solve.list <- list()
solved <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca, bpparam = BiocParallel::MulticoreParam())
saveRDS(solved, "data/Interpretation/jointly_immune.Rds")

## Embed using the H matrix
H <- t(do.call("cbind",solved$Hmat))
H <- H[ match(colnames(dataset.subset), rownames(H)),]
H <- scale(H)
H <- t(scale(t(H)))
colnames(H) <- paste("JOINTLY", 1:ncol(H), sep="_")
dataset.subset[["jointly"]] <- CreateDimReducObject(as.matrix(H), assay = "RNA")
dataset.subset <- RunUMAP(dataset.subset, reduction = "jointly", dims = 1:15)
dend <- HGC.dendrogram(G = HGC::SNN.Construction(dataset.subset@reductions$jointly@cell.embeddings))
clusters <- cutree(dend, k = 100)
names(clusters) <- colnames(dataset.subset)
dataset.subset$clusters <- clusters

dataset.subset$label <- "NULL"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("20"), "label"] <- "Bcell"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("18"), "label"] <- "DC"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("41", "62"), "label"] <- "DC"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("26", "44", "3", "10"), "label"] <- "Monocyte"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("51","15"), "label"] <- "Mast"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("68","84", "100", "79"), "label"] <- "NK"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("53","90","57","96","81","92","59","8","11","52","88","82","72","47","34","76","60","97","78","74","89","99","87","19"), "label"] <- "Tcell"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("49","4","23","5","14","93","80","65","98","40","73","67","37","58","70","32","63","42","64","24","45","75","30","28","43","6","12","13","16","27","94","83","61","21","49","66"), "label"] <- "Macrophage"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("33","46"), "label"] <- "Macrophage"
dataset.subset@meta.data[ dataset.subset@meta.data$clusters %in% c("1","2","7","9","17","22","25","29","31","35","36","38","39","48","50","54","55","56","69","71","77","85","86","91","95"), "label"] <- "Macrophage"

## Transfer labels to main object
for (lbl in unique(dataset.subset$label)) {
  dataset@meta.data[ rownames(dataset@meta.data) %in% rownames(dataset.subset@meta.data[dataset.subset@meta.data$label == lbl,]), "label"] <- lbl
}

# Save the object
saveRDS(dataset, "data/Interpretation/adipose.rds")

## UMAP
DimPlot(dataset, group.by="label", label = TRUE, repel = TRUE)

## Integration
DimPlot(dataset, split.by= "individual")

## Markers
DotPlot_scCustom(dataset, features = c("F13A1","SIGLEC1","CSF1R","FCN1","S100A9","COTL1","THEMIS","IL7R","BCL11B","CD1C","FCER1A","CLEC10A","MS4A1","BCL11A","SEL1L3","PDGFRA", "C7", "COL6A3", "PPARG", "PLIN1", "GPAM", "VWF", "PECAM1", "MECOM", "CPA3", "KIT","MS4A2", "GNLY", "CD247", "NKG7", "ACTA2", "MYH11", "MYO1B", "BNC1", "UPK3B", "MSLN", "PROX1","PDPN","FLT4"), group.by="label", x_lab_rotate = TRUE, flip_axes = TRUE, cluster.idents = TRUE)

## Factors
agg <- aggregate(dataset@meta.data[,c("factor_1_UCell","factor_2_UCell","factor_3_UCell","factor_4_UCell","factor_5_UCell","factor_6_UCell","factor_7_UCell","factor_8_UCell","factor_9_UCell","factor_10_UCell","factor_11_UCell","factor_12_UCell","factor_13_UCell","factor_14_UCell","factor_15_UCell")], by = list(dataset$label), FUN="median")
rownames(agg) <- agg[,1]
agg <- agg[,-1]
Heatmap(as.matrix(agg))

## Label using enrichR
dbs <- c("Azimuth_Cell_Types_2021", "CellMarker_Augmented_2021")

# Labelling
for (i in 1:15) {
  enriched <- enrichr(modules[[i]], dbs)
  enriched[[2]] <- enriched[[2]][grep("Adipose", enriched[[2]]$Term),]
  enriched[[2]] <- enriched[[2]][grep("Brown", invert = TRUE, enriched[[2]]$Term),]
  enriched[[2]] <- enriched[[2]][grep("Beige", invert = TRUE, enriched[[2]]$Term),]
  enriched[[2]] <- enriched[[2]][!is.na(enriched[[1]]$Term),]
  enriched[[1]] <- enriched[[1]][grep("Neuron", invert = TRUE, enriched[[1]]$Term),]
  enriched[[1]] <- enriched[[1]][grep("neuron", invert = TRUE, enriched[[1]]$Term),]
  enriched[[1]] <- enriched[[1]][!is.na(enriched[[1]]$Term),]
  enriched <- do.call("rbind", enriched)
  enriched <- enriched[ enriched$Adjusted.P.value <= 0.1,]
  enriched <- enriched[order(enriched$Combined.Score, decreasing = TRUE),]
  if (i == 1) {
    ct_ranks <- data.frame(factor = i, rank = seq(1,5,1), celltype = enriched[ 1:5,"Term"])
  } else {
    ct_ranks <- rbind(ct_ranks, data.frame(factor = i, rank = seq(1,5,1), celltype = enriched[ 1:5,"Term"]))
  }
}
colnames(ct_ranks)[3] <- "Predicted_CellType"
ct_ranks$Label <- "NA"
ct_ranks[ ct_ranks$factor == 1, "Label"] <- "NK/T"
ct_ranks[ ct_ranks$factor == 2, "Label"] <- "Mesothelial"
ct_ranks[ ct_ranks$factor == 3, "Label"] <- "Mesothelial"
ct_ranks[ ct_ranks$factor == 4, "Label"] <- "DC/Monocyte/Macrophage"
ct_ranks[ ct_ranks$factor == 5, "Label"] <- "Mural"
ct_ranks[ ct_ranks$factor == 6, "Label"] <- "FAP"
ct_ranks[ ct_ranks$factor == 7, "Label"] <- "LEC"
ct_ranks[ ct_ranks$factor == 8, "Label"] <- "Endothelial"
ct_ranks[ ct_ranks$factor == 9, "Label"] <- "Mesothelial"
ct_ranks[ ct_ranks$factor == 10, "Label"] <- "Adipocyte/FAP"
ct_ranks[ ct_ranks$factor == 11, "Label"] <- "Mesothelial"
ct_ranks[ ct_ranks$factor == 12, "Label"] <- "Endothelial"
ct_ranks[ ct_ranks$factor == 13, "Label"] <- "Adipocyte"
ct_ranks[ ct_ranks$factor == 14, "Label"] <- "Endothelial"
ct_ranks[ ct_ranks$factor == 15, "Label"] <- "DC/Monocyte/Macrophage"
ct_ranks$Hit <- 0
ct_ranks[c(1,6,14,16,21,30,31,39,42,49,60,63,67,71),"Hit"] <- 1
write.table(ct_ranks, "data/Interpretation/celltype_ranks.txt", quote = F, row.names = FALSE, col.names = FALSE)

## Highlight example of factors
# Feature Plot
FeaturePlot(dataset, "factor_10_UCell", min.cutoff = "q80")

# Factor 10 is associated to insulin signaling + PI3K-Akt
dbs <- c("WikiPathway_2021_Human", "KEGG_2021_Human","Reactome_2022")
enriched <- enrichr(modules[[10]], dbs)
enriched <- do.call("rbind", enriched)
enriched$FDR <- p.adjust(enriched$P.value, method = "fdr")
enriched <- enriched[ enriched$FDR <= 0.05,]
enriched <- enriched[order(enriched$FDR, decreasing = FALSE),]
enriched <- enriched[c(1,8,18),]
barplot(-log10(enriched$FDR), las=2, ylab = "-log10 FDR", names = enriched$Term)

# Compare genders
par(mfcol=c(1,2))
adipocytes <- subset(dataset, label == "Adipocyte")@meta.data
FAP <- subset(dataset, label == "FAP")@meta.data

adipocytes[ adipocytes$factor_10_UCell >= quantile(adipocytes$factor_10_UCell, quant), "positive"] <- 1
agg <- aggregate(adipocytes$factor_10_UCell, by = list(adipocytes$individual), FUN="mean")
freq <- agg$x
names(freq) <- agg$Group.1
md <- adipocytes[ duplicated(adipocytes$individual)==FALSE,]
md <- md[ match(names(freq), md$individual),]
md$freq <- freq
boxplot(freq ~ sex, data = md, ylab="Mean factor 10 signal", names = c("Female","Male"), xlab="Gender", las = 1, main = "Adipocyte")

FAP[ FAP$factor_10_UCell >= quantile(adipocytes$factor_10_UCell, quant), "positive"] <- 1
agg <- aggregate(FAP$factor_10_UCell, by = list(FAP$individual), FUN="mean")
freq <- agg$x
names(freq) <- agg$Group.1
md <- FAP[ duplicated(FAP$individual)==FALSE,]
md <- md[ match(names(freq), md$individual),]
md$freq <- freq
boxplot(freq ~ sex, data = md, ylab="Mean factor 10 signal", names = c("Female","Male"), xlab="Gender", las = 1, main = "FAP")

#### Compare interpretable factors with DE and between clusters
# Overall test
DE <- wilcoxauc(dataset, group_by = "label")
DE.list <- split(DE, f = DE$group)

# Batch consistency
consistency <- c()
for (ds in unique(dataset$individual)) {
  DE.batch <- wilcoxauc(subset(dataset, individual == ds), group_by = "label")
  DE.batch <- DE.batch[ DE.batch$logFC > 0 & DE.batch$padj <= 0.05,]
  consistency <- c(consistency, paste(DE.batch$group, DE.batch$feature, sep="_"))
}
consistency <- as.data.frame(table(consistency))
consistency$group <- substr(consistency$consistency, 0, regexpr("_", consistency$consistency)-1)
consistency$feature <- substr(consistency$consistency, regexpr("_", consistency$consistency)+1, nchar(as.character(consistency$consistency)))
consistency$Freq <- consistency$Freq / length(unique(dataset$individual))
for (lbl in unique(dataset$label)) {
  idx <- which(names(DE.list) == lbl)
  cons <- consistency[ consistency$group == lbl,]
  cons <- cons[ match(DE.list[[idx]]$feature, cons$feature),]
  DE.list[[idx]]$batch_consistency <- cons$Freq
}
DE.list <- lapply(DE.list, FUN = function(x) { x$marker <- ifelse(x$logFC > 0 & x$padj <= 0.05, "yes", "no"); return(x) })
DE.list <- lapply(DE.list, FUN = function(x) { x$module <- ifelse(x$feature %in% unlist(modules), "in", "out"); return(x) } )
DE.list <- lapply(DE.list, FUN = function(x) { x$batch_consistency <- ifelse(is.na(x$batch_consistency), 0, x$batch_consistency); return(x)})

## Overlap between modules and markers
# Calculate Jaccard index
jaccard.matrix <- as.data.frame(matrix(ncol=13, nrow=15))
for (m in 1:15) {
  for (j in 1:13) {
    feats <- DE.list[[j]][ DE.list[[j]]$marker == "yes" & DE.list[[j]]$module == "in" ,"feature"]
    jaccard.matrix[m,j] <- sum(feats %in% modules[[m]]) / (length(modules[[m]]) + length(feats) - sum(feats %in% modules[[m]]))
    colnames(jaccard.matrix)[j] <- unique(DE.list[[j]]$group)
    rownames(jaccard.matrix)[m] <- names(modules)[m]
  }
}

# Plot it
Heatmap(jaccard.matrix, cluster_columns = FALSE, cluster_rows = FALSE)

# Calculate average AUCs
marker.module.auc <- marker.notmodule.auc <- notmarker.module.auc <- notmarker.notmodule.auc <- marker.module.cons <- marker.notmodule.cons <- c() 
for (m in 1:13) {
  j <- which.max(jaccard.matrix[,m])
  marker.module.auc <- c(marker.module.auc, median(DE.list[[m]][ DE.list[[m]]$marker == "yes" & DE.list[[m]]$feature %in% modules[[j]], "auc"]))
  marker.notmodule.auc <- c(marker.notmodule.auc, median(DE.list[[m]][ DE.list[[m]]$marker == "yes" & !(DE.list[[m]]$feature %in% modules[[j]]), "auc"]))
  notmarker.module.auc <- c(notmarker.module.auc, median(DE.list[[m]][ DE.list[[m]]$marker == "no" & DE.list[[m]]$feature %in% unlist(modules[seq(1,15,1)[-j]]), "auc"]))
  notmarker.notmodule.auc <- c(notmarker.notmodule.auc, median(DE.list[[m]][ DE.list[[m]]$marker == "no" & !(DE.list[[m]]$feature %in% unlist(modules[seq(1,15,1)[-j]])), "auc"]))
  marker.module.cons <- c(marker.module.cons, median(DE.list[[m]][ DE.list[[m]]$marker == "yes" & DE.list[[m]]$feature %in% modules[[j]], "batch_consistency"]))
  marker.notmodule.cons <- c(marker.notmodule.cons, median(DE.list[[m]][ DE.list[[m]]$marker == "yes" & !(DE.list[[m]]$feature %in% modules[[j]]), "batch_consistency"]))
}

# Plot it
par(mfcol=c(1,2))
boxplot(marker.module.auc, marker.notmodule.auc, notmarker.module.auc, notmarker.notmodule.auc, las = 1, ylab = "Median AUC per cell type", names = c("A","B","C","D"))
boxplot(marker.module.cons, marker.notmodule.cons, las = 1, ylab = "Median AUC per cell type", names = c("A","B"), ylim=c(0,1))
dev.off()

### Time-courses
path <- "data/Differentiation/gastrulation.h5ad"
root <- rhdf5::H5Fopen(path)
target <- "/X"
i_path <- paste0(target,"/indices")
p_path <- paste0(target,"/indptr")
x_path <- paste0(target,"/data")
i <- as.vector(unlist(rhdf5::h5read(root, i_path)))
p <- as.vector(unlist(rhdf5::h5read(root, p_path)))
x <- as.vector(unlist(rhdf5::h5read(root, x_path)))
o <- as.vector(rhdf5::h5read(root, "/obs")$"_index")
v <- as.vector(rhdf5::h5read(root, "/var")$index)
lt <- as.vector(rhdf5::h5read(root, "/obs")$latent_time)
ps <- as.vector(rhdf5::h5read(root, "/obs")$dynamical_velocity_pseudotime)
cl <- as.vector(rhdf5::h5read(root, "/obs")$celltype)
sm <- as.vector(rhdf5::h5read(root, "/obs")$sample)
stage <- as.vector(rhdf5::h5read(root, "/obs")$stage)
batch <- as.vector(rhdf5::h5read(root, "/obs")$sequencing.batch)
theiler <- as.vector(rhdf5::h5read(root, "/obs")$theiler)
pca <- t(rhdf5::h5read(root, "/obsm")$X_pca)
umap <- t(rhdf5::h5read(root, "/obsm")$X_umap)
dims <- c(length(v), length(o))
rhdf5::H5Fclose(root)
m <- Matrix::sparseMatrix(i = i, p = p, x = x, index1 = FALSE, dims = dims)
rownames(m) <- gsub("_", "-", v)
colnames(m) <- o
seu <- CreateSeuratObject(m, min.cells = 0, min.features = 0)
seu$latent_time <- lt
seu$pseudo_time <- ps
seu$celltype <- factor(cl$codes, labels = cl$categories)
seu$sample <- factor(sm$codes, labels = sm$categories)
seu$stage <- factor(stage$codes, labels = stage$categories)
seu$sequencing.batch <- factor(batch$codes, labels = batch$categories)
seu$theiler <- factor(theiler$codes, labels = theiler$categories)
rownames(umap) <- rownames(pca) <- colnames(seu)
colnames(umap) <- c("UMAP_1","UMAP_2")
colnames(pca) <- paste("PCA_", 1:50, sep="")
seu[["pca"]] <- CreateDimReducObject(pca)
seu[["umap"]] <- CreateDimReducObject(umap)

## JOINTLY
proc <- preprocess(seu, "sequencing.batch", verbose = FALSE)
cpca <- cpca(proc, verbose = FALSE, nfeat = 2000)
inputs <- prepareData(cpca$cpca, verbose = FALSE)
for (rep in 1:5) {
  # Embed and evaluate
  solved <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca, bpparam = BiocParallel::MulticoreParam(), progressbar = TRUE, share.objects = FALSE, verbose = TRUE)
  H <- t(do.call("cbind",solved$Hmat))
  H <- H[ match(colnames(seu), rownames(H)),]
  H <- scale(H)
  H <- t(scale(t(H)))
  colnames(H) <- paste("JOINTLY_rep", rep,"_",1:ncol(H), sep="")
  seu[[paste("JOINTLY_rep",rep, sep="")]] <- CreateDimReducObject(H, assay = "RNA")
  rm(H)
}
seu <- RunUMAP(seu, reduction = "JOINTLY_rep5", dims = 1:15, reduction.name = "JOINTLY_umap")
saveRDS(seu, "data/Differentiation/gastrulation.Rds")

# Scale the factors
W <- solved$Wmat
for (i in 1:length(W)) { 
  W.tmp <- W[[i]]
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W[[i]] <- W.tmp
}
W.sum <- W[[1]]
for (i in 2:length(W)) { W.sum <- W.sum + W[[i]]}
W.sum <- W.sum / length(solved$Wmat)
rownames(W.sum) <- rownames(cpca$normalized[[1]])

# Scale the sum matrix
W.tmp <- W.sum
W.tmp <- scale(W.tmp)
W.tmp <- t(W.tmp)
W.tmp <- scale(W.tmp)
W.tmp <- t(W.tmp)
W.sum <- W.tmp

# Get modules
modules <- list()
for (i in 1:15) {
  modules[[length(modules)+1]] <- names(sort(W.sum[,i], decreasing = TRUE))[1:inflection::uik(y = sort(W.sum[,i], decreasing = TRUE), x = seq(1, nrow(W.sum),1))]
  names(modules)[length(modules)] <- paste("factor_", i, sep="")
}
saveRDS(modules, "data/Differentiation/modules.Rds")

# Calculate scores
seu <- NormalizeData(seu)
seu <- UCell::AddModuleScore_UCell(seu, features = modules)
saveRDS(seu, "data/Interpretation/gastrulation.Rds")

FeaturePlot(seu, reduction = "JOINTLY_umap", "latent_time") + DimPlot(seu, reduction = "JOINTLY_umap", group.by="celltype") + FeaturePlot(seu, reduction = "umap", "latent_time") + DimPlot(seu, reduction = "umap", group.by="celltype")

agg <- aggregate(seu@meta.data[,11:25], by =list(cut(seu@meta.data$latent_time, breaks = 50)), FUN="mean")
for (i in 2:ncol(agg)) { agg[,i] <- (agg[,i] - max(agg[,i])) / (min(agg[,i]) - max(agg[,i])) }
pdf("/work/NMF_project/Figures/V4/Differentiation_Heatmap.pdf", width = 10, height = 10, useDingbats = FALSE)
Heatmap(t(agg[,c(12,2,10,14,5,8,9,4,6,3,7,16,11,13,15)]), cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

early <- unique(unlist(modules[c(11,1,9,13)]))
mid_early <- unique(unlist(modules[c(4,7,8)]))
mid_late <- unique(unlist(modules[c(3)]))
late <- unique(unlist(modules[c(2,5,6)]))
temp <- unique(unlist(modules[c(15,10,12,14)]))

dbs <- c("Reactome_2022", "WikiPathways_2019_Mouse", "KEGG_2019_Mouse", "GO_Biological_Process_2023","Azimuth_Cell_Types_2021", "CellMarker_Augmented_2021")

early_enriched <- enrichr(early, dbs)
mid_early_enriched <- enrichr(mid_early, dbs)
mid_late_enriched <- enrichr(mid_late, dbs)
late_enriched <- enrichr(late, dbs)
temp_enriched <- enrichr(temp, dbs)

# Extract
par(mfcol=c(1,5))
term <- "Insulin Signaling WP65"
idx <- 2
barplot(c(
  ifelse(identical(numeric(0), early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),0,early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),0,temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score),0,late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score)), las = 1, ylab="Combined score", main = "Insulin Signaling", names = c("Early","Temporary","Early mid", "Late mid","Late"))


term <- "IL-6 signaling Pathway WP387"
idx <- 2
barplot(c(
  ifelse(identical(numeric(0), early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),0,early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),0,temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score),0,late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score)), las = 1, ylab="Combined score", main = "IL-6 Signaling", names = c("Early","Temporary","Early mid", "Late mid","Late"))


term <- "Kit Receptor Signaling Pathway WP407"
idx <- 2
barplot(c(
  ifelse(identical(numeric(0), early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),0,early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),0,temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score),0,late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score)), las = 1, ylab="Combined score", main = "Kit Signaling", names = c("Early","Temporary","Early mid", "Late mid","Late"))


term <- "TGF Beta Signaling Pathway WP113"
idx <- 2
barplot(c(
  ifelse(identical(numeric(0), early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),0,early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),0,temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score),0,late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score)), las = 1, ylab="Combined score", main = "TGFb Signaling", names = c("Early","Temporary","Early mid", "Late mid","Late"))

term <- "Heme Biosynthesis WP18"
idx <- 2
barplot(c(
  ifelse(identical(numeric(0), early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),0,early_enriched[[idx]][ early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),0,temp_enriched[[idx]][ temp_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_early_enriched[[idx]][ mid_early_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),0,mid_late_enriched[[idx]][ mid_late_enriched[[idx]][,1] == term,]$Combined.Score),
  ifelse(identical(numeric(0), late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score),0,late_enriched[[idx]][ late_enriched[[idx]][,1] == term,]$Combined.Score)), las = 1, ylab="Combined score", main = "Heme Biosynthesis", names = c("Early","Temporary","Early mid", "Late mid","Late"))
