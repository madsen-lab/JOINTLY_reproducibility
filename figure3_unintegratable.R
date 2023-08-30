### Setup
## Options
options(future.globals.maxSize= 40000*1024^2)

## Dependencies
source("dependencies.R")

## Libraries
library(Seurat)
library(SeuratDisk)
library(BiocParallel)
library(enrichR)
library(ggplot2)
library(patchwork)
library(data.table)
library(HGC)
library(reticulate)
use_virtualenv("PythonEnv/scVI_new/")

## Functions
source("functions.R")

## Import datasets
pancreas <- readRDS("data/Pancreas/Human_Pancreas.rds")
pancreas$tissue <- "Pancreas"
lung <- readRDS("data/Lung/Human_Lung_sub.rds")
lung$tissue <- "Lung"

## Combine the datasets
pancreas.counts <- GetAssayData(pancreas, "count")
lung.counts <- GetAssayData(lung, "count")
common <- intersect(rownames(pancreas.counts), rownames(lung.counts))
pancreas.counts <- pancreas.counts[ rownames(pancreas.counts) %in% common,]
lung.counts <- lung.counts[ rownames(lung.counts) %in% common,]
lung.counts <- lung.counts[ match(rownames(pancreas.counts), rownames(lung.counts)),]
pancreas[["RNA"]] <- CreateAssayObject(counts = pancreas.counts)
lung[["RNA"]] <- CreateAssayObject(counts = lung.counts)
mix <- merge(pancreas, lung)

## Prepare the datasets
dataset <- mix
dataset$batch_label <- as.character(dataset$batch_label)
assays.remove <- names(dataset@assays)[names(dataset@assays) != "RNA"]
if (length(assays.remove) > 0) { for (assay in assays.remove) { dataset[[assay]] <- NULL } }

## Output for scGPT (for review)
prepareGPT(dataset = dataset, batch_var = "batch_label", load = FALSE, outpath = "data/Mixture/mix.rds")

## Run scGPT (see scGPT.ipynb)

## Define dimensions to use for each method
dim.list <- list(scVI = 10, Harmony = 20, JOINTLY = 15, RPCA = 30, FastMNN = 20, LIGER = 20, Scanorama = 50, Unintegrated = 20)

## Embed
outfile <- "data/Mixture/mix_embeddings.rds"
embed <- embedMethods(data = dataset, batch_var = "batch_label", outpath = outfile)

## Add scGPT to embedding list
for (rep in 0:4) {
  # Import embeddings
  embedding <- read.delim(paste("data/Mixture/mix_scGPT_", rep, ".txt", sep=""), header=TRUE, sep = ",")
  rownames(embedding) <- embedding[,1]
  embedding <- embedding[,-1]
  colnames(embedding) <- paste("scGPT_", 1:ncol(embedding), sep="")
  embedding <- embedding[ match(colnames(dataset), rownames(embedding)),]
    
  # Insert into embed object
  embed$embeddings[[length(embed$embeddings)+1]] <- embedding
  names(embed$embeddings)[length(embed$embeddings)] <- paste("scGPT_rep", rep+1, sep="")
  print(rep)
}

## Save the embedding list
saveRDS(embed, "data/Mixture/mix_scGPT_embeddings.rds")

## Add labels
pancreas_labels <- read.delim("data/Pancreas/Transferred_labels.tsv", header = TRUE)
lung_labels <- read.delim("data/Lung/Transferred_labels.tsv", header = TRUE)
labels <- rbind(lung_labels, pancreas_labels)
labels <- labels[ match(colnames(dataset), labels$X),]
dataset$transfer_label <- labels[,2]
dataset <- subset(dataset, transfer_label != "")
dataset <- subset(dataset, transfer_label %in% names(which(table(dataset$transfer_label) >= 10)))

## Calculate best ARI for each embedding
metadata <- dataset@meta.data
summary <- as.data.frame(matrix(ncol=2, nrow=40))
for (embed.idx in 1:length(embed$embeddings)) {
  embed.mat <- as.matrix(embed$embeddings[[embed.idx]])
  embed.mat <- embed.mat[ rownames(embed.mat) %in% colnames(dataset),]
  embed.mat <- embed.mat[ match(colnames(dataset), rownames(embed.mat)),]
  method <- substr(names(embed$embeddings)[embed.idx],0,regexpr("_",names(embed$embeddings)[embed.idx])-1)
  dims <- dim.list[[which(names(dim.list) == method)]]

  # Setup to capture clustering metrics
  label_var <- "transfer_label"
  global_aris <- c()
  cl.min = 1; cl.max = 50
  cluster_list <- list()
  
  # Evaluate clustering metrics
  dend <- HGC.dendrogram(G = HGC::SNN.Construction(embed.mat[,1:dims]))
  for (cl in seq(cl.min, cl.max, 1)) {
    # Cluster
    if (cl == 1) {
      clusters <- rep(1, nrow(metadata))
      names(clusters) <- rownames(metadata)
    } else {
      clusters <- cutree(dend, k = cl)
      names(clusters) <- rownames(metadata)
    }
    cluster_list[[length(cluster_list)+1]] <- clusters
    
    # Capture global metrics for labels
    global_aris <- c(global_aris, aricode::ARI(clusters, factor(metadata[,label_var])))
  }

  # Save
  summary[embed.idx,1] <- embed.idx
  summary[embed.idx,2] <- max(global_aris)
}

## For each method, select the run with the best global ARI and create a UMAP
summary[,3] <- names(embed$embeddings)
summary[,4] <- substr(summary[,3], 0, regexpr("_", summary[,3])-1)
for (method in unique(summary[,4])) {
  embed.mat <- embed$embeddings[[summary[ summary[,4] %in% method,][which.max(summary[ summary[,4] %in% method,2]),1]]]
  embed.mat <- embed.mat[ rownames(embed.mat) %in% colnames(dataset),]
  embed.mat <- embed.mat[ match(colnames(dataset), rownames(embed.mat)),]
  dims <- dim.list[[which(names(dim.list) == method)]]
  dataset[[paste(method, "_raw",sep="")]] <- suppressWarnings(CreateDimReducObject(as.matrix(embed.mat)))
  dataset <- RunUMAP(dataset, dims = 1:dims, verbose = FALSE, reduction = paste(method, "_raw",sep=""), reduction.name = paste(method, "_umap",sep=""))
}

# Save object
saveRDS(dataset, "data/Mixture/mix_seurat_umaps.rds")

## Add scGPT to the results
dataset <- readRDS("data/Mixture/mix_seurat_umaps.rds")
aris <- c()
for (embed.idx in 41:length(embed$embeddings)) {
  embed.mat <- as.matrix(embed$embeddings[[embed.idx]])
  embed.mat <- embed.mat[ rownames(embed.mat) %in% colnames(dataset),]
  embed.mat <- embed.mat[ match(colnames(dataset), rownames(embed.mat)),]
  method <- substr(names(embed$embeddings)[embed.idx],0,regexpr("_",names(embed$embeddings)[embed.idx])-1)
  dims <- dim.list[[which(names(dim.list) == method)]]
  
  # Setup to capture clustering metrics
  label_var <- "transfer_label"
  global_aris <- c()
  cl.min = 1; cl.max = 50
  cluster_list <- list()
  
  # Evaluate clustering metrics
  dend <- HGC.dendrogram(G = HGC::SNN.Construction(embed.mat[,1:dims]))
  for (cl in seq(cl.min, cl.max, 1)) {
    # Cluster
    if (cl == 1) {
      clusters <- rep(1, nrow(metadata))
      names(clusters) <- rownames(metadata)
    } else {
      clusters <- cutree(dend, k = cl)
      names(clusters) <- rownames(metadata)
    }
    cluster_list[[length(cluster_list)+1]] <- clusters
    
    # Capture global metrics for labels
    global_aris <- c(global_aris, aricode::ARI(clusters, factor(metadata[,label_var])))
  }
  
  # Record best ARIs
  aris <- c(aris, max(global_aris))
}

## For each method, select the run with the best global ARI and create a UMAP
embed.mat <- embed$embeddings[[(which.max(aris)-1)+41]]
embed.mat <- embed.mat[ rownames(embed.mat) %in% colnames(dataset),]
embed.mat <- embed.mat[ match(colnames(dataset), rownames(embed.mat)),]
dims <- dim.list[[which(names(dim.list) == method)]]
dataset[[paste(method, "_raw",sep="")]] <- suppressWarnings(CreateDimReducObject(as.matrix(embed.mat)))
dataset <- RunUMAP(dataset, dims = 1:dims, verbose = FALSE, reduction = paste(method, "_raw",sep=""), reduction.name = paste(method, "_umap",sep=""))

# Save object
saveRDS(dataset, "data/Mixture/mix_seurat_scGPT_umaps.rds")

## Generate UMAPs of cell types with consistent colors
ds <- dataset
ds$new_label <- paste(ds$tissue, ds$transfer_label)
color_map = c(
  'Pancreas alpha'='#F27071',
  'Pancreas beta'='#F58771',
  'Pancreas delta'='#EF8F37',
  'Pancreas gamma'='#A55A26',
  'Pancreas epsilon'='#EC3132',
  'Pancreas acinar'='#55A445',
  'Pancreas ductal'='#59B847',
  'Pancreas activated_stellate'='#94BF3D',
  'Pancreas quiescent_stellate'='#5A7635',
  'Pancreas schwann'='#9CBC42',
  'Pancreas endothelial'='#35B34A',
  'Pancreas macrophage'='#0C8843',
  'Lung respiratory goblet cell'='#4D429A',
  'Lung lung ciliated cell'='#634EA0',
  'Lung adventitial cell'='#7E64AB',
  'Lung basal cell'='#9B8CB6',
  'Lung type I pneumocyte'='#C978B1',
  'Lung CD8-positive, alpha-beta T cell'='#BB6CAC',
  'Lung B cell'='#C273AD',
  'Lung plasma cell'='#AE7EB6',
  'Lung macrophage'='#954D9E',
  'Lung classical monocyte'='#AB6AAB',
  'Lung dendritic cell'='#C597C5',
  'Lung fibroblast'='#1FB1E7',
  'Lung blood vessel endothelial cell'='#065C9A',
  'Lung lung microvascular endothelial cell'='#3480C3',
  'Lung vein endothelial cell'='#7EB2E0'
)

# Generate the plots
methods <- c("scVI","Scanorama","Unintegrated", "FastMNN", "LIGER", "RPCA", "Harmony", "JOINTLY", "scGPT")
plt.list <- list()
for (method in methods) {
  if (method != methods[length(methods)]) {
    plt.list[[length(plt.list) + 1]] <- DimPlot(ds, reduction = paste(method, "_umap",sep=""), group.by = "new_label", cols = color_map) + ggtitle(method) & NoLegend()
  } else {
    plt.list[[length(plt.list) + 1]] <- DimPlot(ds, reduction = paste(method, "_umap",sep=""), group.by = "new_label", cols = color_map) + ggtitle(method)
  }
}

# Plot it
plt.list[[1]] + plt.list[[2]] + plt.list[[3]] + plt.list[[4]] + plt.list[[5]] + plt.list[[6]] + plt.list[[7]] + plt.list[[8]] + plt.list[[9]]

## Generate UMAPs of donors with consistent colors for each batch
color_map = c(
  'Donor_1'='#F27071',
  'Donor_2'='#F58771',
  'Donor_3'='#EF8F37',
  'Donor_4'='#A55A26',
  'Donor_5'='#EC3132',
  'Donor_6'='#55A445',
  'Donor_7'='#59B847',
  'Donor_8'='#94BF3D',
  'Donor_9'='#5A7635',
  'Donor_10'='#9CBC42',
  'Donor_11'='#35B34A',
  'Donor_12'='#0C8843',
  'A1'='#4D429A',
  'A2'='#634EA0',
  'A3'='#7E64AB',
  'A4'='#9B8CB6',
  'A5'='#C978B1',
  'A6'='#BB6CAC'
)

# Generate the plots
methods <- c("scVI","Scanorama","Unintegrated", "FastMNN", "LIGER", "RPCA", "Harmony", "JOINTLY", "scGPT")
plt.list <- list()
for (method in methods) {
  if (method != methods[length(methods)]) {
    plt.list[[length(plt.list) + 1]] <- DimPlot(ds, reduction = paste(method, "_umap",sep=""), group.by = "batch_label", cols = color_map) + ggtitle(method) & NoLegend()
  } else {
    plt.list[[length(plt.list) + 1]] <- DimPlot(ds, reduction = paste(method, "_umap",sep=""), group.by = "batch_label", cols = color_map) + ggtitle(method)
  }
}

# Plot it
plt.list[[1]] + plt.list[[2]] + plt.list[[3]] + plt.list[[4]] + plt.list[[5]] + plt.list[[6]] + plt.list[[7]] + plt.list[[8]] + plt.list[[9]]

## Differential expression within mixed cell types
# Set labels
ds@meta.data[ grep("endothelial cell", ds@meta.data$transfer_label),"transfer_label"] <- "endothelial"
ds$ct <- paste(ds$transfer_label, ds$tissue, sep="_")
ds$ps <- paste(ds$ct, ds$batch_label, sep="_")

# Calculate pseudo-bulk expression levels
psb <- matrix(ncol=length(unique(ds$ps)), nrow = nrow(ds))
colnames(psb) <- unique(ds$ps)
rownames(psb) <- rownames(ds@assays$RNA@counts)
for (ps in unique(ds$ps)) {
  psb[,colnames(psb) == ps] <- rowSums(ds@assays$RNA@counts[, colnames(ds@assays$RNA@counts) %in% rownames(ds@meta.data[ ds@meta.data$ps == ps,]), drop = FALSE])
}
psb <- psb[ which(rowSums(psb) >= 10), ]
md <- data.frame(sample = colnames(psb))
md$label <- "NA"
md$donor <- "NA"
md$celltype <- "NA"
md$tissue <- "NA"
for (dnr in unique(ds$batch_label)) {  md[ grep(dnr, md$sample), "donor"] <- dnr }
for (ct in unique(ds$ct)) {  md[ grep(ct, md$sample), "label"] <- ct }
for (ct in unique(ds$transfer_label)) {  md[ grep(ct, md$label), "celltype"] <- ct }
for (ct in unique(ds$tissue)) {  md[ grep(ct, md$label), "tissue"] <- ct }

# Setup DE tests
DGE <- DGEList(counts = psb, group = md$label)
DGE <- calcNormFactors(DGE)
design <- model.matrix(~ md$label + 0)
colnames(design) <- sort(unique(md$label))
DGE <- estimateDisp(DGE, design)
fit <- glmQLFit(DGE, design)

# Perform all pairwise tests involving endothelial cells
tests <- expand.grid(1:25,1:25)
tests <- tests[ tests[,1] != tests[,2],]
pvals <- fdrs <- logFCs <- as.data.frame(matrix(100, nrow=nrow(psb), ncol=nrow(tests)))
counter <- 1
for (test in 1:nrow(tests)) {
  positive <- tests[test,1]
  negative <- tests[test,2]
  name <- colnames(design)[positive]
  if (grepl("ndothelial", name)) {
    coef <- rep(0, 25)
    coef[positive] <- 1
    coef[negative] <- -1
    markers <- glmQLFTest(fit, contrast = coef)
    markers <- as.data.frame(topTags(markers, n = nrow(psb), sort.by = "none"))
    pvals[,counter] <- markers[,4]
    fdrs[,counter] <- markers[,5]
    logFCs[,counter] <- markers[,1]
    colnames(pvals)[counter] <- colnames(fdrs)[counter]  <- colnames(logFCs)[counter] <- paste0(colnames(design)[c(positive, negative)], collapse="-")
    rownames(pvals) <- rownames(fdrs) <- rownames(logFCs) <- rownames(markers)
    counter <- counter + 1
  }
}

# Select endothelial cell markers (removing genes with values 1364600)
pvals <- pvals[ ,colSums(pvals) != 1364600]
fdrs <- fdrs[ ,colSums(fdrs) != 1364600]
logFCs <- logFCs[ ,colSums(logFCs) != 1364600]
states <- c(grep("endothelial_Lung-endothelial_Pancreas", colnames(pvals)),grep("endothelial_Pancreas-endothelial_Lung", colnames(pvals)))
lung.idx <- grep("Lung", colnames(pvals))
lung.idx <- lung.idx[!(lung.idx %in% states)]
pancreas.idx <- grep("Pancreas", colnames(pvals))
pancreas.idx <- pancreas.idx[!(pancreas.idx %in% states)]
HMP.lung <- apply(pvals[,lung.idx],1,FUN= function(x) { 1/mean(1/x) })
HMP.lung.FDR <- p.adjust(HMP.lung, "fdr")
HMP.pancreas <- apply(pvals[,pancreas.idx],1,FUN= function(x) { 1/mean(1/x) })
HMP.pancreas.FDR <- p.adjust(HMP.pancreas, "fdr")
genes.foldchange.include <- c()
for (i in 1:nrow(logFCs)) {
  if (min(logFCs[i,-states]) >= log2(1.5)) {
    genes.foldchange.include <- c(genes.foldchange.include, rownames(logFCs)[i])
  }
}
genes.include <- intersect(genes.foldchange.include, names(which(HMP.lung.FDR <= 0.01)))
genes.include <- intersect(genes.include, names(which(HMP.pancreas.FDR <= 0.01)))

# Generate a heatmap of endothelial cell markers
cnts <- cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 5)
ord <- order(md$tissue, md$celltype)
md <- md[ ord,]
cnts <- cnts[,match(md$sample, colnames(cnts))]
annot <- HeatmapAnnotation(Tissue = md$tissue, Celltype = md$celltype)
Heatmap(t(scale(t(cnts[ rownames(cnts) %in% genes.foldchange.include,]))), bottom_annotation = annot, show_column_dend = FALSE, show_row_dend = FALSE, show_column_names = FALSE, cluster_columns = FALSE)

# Select endothelial state markers
genes.fdr.include <- names(which(rowSums(fdrs[, states] <= 0.01) > 0))
genes.logFC.include <- names(which(apply(logFCs[, states],1,FUN="max") >= log2(1.5)))
genes.include <- intersect(genes.logFC.include, genes.fdr.include)

# Generate a heatmap of endothelial state markers
cnts.plot <- t(scale(t(cnts[rownames(cnts) %in% genes.include,grep("ndothelial", colnames(cnts))])))
annot <- HeatmapAnnotation(tissue = c(rep("Lung",6), rep("Pancreas",12)))
set.seed(42)
htm <- Heatmap(cnts.plot, bottom_annotation = annot, show_row_dend = FALSE, show_row_names = FALSE, show_column_names = FALSE, km = 2)
htm <- draw(htm)

### Pathway analysis
## Define gene sets
pancreas.genes <- rownames(cnts.plot)[row_order(htm)[[1]]]
lung.genes <- rownames(cnts.plot)[row_order(htm)[[2]]]

## Enriched TFs (with enrichR)
pancreas.tfs <- enrichR::enrichr(pancreas.genes, "ChEA_2022")[[1]]
pancreas.tfs$tissue <- "Pancreas"
pancreas.tfs$TF <- substr(pancreas.tfs$Term, 0, regexpr(" ", pancreas.tfs$Term)-1)
pancreas.tfs <- pancreas.tfs[ order(pancreas.tfs$TF, pancreas.tfs$Adjusted.P.value),]
pancreas.tfs <- pancreas.tfs[ duplicated(pancreas.tfs$TF)==FALSE,]
pancreas.tfs <- pancreas.tfs[ order(pancreas.tfs$Adjusted.P.value),]
pancreas.tfs <- pancreas.tfs[,c("TF","tissue","Adjusted.P.value")]
lung.tfs <- enrichR::enrichr(lung.genes, "ChEA_2022")[[1]]
lung.tfs$tissue <- "Lung"
lung.tfs$TF <- substr(lung.tfs$Term, 0, regexpr(" ", lung.tfs$Term)-1)
lung.tfs <- lung.tfs[ order(lung.tfs$TF, lung.tfs$Adjusted.P.value),]
lung.tfs <- lung.tfs[ duplicated(lung.tfs$TF)==FALSE,]
lung.tfs <- lung.tfs[ order(lung.tfs$Adjusted.P.value),]
lung.tfs <- lung.tfs[,c("TF","tissue","Adjusted.P.value")]
tfs.include <- c(lung.tfs[1:10,"TF"],pancreas.tfs[1:10,"TF"])
tfs <- rbind(lung.tfs[lung.tfs$TF %in% tfs.include,], pancreas.tfs[pancreas.tfs$TF %in% tfs.include,])
colnames(tfs)[3] <- "padj"
tfs[,3] <- -log10(tfs[,3])
tfs <- tfs[ order(tfs$TF),]
tfs$ordering <- 0
for (i in seq(2, nrow(tfs),2)) { tfs[i,4] <- tfs[(i-1),4] <- tfs[i,3] - tfs[(i-1),3]}
ggplot(data = tfs, aes(y = reorder(TF, -ordering), x = padj, color = tissue)) + geom_point() + theme_minimal()

## IFNy genes
# Extract IFNy genes from all genes and Reactome
IFNy <- enrichR::enrichr(rownames(ds), "Reactome_2022")[[1]]
IFNy <- strsplit(IFNy[ IFNy$Term == "Interferon Gamma Signaling R-HSA-877300","Genes"], ";")[[1]]
IFNy.subset <- IFNy[IFNy %in% genes.include]

# Generate a heatmap of endothelial state markers
cnts.plot <- t(scale(t(cnts[rownames(cnts) %in% IFNy.subset,grep("ndothelial", colnames(cnts))])))
annot <- HeatmapAnnotation(tissue = c(rep("Lung",6), rep("Pancreas",12)))
htm <- Heatmap(cnts.plot, bottom_annotation = annot, show_row_dend = FALSE, show_row_names = TRUE, show_column_names = FALSE, km = 2)
htm <- draw(htm)

### Over-correction evaluation with LISI
## Calculate iLISI and cLISI within and across tissues
lisi.res <- as.data.frame(matrix(ncol=6,nrow=1))
lisi.counter <- 1
dataset$new_label <- paste(dataset$tissue, dataset$transfer_label)
md <- dataset@meta.data
for (method in names(dim.list)[names(dim.list) != "Unintegrated"]) {
  lisi.res[lisi.counter,1] <- method
  overall.lisi <- compute_lisi(dataset[[paste(method, "_raw", sep="")]]@cell.embeddings[,1:dim.list[[which(names(dim.list) == method)]]], meta_data = md, label_colnames = c("tissue", "batch_label"))
  lisi.res[lisi.counter,2] <- mean((overall.lisi[,1]-1) / (2-1))
  lung.lisi <- compute_lisi(dataset[[paste(method, "_raw", sep="")]]@cell.embeddings[which(dataset$tissue == "Lung"),1:dim.list[[which(names(dim.list) == method)]]], meta_data = md[which(dataset$tissue == "Lung"),], label_colnames = c("new_label","batch_label"))
  lisi.res[lisi.counter,3] <- mean((lung.lisi[,1]-1) / (length(unique(md[which(dataset$tissue == "Lung"),"batch_label"]))-1))
  lisi.res[lisi.counter,4] <- mean((length(unique(md[which(dataset$tissue == "Lung"),"new_label"]))-lung.lisi[,2]) / (length(unique(md[which(dataset$tissue == "Lung"),"new_label"]))-1))
  pancreas.lisi <- compute_lisi(dataset[[paste(method, "_raw", sep="")]]@cell.embeddings[which(dataset$tissue == "Pancreas"),1:dim.list[[which(names(dim.list) == method)]]], meta_data = md[which(dataset$tissue == "Pancreas"),], label_colnames = c("new_label","batch_label"))
  lisi.res[lisi.counter,5] <- mean((pancreas.lisi[,1]-1) / (length(unique(md[which(dataset$tissue == "Pancreas"),"batch_label"]))-1))
  lisi.res[lisi.counter,6] <- mean((length(unique(md[which(dataset$tissue == "Pancreas"),"new_label"]))-pancreas.lisi[,2]) / (length(unique(md[which(dataset$tissue == "Pancreas"),"new_label"]))-1))
  lisi.counter <- lisi.counter + 1
}

## Extract ranks for iLISI across and within tissues, and cLISI within tissues
ranks <- rbind(
  rank(-signif(lisi.res[,3],3), ties.method = 'min'),
  rank(-signif(lisi.res[,5],3), ties.method = 'min'),
  rank(-signif(lisi.res[,4],3), ties.method = 'min'),
  rank(-signif(lisi.res[,6],3),ties.method = 'min'),
  rank(signif(lisi.res[,2],3),ties.method = 'min'))

## Calculate an overall rank
ranks <- rbind(ranks, rank(colMeans(rbind(rank(colMeans(ranks[c(1,2),]),ties.method = 'min'),rank(colMeans(ranks[c(3,4),]),ties.method = 'min'), ranks[5,])), ties.method = 'min'))

## Setup the data
ranks <- as.data.frame(ranks)
colnames(ranks) <- names(dim.list)[names(dim.list) != "Unintegrated"]
rownames(ranks) <- c("iLISI_Lung","iLISI_Pancreas","cLISI_Lung","cLISI_Pancreas","Over-correction", "Overall")

## Plot the heatmap
Heatmap(ranks ,
        show_heatmap_legend = F,
        heatmap_legend_param = list(color_bar = "discrete"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "black", lwd = 1),
        col = colorRamp2(c(1, 7), c("#EFF3FF", "#2171B5")),
        row_names_side = "left")

## Write the metrics to a spreadsheet
colnames(lisi.res) <- c("Method","Overcorrection","iLISI_Lung","iLISI_Pancreas","cLISI_Lung","cLISI_Pancreas")
writexl::write_xlsx(lisi.res, path = "overcorrection_supplementary_table.xlsx")
