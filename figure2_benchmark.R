### Setup
## Options
options(future.globals.maxSize= 20000*1024^2)

## Dependencies
source("dependencies.R")

## Libraries
library(Seurat)
library(SeuratDisk)
library(BiocParallel)
library(patchwork)
library(data.table)
library(HGC)
library(ggplot2)
library(JOINTLY)
library(circlize)
library(lisi)
library(cowplot)
library(writexl)
library(reticulate)
use_virtualenv("/work/PythonEnv/scVI_new/")

## Functions
source("functions.R")

## Define dimensions to use for each method
dim.list <- list(scVI = 10, Harmony = 20, JOINTLY = 20, RPCA = 30, FastMNN = 20, LIGER = 20, Scanorama = 50, Unintegrated = 20, scGPT = 512)

## Define file list
files <- c(
  "data/Liver/Human_Liver.rds",
  "data/Pancreas/Human_Pancreas.rds",
  "data/Kidney/Human_Kidney_sub.rds",
  "data/Lung/Human_Lung_sub.rds",
  "data/PBMC/PBMC.rds")

## Loop across files, embed and evaluate
for (file in files) {
  dataset <- readRDS(file)
  dataset$batch_label <- as.character(dataset$batch_label)
  assays.remove <- names(dataset@assays)[names(dataset@assays) != "RNA"]
  if (length(assays.remove) > 0) { for (assay in assays.remove) { dataset[[assay]] <- NULL } }
  outfile <- gsub(".rds", "_embeddings.rds", file)
  embed <- embedMethods(data = dataset, batch_var = "batch_label", outpath = outfile)
  for (i in 1:5) { colnames(embed$embeddings[[i]]) <- paste("JOINTLY", 1:ncol(embed$embeddings[[i]]), sep="_") }
  labels <- read.delim(paste(substr(file, 0, nchar(file)- regexpr("/", intToUtf8(rev(utf8ToInt(file))))+1), "Transferred_labels.tsv", sep=""), header=TRUE)
  labels <- labels[ match(colnames(dataset), labels$X),]
  dataset$transfer_label <- labels[,2]
  dataset <- subset(dataset, transfer_label != "")
  dataset <- subset(dataset, transfer_label %in% names(which(table(dataset$transfer_label) >= 10)))
  metadata <- dataset@meta.data
  outfile <- gsub(".rds", "_result.rds", file)
  result <- bplapply(X = embed$embeddings, FUN = evaluateEmbedding, BPPARAM = SerialParam(), metadata = metadata, batch_var = "batch_label", label_var = "transfer_label", dim.list = dim.list, cl.min = 1, cl.max = 50)
  saveRDS(result, outfile)
  outfile <- gsub(".rds", "_summary.rds", file)
  summary <- summarizeResults(result = result, data = dataset, embeddings = embed$embeddings, dim.list = dim.list, selectRun = "globalARI")
  saveRDS(summary, outfile)
  rm(dataset)
  rm(embed)
  gc()
}

## Add scGPT to evaluation (review question)
# Export adata objects
for (file in files) { prepareGPT(file, "batch_label") }

# Run fine-tuning (see scGPT.ipynb)

# Add results
for (file in files) {
  # Import the embeddings
  outfile <- gsub(".rds", "_embeddings.rds", file)
  newembed <- gsub(".rds", "_scGPT_embeddings.rds", file)
  embed <- readRDS(outfile)
  
  # Import the results from other methods
  outfile <- gsub(".rds", "_result.rds", file)
  newresult <- gsub(".rds", "_scGPT_result.rds", file)
  result <- readRDS(outfile)
  
  # Define new output path for summary
  outfile <- gsub(".rds", "_summary.rds", file)
  newsummary <- gsub(".rds", "_scGPT_summary.rds", file)
  summary <- readRDS(outfile)
  
  # Import and process the dataset
  dataset <- readRDS(file)
  dataset$batch_label <- as.character(dataset$batch_label)
  assays.remove <- names(dataset@assays)[names(dataset@assays) != "RNA"]
  if (length(assays.remove) > 0) { for (assay in assays.remove) { dataset[[assay]] <- NULL } }
  labels <- read.delim(paste(substr(file, 0, nchar(file)- regexpr("/", intToUtf8(rev(utf8ToInt(file))))+1), "Transferred_labels.tsv", sep=""), header=TRUE)
  labels <- labels[ match(colnames(dataset), labels$X),]
  dataset$transfer_label <- labels[,2]
  dataset <- subset(dataset, transfer_label != "")
  dataset <- subset(dataset, transfer_label %in% names(which(table(dataset$transfer_label) >= 10)))
  metadata <- dataset@meta.data

  # Get the correct path
  file <- gsub("data/", "", gsub(".rds", "", gsub("_sub", "", file)))
  flder <- substr(file, 0, regexpr("/", file)-1)
  file <- substr(file, regexpr("/", file)+1, nchar(file))
  
  # Import and evaluate the embeddings from scGPT
  for (rep in 0:4) {
    # Import embeddings
    embedding <- read.delim(paste("data/", flder, "/", file, "_scGPT_", rep, ".txt", sep=""), header=TRUE, sep = ",")
    rownames(embedding) <- embedding[,1]
    embedding <- embedding[,-1]
    colnames(embedding) <- paste("scGPT_", 1:ncol(embedding), sep="")
    embedding <- embedding[ match(colnames(dataset), rownames(embedding)),]
    
    # Insert into embed object
    embed$embeddings[[length(embed$embeddings)+1]] <- embedding
    names(embed$embeddings)[length(embed$embeddings)] <- paste("scGPT_rep", rep+1, sep="")
    
    # Evaluate the embedding
    eval.res <- evaluateEmbedding(embedding, metadata = metadata, batch_var = "batch_label", label_var = "transfer_label", dim.list = dim.list, cl.min = 1, cl.max = 50)
    result[[length(result)+1]] <- eval.res
    names(result)[length(result)] <- paste("scGPT_rep", rep+1, sep="")
  }
  
  # Add to summary
  data <- summary$object
  summary.tmp <- do.call("rbind", lapply(result, FUN = function(x) { x$summary}))
  summary.tmp$Method <- substr(rownames(summary.tmp),0,regexpr("_", rownames(summary.tmp))-1)
  idx <- which(summary.tmp$Method == "scGPT" & summary.tmp[,"globalARI"] == max(summary.tmp[ summary.tmp$Method == "scGPT", "globalARI"]))[1]
  embed.tmp <- embed$embeddings[[idx]]
  embed.tmp <- embed.tmp[ rownames(embed.tmp) %in% colnames(data),]
  embed.tmp <- embed.tmp[ match(colnames(data), rownames(embed.tmp)),]
  data[["scGPT"]] <- Seurat::CreateDimReducObject(as.matrix(embed.tmp), assay = "RNA")
  clusters <- result[[idx]]$clusters[[which.max(result[[idx]]$cluster_metrics[,"globalARI"])]]
  clusters <- clusters[ names(clusters) %in% colnames(data)]
  clusters <- clusters[ match(colnames(data), names(clusters))]
  data@meta.data[,(ncol(data@meta.data)+1)] <- clusters
  colnames(data@meta.data)[ncol(data@meta.data)] <- paste("scGPT", "_clusters", sep="")
  dims <- dim.list[[which(names(dim.list) == "scGPT")]]
  data <- RunUMAP(data, dims = 1:dims, reduction = "scGPT", reduction.name = paste("scGPT", "_UMAP", sep=""), verbose = FALSE)
  summary$object <- data
  summary$summary <- rbind(summary$summary, summary.tmp[idx,])

  # Save all objects
  saveRDS(embed, newembed)
  saveRDS(result, newresult)
  saveRDS(summary, newsummary)
}
  
##### HEATMAPS
# Import results
lung <- readRDS("data/Lung/Human_Lung_sub_scGPT_summary.rds")
liver <- readRDS("data/Liver/Human_Liver_scGPT_summary.rds")
pancreas <- readRDS("data/Pancreas/Human_Pancreas_scGPT_summary.rds")
kidney <- readRDS("data/Kidney/Human_Kidney_sub_scGPT_summary.rds")
pbmc <- readRDS("data/PBMC/PBMC_scGPT_summary.rds")

# Setup 
method_order = c('JOINTLY' ,'scVI', 'Harmony', 'FastMNN', 'Scanorama', 'RPCA', 'LIGER', 'scGPT', 'Unintegrated')
data_order = c('All global', 'pancreas global', 'lung global', 'liver global', 'pbmc global', 'kidney global',
               'All worst', 'pancreas worst', 'lung worst', 'liver worst', 'pbmc worst', 'kidney worst')
metric = colnames(kidney$summary)
metric_select_global = c(1,3,5, 11, 13, 15,17)
metric_select_worst = c(2,4,6, 12, 14, 16,18)
metric_select_global_objects = list()
metric_select_worst_objects = list()

## Run
# Global
for (i in 1:length(metric_select_global)){
  metric_select_global_objects[[i]] = process_summary(metric_select_global[i], label = 'global')
}

# Worst
for (i in 1:length(metric_select_worst)){
  metric_select_worst_objects[[i]] = process_summary(metric_select_worst[i], label = 'worst')
}

# Create heatmaps
ht_list = NULL
metric_sum_g = 0
metric_sum_w = 0
n_metrics = length(metric_select_global)
for (i in 1:n_metrics){
  title = metric[metric_select_global[i]]
  title = substr(title, start=7, stop = nchar(title))
  toplot = rbind(metric_select_global_objects[[i]]$rank,
                 metric_select_worst_objects[[i]]$rank)[data_order,method_order]
  ht_list = ht_list + Heatmap(toplot ,
                              row_split = c(rep(c("Global \nperformance"), 6),  rep(c("Worst \nperformance"), 6)),
                              show_heatmap_legend = F,
                              heatmap_legend_param = list(color_bar = "discrete"),
                              column_title =title,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              rect_gp = gpar(col = "black", lwd = 1),
                              col = colorRamp2(c(1, 9), c("#EFF3FF", "#2171B5")),
                              row_names_side = "left",
                              
  )
}

# Plot it
ht_list

##### UMAPs
## Kidney
best_kidney = rownames(kidney$summary)
DimPlot(kidney$object, group.by = c('transfer_label', 'JOINTLY_clusters'), reduction = 'JOINTLY_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(kidney$object, group.by = c('transfer_label', 'Harmony_clusters'), reduction = 'Harmony_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(kidney$object, group.by = c('transfer_label', 'RPCA_clusters'), reduction = 'RPCA_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(kidney$object, group.by = c('transfer_label', 'LIGER_clusters'), reduction = 'LIGER_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(kidney$object, group.by = c('transfer_label', 'FastMNN_clusters'), reduction = 'FastMNN_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(kidney$object, group.by = c('transfer_label', 'scVI_clusters'), reduction = 'scVI_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(kidney$object, group.by = c('transfer_label', 'Scanorama_clusters'), reduction = 'Scanorama_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(kidney$object, group.by = c('transfer_label', 'Unintegrated_clusters'), reduction = 'Unintegrated_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(kidney$object, group.by = c('transfer_label', 'scGPT_clusters'), reduction = 'scGPT_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()

## Liver
best_liver = rownames(liver$summary)
DimPlot(liver$object, group.by = c('transfer_label', 'JOINTLY_clusters'), reduction = 'JOINTLY_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(liver$object, group.by = c('transfer_label', 'Harmony_clusters'), reduction = 'Harmony_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(liver$object, group.by = c('transfer_label', 'RPCA_clusters'), reduction = 'RPCA_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(liver$object, group.by = c('transfer_label', 'LIGER_clusters'), reduction = 'LIGER_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(liver$object, group.by = c('transfer_label', 'FastMNN_clusters'), reduction = 'FastMNN_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(liver$object, group.by = c('transfer_label', 'scVI_clusters'), reduction = 'scVI_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(liver$object, group.by = c('transfer_label', 'Scanorama_clusters'), reduction = 'Scanorama_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(liver$object, group.by = c('transfer_label', 'Unintegrated_clusters'), reduction = 'Unintegrated_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(liver$object, group.by = c('transfer_label', 'scGPT_clusters'), reduction = 'scGPT_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()

## Pancreas
best_pancreas = rownames(pancreas$summary)
DimPlot(pancreas$object, group.by = c('transfer_label', 'JOINTLY_clusters'), reduction = 'JOINTLY_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pancreas$object, group.by = c('transfer_label', 'Harmony_clusters'), reduction = 'Harmony_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pancreas$object, group.by = c('transfer_label', 'RPCA_clusters'), reduction = 'RPCA_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pancreas$object, group.by = c('transfer_label', 'LIGER_clusters'), reduction = 'LIGER_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pancreas$object, group.by = c('transfer_label', 'FastMNN_clusters'), reduction = 'FastMNN_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pancreas$object, group.by = c('transfer_label', 'scVI_clusters'), reduction = 'scVI_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pancreas$object, group.by = c('transfer_label', 'Scanorama_clusters'), reduction = 'Scanorama_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pancreas$object, group.by = c('transfer_label', 'Unintegrated_clusters'), reduction = 'Unintegrated_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pancreas$object, group.by = c('transfer_label', 'scGPT_clusters'), reduction = 'scGPT_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()

## Lung
best_lung = rownames(lung$summary)
lung_colors = c('blood vessel endothelial cell' = '#FF6166',
                'vein endothelial cell' = '#FF3339',
                'lung microvascular endothelial cell' = '#FF050D',
                
                'B cell' = '#4D96FF',
                'plasma cell' = '#267AF2',
                'CD8-positive, alpha-beta T cell' = '#005EE5',
                
                'classical monocyte' = '#00EB4B',
                'macrophage' = '#00CF42',
                'dendritic cell' = '#00B339',
                
                'adventitial cell' = '#E1B10E',
                'respiratory goblet cell' = '#DF01C1',
                'type I pneumocyte' = '#FF4DED',
                'basal cell' = '#8E00B3',
                'fibroblast' = '#0EE184',
                'lung ciliated cell' = '#00bad6',
                
                '1' = '#f8766d',
                '2' = '#b89000',
                '3' = '#53b300',
                '4' = '#00b88d',
                '5' = '#00b4eb',
                '6' = '#a58aff',
                '7' = '#fb60d7'
                
)
DimPlot(lung$object, group.by = c('transfer_label', 'JOINTLY_clusters'), cols = lung_colors, reduction = 'JOINTLY_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(lung$object, group.by = c('transfer_label', 'Harmony_clusters'), cols = lung_colors, reduction = 'Harmony_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(lung$object, group.by = c('transfer_label', 'RPCA_clusters'), cols = lung_colors, reduction = 'RPCA_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(lung$object, group.by = c('transfer_label', 'LIGER_clusters'), cols = lung_colors, reduction = 'LIGER_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(lung$object, group.by = c('transfer_label', 'FastMNN_clusters'), cols = lung_colors, reduction = 'FastMNN_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(lung$object, group.by = c('transfer_label', 'scVI_clusters'), cols = lung_colors, reduction = 'scVI_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(lung$object, group.by = c('transfer_label', 'Scanorama_clusters'), cols = lung_colors, reduction = 'Scanorama_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(lung$object, group.by = c('transfer_label', 'Unintegrated_clusters'), cols = lung_colors, reduction = 'Unintegrated_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(lung$object, group.by = c('transfer_label', 'scGPT_clusters'), reduction = 'scGPT_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()

## PBMC
best_pbmc = rownames(pbmc$summary)
DimPlot(pbmc$object, group.by = c('transfer_label', 'JOINTLY_clusters'), reduction = 'JOINTLY_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pbmc$object, group.by = c('transfer_label', 'Harmony_clusters'), reduction = 'Harmony_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pbmc$object, group.by = c('transfer_label', 'RPCA_clusters'), reduction = 'RPCA_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pbmc$object, group.by = c('transfer_label', 'LIGER_clusters'), reduction = 'LIGER_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pbmc$object, group.by = c('transfer_label', 'FastMNN_clusters'), reduction = 'FastMNN_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pbmc$object, group.by = c('transfer_label', 'scVI_clusters'), reduction = 'scVI_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pbmc$object, group.by = c('transfer_label', 'Scanorama_clusters'), reduction = 'Scanorama_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pbmc$object, group.by = c('transfer_label', 'Unintegrated_clusters'), reduction = 'Unintegrated_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()
DimPlot(pbmc$object, group.by = c('transfer_label', 'scGPT_clusters'), reduction = 'scGPT_UMAP', label = TRUE, ncol = 1)  &  theme_classic(base_size = 4) & guides(color = guide_legend(override.aes = list(size=3), ncol=1) ) & NoAxes()

##### Supplementary table
## Import data
lung <- readRDS("data/Lung/Human_Lung_sub_scGPT_result.rds")
liver <- readRDS("data/Liver/Human_Liver_scGPT_result.rds")
pancreas <- readRDS("data/Pancreas/Human_Pancreas_scGPT_result.rds")
kidney <- readRDS("data/Kidney/Human_Kidney_sub_scGPT_result.rds")
pbmc <- readRDS("data/PBMC/PBMC_scGPT_result.rds")
data.list <- list(lung, liver, pancreas, kidney, pbmc)
names(data.list) <- c("Lung","Liver","Pancreas","Kidney","PBMC")

## All results
results <- data.frame(matrix(nrow = 40, ncol = 17))
counter <- 1
for (m in 1:5) {
  ds <- names(data.list)[m]
  tmp <- data.list[[m]]
  for (i in 1:45) {
    results[counter,4:17] <- tmp[[i]]$summary[,c(1:6,11:18)]
    results[counter,1] <- ds
    results[counter,2] <- substr(names(tmp)[i],0,regexpr("_", names(tmp)[i])-1)
    results[counter,3] <- substr(names(tmp)[i],regexpr("_", names(tmp)[i])+4, nchar(names(tmp)[i]))
    counter <- counter + 1
  }
}
colnames(results) <- c("Dataset", "Method", "Replicate", colnames(tmp[[1]]$summary)[c(1:6,11:18)])
writexl::write_xlsx(results, path = "benchmarking_supplementary_table.xlsx")

#### Subsetting analysis on human liver
### Subset and integrate
## Load the dataset and clean it
dataset <- readRDS(files[1])
dataset$batch_label <- as.character(dataset$batch_label)
assays.remove <- names(dataset@assays)[names(dataset@assays) != "RNA"]
if (length(assays.remove) > 0) { for (assay in assays.remove) { dataset[[assay]] <- NULL } }

## Calculate how many batches to maximally remove
nbatches <- length(unique(dataset$batch_label))
max.remove <- nbatches - 2

## Setup to keep track of removed datasets
counter <- 1
tracking <- as.data.frame(matrix(ncol=4, nrow = 1))

## Loop across all possible combinations to remove
for (remove in 1:max.remove) {
  remove.ds <- as.data.frame(t(combn(unique(dataset$batch_label), remove)))
  for (rm.idx in 1:nrow(remove.ds)) {
    # Remove datasets
    ds.remove <- c()
    for (col in 1:ncol(remove.ds)) { ds.remove <- c(ds.remove, remove.ds[rm.idx,col])}
    dataset.subset <- subset(dataset, cells = rownames(dataset@meta.data[ !(dataset@meta.data$batch_label %in% ds.remove),]))
    if (!identical(character(0), names(dataset.subset@reductions))) {
      for (red.remove in names(dataset.subset@reductions)) {
        dataset.subset[[red.remove]] <- NULL
      }
    }
    
    # Run JOINTLY 5 times
    proc <- preprocess(dataset.subset, "batch_label", verbose = FALSE)
    cpca <- cpca(proc, verbose = FALSE)
    inputs <- prepareData(cpca$cpca, verbose = FALSE)
    for (rep in 1:5) {
      # Embed and evaluate
      solved <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca, bpparam = BiocParallel::MulticoreParam(), progressbar = FALSE, share.objects = FALSE, verbose = FALSE)
      H <- t(do.call("cbind",solved$Hmat))
      H <- H[ match(colnames(dataset.subset), rownames(H)),]
      H <- scale(H)
      H <- t(scale(t(H)))
      colnames(H) <- paste("JOINTLY_rep", rep,"_",1:ncol(H), sep="")
      dataset.subset[[paste("JOINTLY_rep",rep, sep="")]] <- CreateDimReducObject(H, assay = "RNA")
      rm(H)
    }
    
    # Save the object
    saveRDS(dataset.subset, paste("data/Subsets/object", counter, ".Rds", sep=""))
    
    # Tracking
    tracking[counter, 1] <- file
    tracking[counter, 2] <- paste("data/Subsets/object", counter, ".Rds", sep="")
    tracking[counter, 3] <- remove
    tracking[counter, 4] <- paste0(ds.remove, collapse = "|")
    counter <- counter + 1
  }
}

## Save the tracking object for future reference
saveRDS(tracking, "data/Subsets/tracking.Rds")

### Evaluate
## Full dataset
# Setup to capture results
results <- as.data.frame(matrix(ncol=7, nrow = 1))
results.counter <- 1

# Embed
file <- files[1]
dataset <- readRDS(file)
dataset$batch_label <- as.character(dataset$batch_label)
assays.remove <- names(dataset@assays)[names(dataset@assays) != "RNA"]
if (length(assays.remove) > 0) { for (assay in assays.remove) { dataset[[assay]] <- NULL } }

proc <- preprocess(dataset, "batch_label", verbose = FALSE)
cpca <- cpca(proc, verbose = FALSE)
inputs <- prepareData(cpca$cpca, verbose = FALSE)
for (rep in 1:6) {
  # Embed and evaluate
  solved <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca, bpparam = BiocParallel::MulticoreParam(), progressbar = TRUE, share.objects = FALSE, verbose = TRUE)
  H <- t(do.call("cbind",solved$Hmat))
  H <- H[ match(colnames(dataset), rownames(H)),]
  H <- scale(H)
  H <- t(scale(t(H)))
  colnames(H) <- paste("JOINTLY_rep", rep,"_",1:ncol(H), sep="")
  dataset[[paste("JOINTLY_rep",rep, sep="")]] <- CreateDimReducObject(H, assay = "RNA")
  rm(H)
}

# Get labels
labels <- read.delim(paste(substr(file, 0, nchar(file)- regexpr("/", intToUtf8(rev(utf8ToInt(file))))+1), "Transferred_labels.tsv", sep=""), header=TRUE)
labels <- labels[ match(colnames(dataset), labels$X),]
dataset$transfer_label <- labels[,2]
dataset <- subset(dataset, transfer_label != "")
dataset <- subset(dataset, transfer_label %in% names(which(table(dataset$transfer_label) >= 10)))

## LISI and ASW
global_cLISI <- global_iLISI <- global_bASW <- global_cASW <- global_ARI <- c()
md.subset <- dataset@meta.data

# LISI
for (i in 1:5) {
  LISIS <- compute_lisi(dataset[[paste("JOINTLY_rep",i, sep="")]]@cell.embeddings[,1:20], dataset@meta.data, label_colnames = c("batch_label", "transfer_label"))
  global_cLISI <- c(global_cLISI, median((length(unique(dataset@meta.data[,"transfer_label"]))-LISIS[,2]) / (length(unique(dataset@meta.data[,"transfer_label"])) - 1)))
  global_iLISI <- c(global_iLISI, median((LISIS[,1]-1) / (length(unique(dataset@meta.data[,"batch_label"])) - 1)))

  # Distance mat
  space <- dataset[[paste("JOINTLY_rep",i, sep="")]]@cell.embeddings[,1:20]
  space <- space[ rownames(space) %in% rownames(md.subset),]
  space <- space[ match(rownames(md.subset), rownames(space)),]
  dist.mat <- Rfast::Dist(space[,1:20])
  
  # Batch ASW
  sil <- cluster::silhouette(as.numeric(factor(md.subset[,"batch_label"], labels = seq(length(unique(md.subset[,"batch_label"]))))), dmatrix = dist.mat, do.col.sort = FALSE)
  sil <- sil[,"sil_width"]
  sil <- abs(sil)
  avg.sil <- c()
  for (cl in unique(md.subset[,"transfer_label"])) {
    avg.sil <- c(avg.sil, sum(1 - sil[ which(md.subset[,"transfer_label"] == cl)]) / length( which(md.subset[,"transfer_label"] == cl)))
  }
  global_bASW <- c(global_bASW, mean(avg.sil))
  
  # Cell-type ASW
  sil <- cluster::silhouette(as.numeric(factor(md.subset[,"transfer_label"], labels = seq(length(unique(md.subset[,"transfer_label"]))))), dmatrix = dist.mat, do.col.sort = FALSE)
  sil <- sil[,"sil_width"]
  global_cASW <- c(global_cASW, mean((sil + 1) / 2))
  
  # ARI
  dend <- HGC.dendrogram(G = HGC::SNN.Construction(dataset[[paste("JOINTLY_rep",i, sep="")]]@cell.embeddings[,1:20]))
  aris <- c()
  for (cl in seq(1, 50, 1)) {
    # Cluster
    if (cl == 1) {
      clusters <- rep(1, nrow(md.subset))
      names(clusters) <- rownames(md.subset)
    } else {
      clusters <- cutree(dend, k = cl)
      names(clusters) <- rownames(md.subset)
    }

    # Capture the ARI
    aris <- c(aris, aricode::ARI(clusters, factor(md.subset[,"transfer_label"])))
  }
  global_ARI <- c(global_ARI, max(aris))
}

# Results
results[results.counter,1] <- 0
results[results.counter,2] <- "None"
results[results.counter,3] <- global_ARI[ which.max(global_ARI)]
results[results.counter,4] <- global_iLISI[ which.max(global_ARI)]
results[results.counter,5] <- global_cLISI[ which.max(global_ARI)]
results[results.counter,6] <- global_bASW[ which.max(global_ARI)]
results[results.counter,7] <- global_cASW[ which.max(global_ARI)]
results.counter <- results.counter + 1

## Subsets
for (obj in tracking[ tracking[,1] == file,2]) {
  dataset.subset <- readRDS(obj)
  
  # Get labels
  labels <- read.delim(paste(substr(file, 0, nchar(file)- regexpr("/", intToUtf8(rev(utf8ToInt(file))))+1), "Transferred_labels.tsv", sep=""), header=TRUE)
  labels <- labels[ match(colnames(dataset.subset), labels$X),]
  dataset.subset$transfer_label <- labels[,2]
  dataset.subset <- subset(dataset.subset, transfer_label != "")
  dataset.subset <- subset(dataset.subset, transfer_label %in% names(which(table(dataset$transfer_label) >= 10)))
  
  ## LISI and ASW
  global_cLISI <- global_iLISI <- global_bASW <- global_cASW <- global_ARI <- c()
  md.subset <- dataset.subset@meta.data
  for (i in 1:5) {
    # LISI
    LISIS <- compute_lisi(dataset.subset[[paste("JOINTLY_rep",i, sep="")]]@cell.embeddings[,1:20], dataset.subset@meta.data, label_colnames = c("batch_label", "transfer_label"))
    global_cLISI <- c(global_cLISI, median((length(unique(dataset.subset@meta.data[,"transfer_label"]))-LISIS[,2]) / (length(unique(dataset.subset@meta.data[,"transfer_label"])) - 1)))
    global_iLISI <- c(global_iLISI, median((LISIS[,1]-1) / (length(unique(dataset.subset@meta.data[,"batch_label"])) - 1)))
    
    # Distance mat
    space <- dataset.subset[[paste("JOINTLY_rep",i, sep="")]]@cell.embeddings[,1:20]
    space <- space[ rownames(space) %in% rownames(md.subset),]
    space <- space[ match(rownames(md.subset), rownames(space)),]
    dist.mat <- Rfast::Dist(space[,1:20])
    
    # Batch ASW
    sil <- cluster::silhouette(as.numeric(factor(md.subset[,"batch_label"], labels = seq(length(unique(md.subset[,"batch_label"]))))), dmatrix = dist.mat, do.col.sort = FALSE)
    sil <- sil[,"sil_width"]
    sil <- abs(sil)
    avg.sil <- c()
    for (cl in unique(md.subset[,"transfer_label"])) {
      avg.sil <- c(avg.sil, sum(1 - sil[ which(md.subset[,"transfer_label"] == cl)]) / length( which(md.subset[,"transfer_label"] == cl)))
    }
    global_bASW <- c(global_bASW, mean(avg.sil))
    
    # Cell-type ASW
    sil <- cluster::silhouette(as.numeric(factor(md.subset[,"transfer_label"], labels = seq(length(unique(md.subset[,"transfer_label"]))))), dmatrix = dist.mat, do.col.sort = FALSE)
    sil <- sil[,"sil_width"]
    global_cASW <- c(global_cASW, mean((sil + 1) / 2))
    
    # ARI
    dend <- HGC.dendrogram(G = HGC::SNN.Construction(dataset.subset[[paste("JOINTLY_rep",i, sep="")]]@cell.embeddings[,1:20]))
    aris <- c()
    for (cl in seq(1, 50, 1)) {
      # Cluster
      if (cl == 1) {
        clusters <- rep(1, nrow(md.subset))
        names(clusters) <- rownames(md.subset)
      } else {
        clusters <- cutree(dend, k = cl)
        names(clusters) <- rownames(md.subset)
      }
      
      # Capture the ARI
      aris <- c(aris, aricode::ARI(clusters, factor(md.subset[,"transfer_label"])))
    }
    global_ARI <- c(global_ARI, max(aris))
  }
  
  # Record results
  results[results.counter,1] <- tracking[ tracking[,2] == obj,3]
  results[results.counter,2] <- tracking[ tracking[,2] == obj,4]
  results[results.counter,3] <- global_ARI[ which.max(global_ARI)]
  results[results.counter,4] <- global_iLISI[ which.max(global_ARI)]
  results[results.counter,5] <- global_cLISI[ which.max(global_ARI)]
  results[results.counter,6] <- global_bASW[ which.max(global_ARI)]
  results[results.counter,7] <- global_cASW[ which.max(global_ARI)]
  results.counter <- results.counter + 1
  print(results.counter)
}

## Neighbor consistency
reference.set <- 2

# Full dataset
overlap.score <- c()
for (q in (1:5)[-reference.set]) {
  full <- dataset[[paste("JOINTLY_rep", reference.set, sep="")]]@cell.embeddings[,1:20]
  subs <- dataset[[paste("JOINTLY_rep", q, sep="")]]@cell.embeddings[,1:20]
  subs <- subs[ match(rownames(full), rownames(subs)),]
  full.nn <- RANN::nn2(full, k= 101)$nn.idx
  subs.nn <- RANN::nn2(subs, k= 101)$nn.idx
  hits <- 0
  for (row in 1:nrow(full.nn)) {
    hits <- hits + sum(full.nn[1,2:101] %in% subs.nn[1,2:101])
  }
  overlap.score <- c(overlap.score, hits / (nrow(subs.nn) * 100))
}
results[1,8] <- max(overlap.score)

# Loop across different reference.sets
full.overlap <- c()
for (ref.set in 1:6) {
  overlap.score <- c()
  for (q in (1:6)[-ref.set]) {
    full <- dataset[[paste("JOINTLY_rep", ref.set, sep="")]]@cell.embeddings[,1:20]
    subs <- dataset[[paste("JOINTLY_rep", q, sep="")]]@cell.embeddings[,1:20]
    subs <- subs[ match(rownames(full), rownames(subs)),]
    full.nn <- RANN::nn2(full, k= 101)$nn.idx
    subs.nn <- RANN::nn2(subs, k= 101)$nn.idx
    hits <- 0
    for (row in 1:nrow(full.nn)) {
      hits <- hits + sum(full.nn[1,2:101] %in% subs.nn[1,2:101])
    }
    overlap.score <- c(overlap.score, hits / (nrow(subs.nn) * 100))
  }
  full.overlap <- c(full.overlap, max(overlap.score))
}
full.mean <- mean(full.overlap)
full.sem <- sd(full.overlap)/sqrt(6)

# Subsamples
for (obj in tracking[ tracking[,1] == file,2]) {
  # Import data
  dataset.subset <- readRDS(obj)
  
  # Get labels
  labels <- read.delim(paste(substr(file, 0, nchar(file)- regexpr("/", intToUtf8(rev(utf8ToInt(file))))+1), "Transferred_labels.tsv", sep=""), header=TRUE)
  labels <- labels[ match(colnames(dataset.subset), labels$X),]
  dataset.subset$transfer_label <- labels[,2]
  dataset.subset <- subset(dataset.subset, transfer_label != "")
  dataset.subset <- subset(dataset.subset, transfer_label %in% names(which(table(dataset$transfer_label) >= 10)))
  
  # Subset full dataset
  dataset.subsetted <- subset(dataset, cells = colnames(dataset.subset))
  
  # Find neighbors in the full dataset
  overlap.score <- c()
  for (q in 1:5) {
    full <- dataset.subsetted[[paste("JOINTLY_rep", reference.set, sep="")]]@cell.embeddings[,1:20]
    subs <- dataset.subset[[paste("JOINTLY_rep", q, sep="")]]@cell.embeddings[,1:20]
    subs <- subs[ match(rownames(full), rownames(subs)),]
    full.nn <- RANN::nn2(full, k= 101)$nn.idx
    subs.nn <- RANN::nn2(subs, k= 101)$nn.idx
    hits <- 0
    for (row in 1:nrow(full.nn)) {
      hits <- hits + sum(full.nn[1,2:101] %in% subs.nn[1,2:101])
    }
    overlap.score <- c(overlap.score, hits / (nrow(subs.nn) * 100))
  }
  results[ results[,2] == tracking[ tracking[,2] == obj,4],8] <- max(overlap.score)
  print(obj)
}

## Calculate averages and SEMs
agg.mean <- aggregate(results[,c(3,4,5,6,7,8)], by = list(results[,1]), FUN="mean")
agg.sem <- aggregate(results[,c(3,4,5,6,7,8)], by = list(results[,1]), FUN="sd")
agg.sem[2,2:7] <- agg.sem[2,2:7] / sqrt(nrow(results[ results[,1] == 1,]))
agg.sem[3,2:7] <- agg.sem[3,2:7] / sqrt(nrow(results[ results[,1] == 2,]))
agg.sem[4,2:7] <- agg.sem[4,2:7] / sqrt(nrow(results[ results[,1] == 3,]))
agg.sem[4,2:7] <- agg.sem[4,2:7] / sqrt(nrow(results[ results[,1] == 3,]))

## Plots
par(mfcol=c(1,6))
plot(x = c(0,1,2,3), y = agg.mean[,2], ylim=c(0,1), type="o", las = 2, ylab="ARI", pch=16, xlab= "# Datasets removed", xaxt = "n")
axis(side = 1, at = c(0,1,2,3), labels = c(0,1,2,3))
arrows(c(0,1,2,3), agg.mean[,2] - agg.sem[,2], c(0,1,2,3), agg.mean[,2] + agg.sem[,2], angle = 90, code = 3, length = 0.05)

plot(x = c(0,1,2,3), y = agg.mean[,3], ylim=c(0,1), type="o", las = 2, ylab="Global iLISI", pch=16, xlab= "# Datasets removed", xaxt = "n")
axis(side = 1, at = c(0,1,2,3), labels = c(0,1,2,3))
arrows(c(0,1,2,3), agg.mean[,3] - agg.sem[,3], c(0,1,2,3), agg.mean[,3] + agg.sem[,3], angle = 90, code = 3, length = 0.05)

plot(x = c(0,1,2,3), y = agg.mean[,5], ylim=c(0,1), type="o", las = 2, ylab="Global bASW", pch=16, xlab= "# Datasets removed", xaxt = "n")
axis(side = 1, at = c(0,1,2,3), labels = c(0,1,2,3))
arrows(c(0,1,2,3), agg.mean[,5] - agg.sem[,5], c(0,1,2,3), agg.mean[,5] + agg.sem[,5], angle = 90, code = 3, length = 0.05)

plot(x = c(0,1,2,3), y = agg.mean[,4], ylim=c(0,1), type="o", las = 2, ylab="Global cLISI", pch=16, xlab = "# Datasets removed", xaxt = "n")
axis(side = 1, at = c(0,1,2,3), labels = c(0,1,2,3))
arrows(c(0,1,2,3), agg.mean[,4] - agg.sem[,4], c(0,1,2,3), agg.mean[,4] + agg.sem[,4], angle = 90, code = 3, length = 0.05)

plot(x = c(0,1,2,3), y = agg.mean[,6], ylim=c(0,1), type="o", las = 2, ylab="Global cASW", pch=16, xlab = "# Datasets removed", xaxt = "n")
axis(side = 1, at = c(0,1,2,3), labels = c(0,1,2,3))
arrows(c(0,1,2,3), agg.mean[,6] - agg.sem[,6], c(0,1,2,3), agg.mean[,6] + agg.sem[,6], angle = 90, code = 3, length = 0.05)

plot(x = c(0,1,2,3), y = c(full.mean, agg.mean[2:4,7]), ylim=c(0,1), type="o", las = 2, ylab="Neighborhood preservation", pch=16, xlab = "# Datasets removed", xaxt = "n")
axis(side = 1, at = c(0,1,2,3), labels = c(0,1,2,3))
arrows(c(0,1,2,3), c(full.mean - full.sem, agg.mean[2:4,7] - agg.sem[2:4,7]), c(0,1,2,3), c(full.mean + full.sem, agg.mean[2:4,7] + agg.sem[2:4,7]), angle = 90, code = 3, length = 0.05)
