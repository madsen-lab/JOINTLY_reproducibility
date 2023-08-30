simulateNoise <- function(data, label_var, mean, sd = 0.5, nbatches = 2) {
  # Extract raw counts
  counts <- GetAssayData(data, "counts")
  counts <- as.matrix(counts)
  
  # Define matrices for noise 
  noise_counts <- counts

  # Generate replicate datasets
  noise <- data

  # Define batches
  data$batch <- "batch1"
  noise$batch <- "batch1"

  # Loop across cell types, introduce batches and noise  
  for (celltype in as.character(unique(data@meta.data[, label_var]))) {
    idx <- which(data@meta.data[,label_var] == celltype)
    for (batch in 2:nbatches) {
      ## Set batch labels
      selected.idx <- sample(idx, floor(length(idx) / nbatches))
      idx <- idx[!(idx %in% selected.idx)]
      noise@meta.data[selected.idx, "batch"] <- paste("batch", batch, sep="")
      data@meta.data[selected.idx, "batch"] <- paste("batch", batch, sep="")
      
      ## Introduce noise
      # Draw a max noise level
      drawn_max <- -100
      while (drawn_max < 0) { 
        drawn_max <- rnorm(n = 1, mean = mean, sd = sd)
      }
      
      # Draw a per-gene noise fraction
      per_gene_fraction <- runif(n = nrow(data), min = 0, max = drawn_max)
      
      # Calculate per-gene noise total
      per_gene_total <- rowMeans(counts) * per_gene_fraction
      
      # Add Poisson noise
      noise <- sapply(per_gene_total,FUN = function(x) { rpois(length(selected.idx),x)})
      noise_counts[,selected.idx] <- noise_counts[,selected.idx] + t(noise)
    }
  }
  
  # Insert noised counts into the Seurat objects
  noise@assays$RNA@counts <- as(noise_counts, "dgCMatrix")
  
  # Return
  return(list(no = data, noise = noise))
}

embedMethods <- function(data, batch_var, outpath, reps = 5, nfeat = 1000) {
  ## Packages
  stop <- "No"
  pkg <- .packages(all.available = TRUE)
  required.packages <- c("Rfast","cluster","aricode","clevr","Seurat","JOINTLY","scry","harmony","rliger","lisi","batchelor", "SeuratWrappers","BiocParallel")
  if (any(!(required.packages %in% pkg))) {
    message(paste("The following packages are missing: ", paste0(required.packages[!(required.packages %in% pkg)], collapse = ", "), sep=""))
    stop <- "Yes"
  }
  
  if (stop == "No") {
    ## Load packages
    out <- sapply(required.packages, FUN = function(x) { suppressMessages(suppressWarnings(require(x, warn.conflicts = FALSE, character.only = TRUE, quietly = TRUE))) } )
    final.data <- data
    
    ## Setup to capture results
    if (file.exists(outpath)) {
      embed.list <- readRDS(outpath)
      sel.features <- embed.list$sel.featured
    } else {
      ## Finding variable genes (JOINTLY-style)
      message("Finding variable genes.")
      dataset.list <- JOINTLY::preprocess(data, batch_var)
      features <- as.data.frame(scry::devianceFeatureSelection(dataset.list[[1]]))
      colnames(features)[1] <- "Dataset1"
      for (ds in 2:length(dataset.list)) {
        feats <- scry::devianceFeatureSelection(dataset.list[[ds]])
        features <- merge(as.data.frame(features), as.data.frame(feats), by=0)
        colnames(features)[ncol(features)] <- paste("Dataset", ds, sep="")
        rownames(features) <- features[,1]
        features <- features[,-1]
      }
      for (col in 1:ncol(features)) { features <- features[ !is.na(features[,col]),] }
      features$average <- apply(features[,c(1:ncol(features))],1,FUN="mean")
      features$count <- 0
      for (col in 1:(ncol(features)-2)) { features[ which(rank(-features[,col]) <= nfeat),"count"] <- features[ which(rank(-features[,col]) <= nfeat),"count"] + 1 }
      features <- features[ order(features$count, features$average, decreasing = TRUE),]
      sel.features <- rownames(features[ features$count > 0,])
      embed.list <- list(embeddings = list(), method = 1, sel.featured = sel.features)
    }
    
    ## Running JOINTLY
    # Setup and message
    if (embed.list$method == 1) {
      message("Running JOINTLY.")
      object.list <- embed.list$embeddings
      
      # Run the initial steps of the pipeline
      seu <- data
      proc <- preprocess(seu, batch_var)
      cpca <- cpca(proc)
      inputs <- prepareData(cpca$cpca)
      
      # Run method and evaluate the embedding
      for (rep in 1:reps) {
        message(paste("\tRunning replicate ", rep, sep=""))
        # Embed and evaluate
        solved <- JOINTLYsolve(inputs$kernels, inputs$snn, inputs$rareity, cpca, bpparam = BiocParallel::MulticoreParam(), progressbar = FALSE)
        H <- t(do.call("cbind",solved$Hmat))
        H <- H[ match(colnames(data), rownames(H)),]
        H <- scale(H)
        H <- t(scale(t(H)))
        colnames(H) <- paste("JOINTLY", 1:ncol(H), sep="_")
        
        # Capture objects
        object.list[[length(object.list)+1]] <- H
        names(object.list)[length(object.list)] <- paste("JOINTLY_rep", rep, sep="")
      }
      
      # Save intermediate results
      embed.list$embeddings <- object.list
      embed.list$method <- 2
      saveRDS(embed.list, outpath)
    }
    
    ## Running Harmony
    # Setup and message
    if (embed.list$method == 2) {
      message("Running Harmony.")
      object.list <- embed.list$embeddings
      
      # Create a Seurat object and run PCA
      seu <- data
      seu <- Seurat::NormalizeData(seu, verbose = FALSE)
      VariableFeatures(seu) <- sel.features
      seu <- Seurat::ScaleData(seu, verbose = FALSE)
      seu <- Seurat::RunPCA(seu, verbose = FALSE, npcs = 20)
      
      # Run method and evaluate the embedding
      for (rep in 1:reps) {
        message(paste("\tRunning replicate ", rep, sep=""))
        # Embed and evaluate
        tmp <- suppressWarnings(harmony::RunHarmony(seu, group.by.vars = batch_var, verbose = FALSE))
        tmp.embedding <- Seurat::Embeddings(tmp, "harmony")
        colnames(tmp.embedding) <- paste("Harmony_", seq(1, ncol(tmp.embedding), 1), sep="")
        
        # Capture objects
        object.list[[length(object.list)+1]] <- tmp.embedding
        names(object.list)[length(object.list)] <- paste("Harmony_rep", rep, sep="")
      }
      
      # Save intermediate results
      embed.list$embeddings <- object.list
      embed.list$method <- 3
      saveRDS(embed.list, outpath)
    }
    
    ## Running Seurat
    # Setup and message
    if (embed.list$method == 3) {
      message("Running Seurat using RPCA")
      object.list <- embed.list$embeddings
      
      # Create a Seurat object and run PCA
      seu <- data
      split.object <- Seurat::SplitObject(seu, batch_var)
      split.object <- lapply(split.object, FUN = function(x) { 
        x <- Seurat::NormalizeData(x, verbose = FALSE)
        x <- Seurat::ScaleData(x, verbose = FALSE, features = sel.features)
        x <- suppressWarnings(Seurat::RunPCA(x, verbose = FALSE, features = sel.features))
      })
      
      # Run method and evaluate the embedding
      for (rep in 1:reps) {
        message(paste("\tRunning replicate ", rep, sep=""))
        # Embed and evaluate
        out <- capture.output(anchors <- Seurat::FindIntegrationAnchors(split.object, anchor.features = sel.features, reduction = "rpca", verbose = FALSE))
        tmp <- IntegrateData(anchorset = anchors, verbose = FALSE, k.weight = ifelse(min(unlist(lapply(split.object, "ncol"))) < 100, 50, 100))
        DefaultAssay(tmp) <- "integrated"
        tmp <- ScaleData(tmp, verbose = FALSE)
        tmp <- RunPCA(tmp, npcs = 30, verbose = FALSE)
        tmp.embedding <- Seurat::Embeddings(tmp, "pca")
        colnames(tmp.embedding) <- paste("RPCA_", seq(1, ncol(tmp.embedding), 1), sep="")
        
        # Capture objects
        object.list[[length(object.list)+1]] <- tmp.embedding
        names(object.list)[length(object.list)] <- paste("RPCA_rep", rep, sep="")
      }
      
      # Save intermediate results
      embed.list$embeddings <- object.list
      embed.list$method <- 4
      saveRDS(embed.list, outpath)
    }
    
    ## Running liger
    # Setup and message
    if (embed.list$method == 4) {
      message("Running LIGER")
      object.list <- embed.list$embeddings
      
      # Create a Seurat object and run PCA
      seu <- data
      seu <- Seurat::NormalizeData(seu, verbose = FALSE)
      VariableFeatures(seu) <- sel.features
      seu <- ScaleData(seu, split.by = batch_var, do.center = FALSE, verbose = FALSE)
      
      # Run method and evaluate the embedding
      for (rep in 1:reps) {
        message(paste("\tRunning replicate ", rep, sep=""))
        # Embed and evaluate
        out <- capture.output(tmp <- suppressWarnings(SeuratWrappers::RunOptimizeALS(seu, k = 20, lambda = 5, split.by = batch_var, verbose = FALSE)))
        tmp <- suppressWarnings(SeuratWrappers::RunQuantileNorm(tmp, split.by = batch_var, verbose = FALSE))
        tmp.embedding <- Seurat::Embeddings(tmp, "iNMF")
        colnames(tmp.embedding) <- paste("LIGER_", seq(1, ncol(tmp.embedding), 1), sep="")
        
        # Capture objects
        object.list[[length(object.list)+1]] <- tmp.embedding
        names(object.list)[length(object.list)] <- paste("LIGER_rep", rep, sep="")
      }
      
      # Save intermediate results
      embed.list$embeddings <- object.list
      embed.list$method <- 5
      saveRDS(embed.list, outpath)
    }
    
    ## Running FastMNN
    # Setup and message
    if (embed.list$method == 5) {
      message("Running FastMNN")
      object.list <- embed.list$embeddings
      
      # Run method and evaluate the embedding
      for (rep in 1:reps) {
        message(paste("\tRunning replicate ", rep, sep=""))
        # Embed and evaluate
        sce <- Seurat::as.SingleCellExperiment(data)
        sce_corrected <- batchelor::fastMNN(sce, batch = data@meta.data[,batch_var], auto.merge = TRUE, subset.row = sel.features)
        tmp <- data
        corrected <- SingleCellExperiment::reducedDim(sce_corrected, "corrected")
        colnames(corrected) <- paste("FastMNN_", seq(1, ncol(corrected)), sep="")
        
        # Capture objects
        object.list[[length(object.list)+1]] <- corrected
        names(object.list)[length(object.list)] <- paste("FastMNN_rep", rep, sep="")
      }
      
      # Save intermediate results
      embed.list$embeddings <- object.list
      embed.list$method <- 6
      saveRDS(embed.list, outpath)
    }
    
    ## Running unintegrated
    # Setup and message
    if (embed.list$method == 6) {
      message("Running unintegrated")
      object.list <- embed.list$embeddings
      
      # Setup
      seu <- data
      seu <- Seurat::NormalizeData(seu, verbose = FALSE)
      Seurat::VariableFeatures(seu) <- sel.features
      seu <- Seurat::ScaleData(seu, verbose = FALSE)
      
      # Run method and evaluate the embedding
      for (rep in 1:reps) {
        message(paste("\tRunning replicate ", rep, sep=""))
        # Embed and evaluate
        tmp <- Seurat::RunPCA(seu, seed.use = NULL, verbose = FALSE)
        tmp <- tmp@reductions$pca@cell.embeddings
        colnames(tmp) <- paste("Unintegrated_", seq(1, ncol(tmp)), sep="")
        
        # Capture objects
        object.list[[length(object.list)+1]] <- tmp
        names(object.list)[length(object.list)] <- paste("Unintegrated_rep", rep, sep="")
      }
      
      # Save intermediate results
      embed.list$embeddings <- object.list
      embed.list$method <- 7
      saveRDS(embed.list, outpath)
    }
    
    ## Running scVI and Scanorama
    # Setup and message
    if (embed.list$method == 7) {
      object.list <- embed.list$embeddings
      
      # Save objects
      tmp <- data
      VariableFeatures(tmp) <- NULL
      tmp$batch <- tmp@meta.data[,batch_var]
      tmp <- Seurat::DietSeurat(tmp)
      suppressMessages(SeuratDisk::SaveH5Seurat(tmp, "tmp.h5seurat", verbose = FALSE, overwrite = TRUE))
      suppressMessages(Convert("tmp.h5seurat", "h5ad", verbose = FALSE, overwrite = TRUE))
      write.table(sel.features, "tmp.features", col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
      
      # Call the python script
      system2(command = "bash", args = "./RunPythonMethods.sh", wait = FALSE)
      
      # Import the results
      done <- 0
      curr.rep <- 0
      first.method <- 1
      first.replicate <- 1
      methods <- c("Scanorama","scVI")
      method <- methods[1]
      methods <- methods[ !(methods %in% method)]
      while (done == 0) {
        # Print method if first time checking
        if (first.method == 1) {
          message(paste("Running ", method, sep=""))
          first.method = 0
        }
        
        # Print replicate if first time check
        if (first.replicate == 1) {
          message(paste("\tRunning replicate ", curr.rep+1, sep=""))
          first.replicate <- 0
        }
        
        # Check for file
        if (file.exists(paste("tmp_", method, "_rep", curr.rep, ".txt", sep=""))) {
          # Import data
          Sys.sleep(10)
          tmp.embed <- read.delim(paste("tmp_", method, "_rep", curr.rep, ".txt", sep=""), sep=",")
          rownames(tmp.embed) <- tmp.embed[,1]
          tmp.embed <- tmp.embed[,-1]
          colnames(tmp.embed) <- paste(method, "_", seq(1, ncol(tmp.embed),1), sep="")
          tmp.embed <- tmp.embed[ match(colnames(data), rownames(tmp.embed)),]
          object.list[[length(object.list)+1]] <- tmp.embed
          names(object.list)[length(object.list)] <- paste(method, "_rep", curr.rep+1, sep="")
          tmp <- file.remove(paste("tmp_", method, "_rep", curr.rep, ".txt", sep=""))
          
          if ((curr.rep+1) < reps) {
            curr.rep <- curr.rep + 1
            first.replicate <- 1
          } else {
            if (length(methods) == 0) {
              done <- 1
            } else {
              method <- methods[1]
              methods <- methods[ !(methods %in% method)]
              first.method <- 1
              curr.rep <- 0
              first.replicate <- 1
            }
          }
        } else {
          Sys.sleep(1)
        }
      }
      
      # Cleanup
      tmp <- file.remove("tmp.h5seurat")
      tmp <- file.remove("tmp.h5ad")
      tmp <- file.remove("tmp.features")
      tmp <- file.remove("tmp.err.log")
      tmp <- file.remove("tmp.out.log")
      
      # Save intermediate results
      embed.list$embeddings <- object.list
      embed.list$method <- 8
      saveRDS(embed.list, outpath)
    }
    
    ## Define output
    output <- embed.list
  } else {
    output <- NULL
  }
  
  # Return
  return(output)
}


summarizeResults <- function(result, data, embeddings, dim.list, selectRun = "globalARI") {
  summary <- do.call("rbind", lapply(result, FUN = function(x) { x$summary}))
  summary$Method <- substr(rownames(summary),0,regexpr("_", rownames(summary))-1)
  idx.keep <- c()
  for (method in c("JOINTLY","Harmony", "RPCA","LIGER","FastMNN","scVI","Scanorama","Unintegrated")) {
    idx <- which(summary$Method == method & summary[,selectRun] == max(summary[ summary$Method == method, selectRun]))[1]
    idx.keep <- c(idx.keep, idx)
    embed <- embeddings[[idx]]
    embed <- embed[ rownames(embed) %in% colnames(data),]
    embed <- embed[ match(colnames(data), rownames(embed)),]
    data[[method]] <- Seurat::CreateDimReducObject(as.matrix(embed), assay = "RNA")
    clusters <- result[[idx]]$clusters[[which.max(result[[idx]]$cluster_metrics[,selectRun])]]
    clusters <- clusters[ names(clusters) %in% colnames(data)]
    clusters <- clusters[ match(colnames(data), names(clusters))]
    data@meta.data[,(ncol(data@meta.data)+1)] <- clusters
    colnames(data@meta.data)[ncol(data@meta.data)] <- paste(method, "_clusters", sep="")
    dims <- dim.list[[which(names(dim.list) == method)]]
    data <- RunUMAP(data, dims = 1:dims, reduction = method, reduction.name = paste(method, "_UMAP", sep=""), verbose = FALSE)
  }
  summary <- summary[idx.keep,]
  out <- list(summary = summary, object = data)
  return(out)
}

evaluateEmbedding <- function(x, metadata, batch_var, label_var, dim.list, cl.min = 1, cl.max = 50) {
  # Insert embedding in the Seurat object
  x <- x[ rownames(x) %in% rownames(metadata),]
  x <- x[ match(rownames(metadata), rownames(x)),]
  x <- as.matrix(x)
  
  # Find which number of dimensions to use
  dims <- dim.list[[which(names(dim.list) == substr(colnames(x)[1], 0, regexpr("_", colnames(x)[1])-1))]]
  
  # Setup to capture clustering metrics
  global_aris <- c()
  global_nmis <- c()
  global_vs <- c()
  global_completeness <- c()
  global_homogeneity <- c()
  dataset_aris <- c()
  dataset_nmis <- c()
  dataset_vs <- c()
  dataset_completeness <- c()
  dataset_homogeneity <- c()
  cluster_list <- list()
  
  # Evaluate clustering metrics
  dend <- HGC.dendrogram(G = HGC::SNN.Construction(x[,1:dims]))
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
    
    # Capture global metrics
    global_aris <- c(global_aris, aricode::ARI(clusters, factor(metadata[,label_var])))
    global_nmis <- c(global_nmis, aricode::NMI(clusters, factor(metadata[,label_var])))
    global_vs <- c(global_vs, clevr::v_measure(clusters, factor(metadata[,label_var])))
    global_completeness <- c(global_completeness, clevr::completeness(clusters, factor(metadata[,label_var])))
    global_homogeneity <- c(global_homogeneity, clevr::homogeneity(clusters, factor(metadata[,label_var])))
    
    # Setup to capture dataset metrics
    ds_aris <- c()
    ds_nmis <- c()
    ds_vs <- c()
    ds_completeness <- c()
    ds_homogeneity <- c()
    
    # Calculate dataset metrics
    for (ds in unique(metadata[,batch_var])) {
      md <- metadata[ metadata[,batch_var] == ds,]
      clusters.ds <- clusters[ names(clusters) %in% rownames(md) ]
      clusters.ds <- clusters.ds[ match(rownames(md), names(clusters.ds))]
      ds_aris <- c(ds_aris, aricode::ARI(clusters.ds, factor(md[,label_var])))
      ds_nmis <- c(ds_nmis, aricode::NMI(clusters.ds, factor(md[,label_var])))
      ds_vs <- c(ds_vs, clevr::v_measure(clusters.ds, factor(md[,label_var])))
      ds_completeness <- c(ds_completeness, clevr::completeness(clusters.ds, factor(md[,label_var])))
      ds_homogeneity <- c(ds_homogeneity, clevr::homogeneity(clusters.ds, factor(md[,label_var])))
    }
    
    # Capture dataset metrics
    dataset_aris <- c(dataset_aris, min(ds_aris))
    dataset_nmis <- c(dataset_nmis, min(ds_nmis))
    dataset_vs <- c(dataset_vs, min(ds_vs))
    dataset_completeness <- c(dataset_completeness, min(ds_completeness))
    dataset_homogeneity <- c(dataset_homogeneity, min(ds_homogeneity))
  }
  
  ## cLISI and iLISI
  space <- x
  LISIS <- lisi::compute_lisi(space[,1:dims], metadata, c(batch_var, label_var))
  global_cLISI <- (length(unique(metadata[,label_var]))-median(LISIS[,2])) / (length(unique(metadata[,label_var])) - 1)
  global_iLISI <- (median(LISIS[,1])-1) / (length(unique(metadata[,batch_var])) - 1)
  dataset_cLISI <- c()
  dataset_iLISI <- c()
  for (ds in unique(metadata[,batch_var])) {
    dataset_cLISI <- c(dataset_cLISI, (length(unique(metadata[ metadata[,batch_var] == ds,label_var]))-median(LISIS[ rownames(LISIS) %in% rownames(metadata[ metadata[,batch_var] == ds,]),2])) / (length(unique(metadata[ metadata[,batch_var] == ds,label_var])) - 1))
    dataset_iLISI <- c(dataset_iLISI, (median(LISIS[ rownames(LISIS) %in% rownames(metadata[ metadata[,batch_var] == ds,]),1])-1) / (length(unique(metadata[,batch_var])) - 1))
  }
  
  ## Label and batch ASW
  # Subsample if too large
  if (dim(space)[1] >= 40000) {
    smps <- c()
    pr_ds_cells <- ceiling(table(metadata[,batch_var]) / (sum(table(metadata[,batch_var])) / 40000))
    for (ds in unique(metadata[,batch_var])) {
      pr_ct_cells <- ceiling(table(metadata[ metadata[,batch_var] == ds, label_var]) / (sum(table(metadata[ metadata[,batch_var] == ds, label_var]))  / pr_ds_cells[which(names(pr_ds_cells) == ds)]))
      ds.idx <- which(metadata[,batch_var] == ds)
      for (ct in unique(metadata[ metadata[,batch_var] == ds,label_var])) {
        ct.idx <- which(metadata[ metadata[,batch_var] == ds,label_var] == ct)
        selected <- sample(ct.idx, size = pr_ct_cells[names(pr_ct_cells) == ct], replace = FALSE)
        smps <- c(smps, ds.idx[selected])
      }
    }
    smps <- unique(smps)
    md.subset <- metadata[ rownames(metadata) %in% rownames(metadata)[smps],]
  } else {
    md.subset <- metadata
  }
  
  # Distance matrix
  space <- x
  space <- space[ rownames(space) %in% rownames(md.subset),]
  space <- space[ match(rownames(md.subset), rownames(space)),]
  dist.mat <- Rfast::Dist(space[,1:dims])
  
  # Batch ASW
  sil <- cluster::silhouette(as.numeric(factor(md.subset[,batch_var], labels = seq(length(unique(md.subset[,batch_var]))))), dmatrix = dist.mat, do.col.sort = FALSE)
  sil <- sil[,"sil_width"]
  sil <- abs(sil)
  avg.sil <- c()
  for (cl in unique(md.subset[,label_var])) {
    avg.sil <- c(avg.sil, sum(1 - sil[ which(md.subset[,label_var] == cl)]) / length( which(md.subset[,label_var] == cl)))
  }
  global.basw <- mean(avg.sil)
  dataset.basw <- c()
  for (ds in unique(md.subset[,batch_var])) {
    avg.sil <- c()
    for (cl in unique(md.subset[md.subset[,batch_var] == ds,label_var])) {
      avg.sil <- c(avg.sil, sum(1 - sil[ which(md.subset[,label_var] == cl & md.subset[,batch_var] == ds)]) / length( which(md.subset[,label_var] == cl & md.subset[,batch_var] == ds)))
    }
    dataset.basw <- c(dataset.basw, mean(avg.sil))
  }
  
  # Cell type ASW
  sil <- cluster::silhouette(as.numeric(factor(md.subset[,label_var], labels = seq(length(unique(md.subset[,label_var]))))), dmatrix = dist.mat, do.col.sort = FALSE)
  sil <- sil[,"sil_width"]
  global.casw <- mean((sil + 1) / 2)
  dataset.casw <- c()
  for (ds in unique(md.subset[,batch_var])) {
    dataset.casw <- c(dataset.casw, mean((sil[ which(md.subset[,batch_var] == ds)] + 1) / 2))
  }
  
  # Format results
  asw <- data.frame(celltype = dataset.casw, batch = dataset.basw)
  lisi <- data.frame(celltype = dataset_cLISI, batch = dataset_iLISI)
  clustering_metrics <- data.frame(globalARI = global_aris, globalNMI = global_nmis, globalVmeasure = global_vs, globalHomogeneity = global_homogeneity, globalCompleteness = global_completeness, datasetARI = dataset_aris, datasetNMI = dataset_nmis, datasetVmeasure = dataset_vs, datasetHomogeneity = dataset_homogeneity, datasetCompleteness = dataset_completeness)
  summary <- data.frame(globalARI = max(global_aris), datasetARI = max(dataset_aris), globalNMI = max(global_nmis), datasetNMI = max(dataset_nmis), globalVmeasure = max(global_vs), datasetVmeasure = max(dataset_vs), globalHomogeneity = max(global_homogeneity), datasetHomogeneity = max(dataset_homogeneity), globalCompleteness = max(global_completeness), datasetCompletness = max(dataset_completeness), globalcLISI = global_cLISI, datasetcLISI = min(dataset_cLISI), globaliLISI = global_iLISI, datasetiLISI = min(dataset_iLISI), globalcASW = global.casw, datasetcASW = min(dataset.casw), globalbASW = global.basw, datasetbASW = min(dataset.basw))
  
  # Return
  list(summary = summary, cluster_metrics = clustering_metrics, ASW = asw, LISI = lisi, clusters = cluster_list)
}

getModules <- function(solved) {
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
  
  # Sum the W matrix
  W.sum <- W[[1]]
  for (i in 2:length(W)) { W.sum <- W.sum + W[[i]]}
  W.sum <- W.sum / 2
  
  # Scale the sum matrix
  W.tmp <- W.sum
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W.tmp <- scale(W.tmp)
  W.tmp <- t(W.tmp)
  W.sum <- W.tmp
  
  # Get modules
  modules <- list()
  for (i in 1:ncol(W.sum)) {
    modules[[length(modules)+1]] <- names(sort(W.sum[,i], decreasing = TRUE))[1:inflection::uik(y = sort(W.sum[,i], decreasing = TRUE), x = seq(1, nrow(W.sum),1))]
    names(modules)[length(modules)] <- paste("factor_", i, sep="")
  }
  
  # Return
  return(modules)
}

process_summary <- function(i, label) {
  kidney_res = cbind(kidney$summary[metric[i]], kidney$summary$Method )
  liver_res = cbind(liver$summary[metric[i]], liver$summary$Method )
  pancreas_res = cbind(pancreas$summary[metric[i]], pancreas$summary$Method )
  lung_res = cbind(lung$summary[metric[i]], lung$summary$Method )
  pbmc_res = cbind(pbmc$summary[metric[i]], pbmc$summary$Method )
  
  kidney_res[,1] <- as.numeric(kidney_res[,1])
  liver_res[,1] = as.numeric(liver_res[,1])
  pancreas_res[,1] = as.numeric(pancreas_res[,1])
  lung_res[,1] = as.numeric(lung_res[,1])
  pbmc_res[,1] = as.numeric(pbmc_res[,1])
  
  colnames(kidney_res) = c('kidney', 'Method')
  colnames(liver_res) = c('liver', 'Method')
  colnames(pancreas_res) = c('pancreas', 'Method')
  colnames(lung_res) = c('lung', 'Method')
  colnames(pbmc_res) = c('pbmc', 'Method')
  
  merge = merge(kidney_res, liver_res, by = 'Method')
  merge = merge(merge, pancreas_res, by = 'Method')
  merge = merge(merge, lung_res, by = 'Method')
  merge = merge(merge, pbmc_res, by = 'Method')
  
  rownames(merge) = merge$Method
  merge = merge[2:6]
  rows = rownames(merge)
  merge = apply(merge, 2, as.numeric) 
  rownames(merge) = rows
  merge = t(merge)
  merge = rbind(merge, All = colSums(merge) / 5)
  merge = signif(merge,digits = 2)
  r_merge = apply(-merge, 1, function(x) rank(x, ties.method = 'min'))  
  r_merge = t(r_merge)
  
  rownames(merge) = paste(rownames(merge), label, sep=" ")
  rownames(r_merge) = paste(rownames(r_merge), label, sep=" ")
  return(list(metric = merge, 
              rank = r_merge))
}

prepareGPT <- function(dataset, batch_var, nfeat = 1000, outpath = NULL, load = TRUE, genes = "data/genes.rda") {
  # Read the dataset
  if (load) {
    message("Loading the dataset")
    dataset <- readRDS(dataset)
    dataset@meta.data[,batch_var] <- as.character(dataset@meta.data[,batch_var])
    assays.remove <- names(dataset@assays)[names(dataset@assays) != "RNA"]
    if (length(assays.remove) > 0) { for (assay in assays.remove) { dataset[[assay]] <- NULL } }
    if (is.null(outpath)) { outpath <- dataset }
  }
  
  # Find variable features
  message("Finding variable genes.")
  dataset.list <- JOINTLY::preprocess(dataset, batch_var, verbose = FALSE)
  features <- as.data.frame(scry::devianceFeatureSelection(dataset.list[[1]]))
  colnames(features)[1] <- "Dataset1"
  for (ds in 2:length(dataset.list)) {
    feats <- scry::devianceFeatureSelection(dataset.list[[ds]])
    features <- merge(as.data.frame(features), as.data.frame(feats), by=0)
    colnames(features)[ncol(features)] <- paste("Dataset", ds, sep="")
    rownames(features) <- features[,1]
    features <- features[,-1]
  }
  for (col in 1:ncol(features)) { features <- features[ !is.na(features[,col]),] }
  features$average <- apply(features[,c(1:ncol(features))],1,FUN="mean")
  features$count <- 0
  for (col in 1:(ncol(features)-2)) { features[ which(rank(-features[,col]) <= nfeat),"count"] <- features[ which(rank(-features[,col]) <= nfeat),"count"] + 1 }
  features <- features[ order(features$count, features$average, decreasing = TRUE),]
  sel.features <- rownames(features[ features$count > 0,])
  
  # Convert gene names if necessary
  load(genes)
  id <- names(which.max(apply(human, 2, FUN = function(x) { sum(rownames(dataset) %in% x)})))
  if (id != "Symbol") {
    counts <- dataset@assays$RNA@counts
    genes <- human[ human[,id] %in% rownames(counts),]
    genes <- genes[ duplicated(genes[,id])==FALSE,]
    genes <- genes[ duplicated(genes[,"Symbol"])==FALSE,]
    genes <- genes[ genes[,"Symbol"] != "",]
    counts <- counts[ rownames(counts) %in% genes[,id],]
    genes <- human[ human[,id] %in% rownames(counts),]
    genes <- genes[ match(rownames(counts), genes[,id]),]
    rownames(counts) <- genes$Symbol
    dataset[["RNA"]] <- CreateAssayObject(counts)
    sel.features <- genes[ genes[,id] %in% sel.features,"Symbol"]
  }
  
  # Export adata object
  message("Exporting anndata object.")
  tmp <- dataset
  VariableFeatures(tmp) <- NULL
  tmp$batch <- tmp@meta.data[,batch_var]
  tmp <- Seurat::DietSeurat(tmp)
  suppressMessages(SeuratDisk::SaveH5Seurat(tmp, gsub(".rds", ".h5Seurat", outpath), verbose = FALSE, overwrite = TRUE))
  suppressMessages(Convert(gsub(".rds", ".h5Seurat", outpath), "h5ad", verbose = FALSE, overwrite = TRUE))
  out <- capture.output(file.remove(gsub(".rds", ".h5Seurat", outpath)))
  write.table(sel.features, gsub(".rds", ".features", outpath), col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
}
