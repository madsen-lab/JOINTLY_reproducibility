##############################
###### INSTALL PACKAGES ######
##############################

pkgs <- as.data.frame(installed.packages())
if (!any(grepl("remotes", pkgs$Package))) { install.packages("remotes") }
if (!any(grepl("clevr", pkgs$Package))) { install.packages("clevr") }
if (!any(grepl("aricode", pkgs$Package))) { install.packages("aricode") }
if (!any(grepl("Seurat", pkgs$Package))) { install.packages("Seurat") }
if (!any(grepl("harmony", pkgs$Package))) { install.packages("harmony") }
if (!any(grepl("Rfast", pkgs$Package))) { install.packages("Rfast") }
if (!any(grepl("inflection", pkgs$Package))) { install.packages("inflection") }
if (!any(grepl("scCustomize", pkgs$Package))) { install.packages("scCustomize") }
if (!any(grepl("enrichR", pkgs$Package))) { install.packages("enrichR") }
if (!any(grepl("writexl", pkgs$Package))) { install.packages("writexl") }
if (!any(grepl("gamlss", pkgs$Package))) { install.packages("gamlss") }
if (!any(grepl("emmeans", pkgs$Package))) { install.packages("emmeans") }
if (!any(grepl("BisqueRNA", pkgs$Package))) { install.packages("BisqueRNA") }
if (!any(grepl("compositions", pkgs$Package))) { install.packages("compositions") }
if (!any(grepl("recount3", pkgs$Package))) { BiocManager::install("recount3", update = FALSE, ask = FALSE) }
if (!any(grepl("batchelor", pkgs$Package))) { BiocManager::install("batchelor", update = FALSE, ask = FALSE) }
if (!any(grepl("limma", pkgs$Package))) { BiocManager::install("limma", update = FALSE, ask = FALSE) }
if (!any(grepl("HGC", pkgs$Package))) { BiocManager::install("HGC", update = FALSE, ask = FALSE) }
if (!any(grepl("UCell", pkgs$Package))) { BiocManager::install("UCell", update = FALSE, ask = FALSE) }
if (!any(grepl("ComplexHeatmap", pkgs$Package))) { BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE) }
if (!any(grepl("transformGamPoi", pkgs$Package))) { BiocManager::install("transformGamPoi", update = FALSE, ask = FALSE) }
if (!any(grepl("rliger", pkgs$Package))) { remotes::install_github("welch-lab/liger", upgrade = FALSE) }
if (!any(grepl("SeuratWrappers", pkgs$Package))) { remotes::install_github("satijalab/seurat-wrappers", upgrade = FALSE) }
if (!any(grepl("lisi", pkgs$Package))) { remotes::install_github("immunogenomics/lisi", upgrade = FALSE) }
if (!any(grepl("JOINTLY", pkgs$Package))) { remotes::install_github("madsen-lab/rJOINTLY", upgrade = FALSE) }
if (!any(grepl("SeuratDisk", pkgs$Package))) { remotes::install_github("mojaveazure/seurat-disk", upgrade = FALSE) }
if (!any(grepl("presto", pkgs$Package))) { remotes::install_github("immunogenomics/presto", upgrade = FALSE) }
rm(pkgs)