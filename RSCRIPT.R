# R Script accompanying the manuscript Li and Morarach et al. 2024:  
# Author: Khomgrit Morarach, Karolinska Institutet, Stockholm.
################################################################################
#                                  LIBRARIES                                   #
################################################################################
library(Seurat)
library(SeuratObject)
library(dplyr)
library(DoubletFinder)
library(clustree)
library(patchwork)
library(ggplot2)
library(reticulate)
library(cols4all)
options(future.globals.maxSize= 89128960000)
################################################################################
#                                 FUNCTIONS                                    #
################################################################################
# Utility function for loading CellRanger outputs (filtered_feature_bc_matrix)
# and creating a list of Seurat objects with percent.mt calculated.
# This function takes 2 arguments:
# (1) A list of paths to filtered_feature_bc_matrix folders,
# (2) seurat_object_version: a string for Seurat object version("v3" or "v5").
# If seurat_object_version = NULL, Assay5 (v5) will be used. 
# e.g.
# folders <- list(
#   "~/sm022/filtered_feature_bc_matrix",
#   "~/sm023/filtered_feature_bc_matrix",
#   "~/sm024/filtered_feature_bc_matrix"
# )
load_data_10x <- function(folders, seurat_object_version = "v5"){
  if (is.null(seurat_object_version)){
    options(Seurat.object.assay.version = "v5")
  }
  else if ((seurat_object_version == "v3") | (seurat_object_version == "v5")) {
    options(Seurat.object.assay.version = seurat_object_version)
  }
  else {
    cat("Set seurat_object_version to either v3 or v5", "\n")
  }
  if (is.null(names(folders))) {
    names(folders) <- paste0("sample_", seq_along(folders))
  }
  seu_list <- list()
  for (i in seq_along(folders)) {
    drdata <- Read10X(data.dir = folders[[i]])
    seuobj <- CreateSeuratObject(counts = drdata, project = names(folders)[[i]]) %>%
      PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt")
    seuobj$source_directory <- folders[[i]]
    seu_list[[i]] <- seuobj
  }
  return(seu_list)
}
# load_sm <- load_data_10x(folders)
################################################################################
# Merging datasets from same experiments.  
# ens_sm <- merge(
#   load_sm[[1]],
#   y = c(load_sm[[2]], load_sm[[3]]),
#   add.cell.ids = c("sm022", "sm023", "sm024"),
#   project = "ENSP24"
# )
# (option) Join Assay5 layers: ens_sm[["RNA"]] <- JoinLayers(ens_sm[["RNA"]])
################################################################################
# Utility function to get basic information about a Seurat object. Taking a 
# Seurat object as input and giving the object's information in the R console
report_seurat_status <- function(seu) {
  if (!"Seurat" %in% class(seu)) {
    stop("The input object is not a Seurat object.")
  }
  metadata <- seu@meta.data
  cat("==========================================================\n")
  cat("------------------ SEURAT OBJECT REPORT ------------------\n")
  cat("==========================================================\n")
  cat("Total Cells:", ncol(seu), "\n")
  cat("Total Features:", nrow(seu), "\n")
  cat("Assay Class:", class(seu[["RNA"]]), "\n")
  if ("SCT" %in% names(seu@assays)) {
    cat("Normalization method: SCTransform\n")
  }
  else if ("RNA" %in% names(seu@assays) &&
           "scale.data" %in% names(seu[["RNA"]])) {
    cat("Normalization method: LogNormalize\n")
  }
  else {
    cat("Normalization status: Unable to identify\n")
  }
  reductions <- names(seu@reductions)
  if (length(reductions) > 0) {
    cat("Dimensionality reductions:",
        paste(reductions, collapse = ", "),
        "\n")
  } else {
    cat("No dimensionality reductions present.\n")
  }
  cat("\n")
  cat("==========================================================\n")
  cat("-------------------- METADATA AVAILABLE ------------------\n")
  cat("==========================================================\n")
  if (ncol(metadata) > 0) {
    cat("Metadata columns (#unique values):\n")
    for (colname in colnames(metadata)) {
      unique_vals <- unique(metadata[[colname]])
      cat(colname, ":", length(unique_vals), "\n")
    }
  } else {
    cat("No metadata found.\n")
  }
  cat("==========================================================\n\n")
}
# report_seurat_status(ens_sm)  
################################################################################
# Given a Seurat object, this function filters cells based on the total number
# of UMI counts (nCount_RNA), the total number of genes (nFeature_RNA) and the
# percentage of genes mapped to mitochondrial genes (percent.mt), per cell.
# This function takes 4 arguments: 
# (1) A Seurat object,
# (2) nfeature_bounds, a vector of lower and upper bounds e.g. c(500, 6500),
# (3) ncount_bounds, a vector of lower and upper bounds e.g. c(200, 40000),
# (4) percent_mt_threshold, a number specifying a percent.mt threshold e.g. 20.
# The function returns a vector containing cells to keep. 
filter_cell_by_threshold <-
  function(seu,
           nfeature_bounds = c(500, 6500),
           ncount_bounds = c(0, 40000),
           percent_mt_threshold = 20) {
    seu <- PercentageFeatureSet(seu, pattern = "^mt-", col.name = "percent.mt")
    cells_to_keep <- WhichCells(
      seu,
      expression = nFeature_RNA > nfeature_bounds[1] &
        nFeature_RNA < nfeature_bounds[2] &
        percent.mt < percent_mt_threshold &
        nCount_RNA > ncount_bounds[1] &
        nCount_RNA < ncount_bounds[2]
    )
    return(cells_to_keep)
  }
# cells_to_keep <- filter_cell_by_threshold(ens_sm)
################################################################################
# Given a Seurat object, this function will cluster cells over a range of 
# resolutions using Leiden algorithm. In doing so, it will first perform  
# SCtransform(), RunPCA() and FindNeighbors() using 30 PCs.
# The function takes 2 arguments: 
# (1) A Seurat object,
# (2) reso, a vector of clustering resolutions e.g. c(0.2, 0.6, 1.0, 1.2, 1.8).
# The function returns the Seurat object with clustering results in the metadata.
get_reso_sweep_cluster <- function(seu, reso  = c(0.2, 0.6, 1.0, 1.2, 1.8)){
  seu <- SCTransform(seu)
  seu <- RunPCA(seu)
  seu <- FindNeighbors(seu, dims = 1:30)
  seu <- FindClusters(seu, algorithm = 4, resolution = reso) #or algo 1
  return(seu)
}
# ens_sm <- get_reso_sweep_cluster(ens_sm)
################################################################################
# This function iterates through the metadata with the prefix "^SCT_snn_res." and  
# For each iteration, it identifies cluster outliers, defined per cluster, as 
# cells with minimum nCount_RNA = 800 or nCount_RNA below the 1st percentile,
# or nCount_RNA above 99th percentile, or nFeature_RNA below the 1st percentile.   
# The function takes 5 arguments: 
# (1) A Seurat object,
# (2) counts_min, an integer for a minimum nCount_RNA e.g. 800,
# (3) counts_low, a number of lower bound percentile for counts   e.g 0.001,
# (4) counts_upp, a number of upper bound percentile for counts   e.g 0.099,
# (5) ngenes_low, a number of lower bound percentile for features e.g 0.001.
# The function returns a vector containing aggregated cluster outliers. 
get_cluster_outlier <- function(seu,
                                counts_min   =  800, 
                                counts_low   = 0.01,
                                counts_upp   = 0.99,
                                ngenes_low   = 0.01) {
  creso <- colnames(seu@meta.data[grep("^SCT_snn_res.", colnames(seu@meta.data))])
  cu <- c()
  for (i in 1:length(creso)) {
    split <- SplitObject(seu, split.by = creso[i])
    for (i in 1:length(split)) {
      mincov = counts_min
      if (min(split[[i]]$nCount_RNA) >= mincov) {
        counts_low = min(split[[i]]$nCount_RNA)
      } else{counts_low = quantile(split[[i]]$nCount_RNA, prob = c(0.01))}
      counts_high = quantile(split[[i]]$nCount_RNA,   prob = 0.99)
      ngenes_low  = quantile(split[[i]]$nFeature_RNA, prob = 0.01)
      ci <- WhichCells(split[[i]], invert = T,
                     expression = nFeature_RNA > ngenes_low & 
                       nCount_RNA > counts_low  &
                       nCount_RNA < counts_high)
      cu <- union(cu, ci)
    }
  }
  return(cu)
}
# outliers <- get_cluster_outlier(ens_sm)
################################################################################
# This function calculates parameters that will be used in DoubletFinder and run  
# DoubletFinder() to identify possible doublets in a Seurat object. It
# is adopted from github.com/chris-mcginnis-uscf/DoubletFinder. We set sct=TRUE 
# and reuse.pANN=FALSE, while allowing it to take arguments as follow.    
# The function takes 6 arguments: 
# (1) A Seurat object,
# (2) expected_prop, a decimal for expected doublets e.g. 0.075,
# (3) pN, a number for a proportion of generated artificial doublets e.g 0.25,
# (4) pK, a number for PC neiborhood size e.g 0.09, or NULL to use estimation
# (5) npcs, an integer of the number of PCs used e.g 30,
# (6) annotation, a name metadata colname used for homotypic estimation.
# If annotation is not provided, we used Leiden clusters (reso=0.2) as default. 
# It returns the Seurat object with doublets classification in the medatadata.
# This function take a single Seutat object that contains one sample, 
# if an Assay5 object is used, the layers will be joined together.
run_doubletFinder <-
  function(seu,
           expected_prop = 0.075,
           pN = 0.25,
           pK = NULL,
           npcs = 30,
           annotation = 'SCT_snn_res.0.2') {
    cat("Run doubletFinder without ground truth on a Seurat with SCTransform\n")
    if (class(seu[["RNA"]]) == "Assay5") {
      if (length(seu@assays[["RNA"]]@layers) > 1) {
        seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
        cat("Joining layers in RNA assay...\n") 
        cat("Consider spliting the object if it contains more than one sample\n")
      }
    }
    if (is.null(annotation)) {
    seu_temp <- SCTransform(seu) %>%
      FindNeighbors(dims = 1:npcs) %>%
      RunUMAP(dims = 1:30, n.epochs = 1000) %>%
      FindClusters(resolution = 0.2, algorithm = 4)
    seu@meta.data[["Leiden0.2"]] <- seu_temp@meta.data[["SCT_snn_res.0.2"]]
    rm(seu_temp)
    annotation <- "Leiden0.2"
    }
    rev.db <- seu@meta.data[grep("^DF.classification", colnames(seu@meta.data))]
    unwanted_metadata <- colnames(rev.db)
    for (meta_col in unwanted_metadata) {
      if (meta_col %in% colnames(seu@meta.data)) {
        seu@meta.data[[meta_col]] <- NULL
      }
    }
    homotypic_prop <- modelHomotypic(seu@meta.data[[annotation]])
    homotypic_prop <- min(1, max(0, homotypic_prop))
    nExp_poi <- round(expected_prop * nrow(seu[[]]))
    cat("run_doubletFinder: Use nExp_poi:", nExp_poi, "\n")
    nExp_poi_adj <- max(0, round(nExp_poi * (1 - homotypic_prop)))
    cat("run_doubletFinder: Use nExp_poi_adj:", nExp_poi_adj, "\n")
    if (is.null(pK)) {
      if (('paramSweep_v3' %in% ls('package:DoubletFinder')) & 
          (packageDescription('Seurat')$Version == '4.1.3')) {
        sweep.res.list <- paramSweep_v3(seu, PCs = 1:npcs, sct = T)
      } else {
        sweep.res.list <- paramSweep(seu, PCs = 1:npcs, sct = T)
      }
      sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
      bcmvn <- find.pK(sweep.stats)
      mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
      print(mpK)
      pK_set  <- mpK
      cat("run_doubletFinder: Use mpK:",mpK,"\n")
      } 
    else {
      pK_set  <- pK
      cat("run_doubletFinder: Use pK:",  pK,"\n")
      }
    #run doubletFinder
    for (nExp in c(nExp_poi, nExp_poi_adj)) {
      if (('paramSweep_v3' %in% ls('package:DoubletFinder')) & 
          (packageDescription('Seurat')$Version == '4.1.3')) {
        seu <- doubletFinder_v3(seu, 
                             PCs = 1:npcs, pN = pN, pK = pK_set, nExp = nExp, 
                             reuse.pANN = FALSE, sct = TRUE)
      } else {
        seu <- doubletFinder(seu, 
                             PCs = 1:npcs, pN = pN, pK = pK_set, nExp = nExp, 
                             reuse.pANN = FALSE, sct = TRUE)
      }
    }
    rep.db <- seu@meta.data[grep("^DF.classification", colnames(seu@meta.data))]
    seu <- AddMetaData(seu, metadata = rep.db[1], col.name = 'dbFinder_poi')
    seu <- AddMetaData(seu, metadata = rep.db[2], col.name = 'dbFinder_adj')
    plt <- DimPlot(seu, group.by = 'dbFinder_adj') + theme_bw() + NoLegend()
    print(plt)
    return(seu)
  }
# ens_sm <- run_doubletFinder(ens_sm, pK = NULL)
################################################################################
# Utility function for the analysis of a Seurat object based on SCT method.
# It takes a Seurat object and returns analyzed object.
# This function takes 11 arguments:
# (1)  A Seurat Object,
# (2)  remove_cell: a vector of cell names to be removed,
# (3)  mask_gene: a vector containing genes to be masked from high variable genes.
# (4)  cluster_reso: passed as "resolution" argument to Seurat::FindClusters()
# (5)  cluster_algo: passed as "algorithm" argument to Seurat::FindClusters()
# (6)  umap_dims: passed as "dims" argument to Seurat::RunUMAP()
# (7)  umap_n_neighbhors: passed as "n.neighbors" argument to Seurat::RunUMAP()
# (8)  umap_min_dist: passed as "min.dist" argument to Seurat::RunUMAP()
# (9)  umap_n_epochs: passed as "n.epochs" argument to Seurat::RunUMAP()
# (10) umap_local_connectivity: passed as "local.connectivity" argument to 
#      Seurat::RunUMAP()
# (11) neighbors_dims: passed as "dims" argument to Seurat::FindClusters()
analyze_sctseurat <- function(seu,
                              mask_gene    = NULL,
                              remove_cell  = NULL,
                              cluster_reso = 1.8,
                              cluster_algo = 4,
                              umap_dims    = 1:20,
                              umap_umap_method = "uwot",
                              umap_metric = "cosine",
                              umap_n_neighbhors = 55L,
                              umap_min_dist = 0.35,
                              umap_n_epochs = 2000,
                              umap_local_connectivity = 1L,
                              neighbors_dims = 1:20) {
  if (!is.null(remove_cell)) {
    seu <- subset(seu, cells = remove_cell, invert = T)
  }
  seu <- SCTransform(seu)
  sex_gene <- c("Xist","Gm13305","Tsix","Gm8730","Eif2s3y","Ddx3y","Uty","Kdm5d")
  ime_gene <- c("Fos", "Jun", "Junb")
  if (is.null(mask_gene)) {
    vg <- VariableFeatures(seu)
  } 
  else if (mask_gene == "sex") {
    vg <- VariableFeatures(seu)[!VariableFeatures(seu) %in% sex_gene]
  }
  else if (mask_gene == "stress") {
    vg <- VariableFeatures(seu)[!VariableFeatures(seu) %in% ime_gene]
  }
  else if (mask_gene == "both") {
    vg <- VariableFeatures(seu)[!VariableFeatures(seu) %in% c(ime_gene,sex_gene)]
  }
  else {
    vg <- VariableFeatures(seu)[!VariableFeatures(seu) %in% mask_gene]
  }
  seu <- RunPCA(seu, features = vg)
  seu <- RunUMAP(seu, 
                 dims = umap_dims, 
                 n.neighbors = umap_n_neighbhors, 
                 min.dist = umap_min_dist, 
                 n.epochs = umap_n_epochs, 
                 local.connectivity = umap_local_connectivity) 
  seu <- FindNeighbors(seu, dims = neighbors_dims)
  seu <- FindClusters(seu, resolution = cluster_reso, algorithm = cluster_algo)
  plt <- DimPlot(seu, reduction = 'umap', label = T, repel = T) + 
    theme_bw() + NoLegend()
  print(plt)
  return(seu)
}
# ens_sm <- analyze_sctseurat(ens_sm)
################################################################################
# Utility function to plot (a dot plot) top two markers of an analyzed object. 
# This function takes 4 arguments:
# (1) A Seurat object,
# (2) marker: output of Seurat::FindAllMarkers(),
# (3) cols: to be passed as "cols" argument to Seurat::DotPlot(),
# (4) size_scale: to be passed as "dot.scale" argument to to Seurat::DotPlot().
plot_top2_marker <- function(seu, 
                             marker = NULL, 
                             cols="RdYlBu", 
                             size_scale=10,
                             gene_angle = 60){
  if (is.null(marker)) {
    if (length(Layers(seu)) > 1) {
      seu_tmp <- seu
      seu_tmp <- PrepSCTFindMarkers(seu)
      marker  <- FindAllMarkers(seu)
      rm(seu_tmp)
    }
    marker <- FindAllMarkers(seu)
  }
  top2_marker <- marker %>% group_by(cluster) %>% top_n(n=2, wt=avg_log2FC)
  dp <- DotPlot(seu, features=unique(top2_marker$gene), cols="RdYlBu", 
                dot.scale=size_scale, scale.by="size")
  plt <- dp + geom_point(aes(size=pct.exp), shape=21, colour="black", stroke=0.25) + 
    theme_bw(base_line_size = 0.25) + theme(axis.title.x=element_blank(), 
                       axis.title.y=element_blank(),
                       axis.text.x = element_text(size = 7, angle = gene_angle, 
                                                  hjust = 1))
  # + NoLegend() 
  print(plt)
  # + guides(size=guide_legend(
  # override.aes=list(shape=21, colour="black", fill="white"))) 
}
################################################################################
# Utility function for simple analysis of a list of Seurat objects (SCT method).
# Used for expolatory purposes.
# It takes a list of Seurat objects and returns a list of analyzed objects.
run_simple_sctseurat_list <- function(seurat_list) {
  analyzed_list <- list()
  for (i in seq_along(seurat_list)) {
    seu  <- seurat_list[[i]]
    seu  <- SCTransform(seu) %>%
      RunPCA() %>%
      FindNeighbors(dims = 1:30) %>%
      RunUMAP(dims = 1:50, n.epochs = 1000) %>%
      FindClusters()
    plot <- DimPlot(seu, reduction = 'umap') + theme_bw() + NoLegend()
    print(plot)
    analyzed_list[[i]] <- seu
  }
  return(analyzed_list)
}
################################################################################
# Utility function for simple clustering a Seurat object (SCT method).
# Use for expolatory purposes.
# It takes a Seurat object and returns the analyzed object.
run_simple_sctseurat <- function(seu) {
  seu  <- SCTransform(seu) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:50, n.epochs = 1000) %>%
    FindClusters()
  plot <- DimPlot(seu, reduction = 'umap') + theme_bw() + NoLegend()
  print(plot)
  return(seu)
}
# ens_sm <- run_simple_sctseurat(ens_sm)
################################################################################
# This function warps pre-processing steps of scRNAseq data combining procedures 
# and parameters described in the manuscript. It has been updated to be 
# compatible with Seurat version 5.1.0, as well as doubletFinder 2.0.4 (orginal 
# analyses was performed with Seurat version 4.1.3, and with paramSweep_v3()/ 
# doubletFinder_v3(), which will also work if these versions are loaded).
# This function takes 5 arguments:
# (1) A list of seurat objects,
# (2) dbfinder: TRUE, if doubletFunder is to be used, 
# (3) dbfinder_pK: passed as "pK" argument to doubletFinder::doubletFinder(),
# (4) dbfinder_annotation: a name metadata colname used for homotypic estimation,
# (5) dbfinder_expected_prop: a decimal for expected doublets (loading density).
pre_process_ens <- function(seurat_list,
                            dbfinder = TRUE,
                            dbfinder_annotation = NULL,
                            dbfinder_expected_prop = 0.075,
                            dbfinder_pK = NULL
                            ) {
  cat("Data preparation according to Li and Morarach et al. 2024\n")
  cat("Starting preprocessing...\n")
  if (!is.list(seurat_list)) {
    seurat_list <- list(seurat_list)
  }
  if (is.null(names(seurat_list))) {
    names(seurat_list) <- paste0("seurat_", seq_along(seurat_list))
  }
  processed_list <- lapply(seurat_list, function(seu) {
    cat("\nProcessing object:", seu@meta.data$sample, "\n")
    report_seurat_status(seu)
    seu_raw <- seu
    cat("Filtering cells based on thresholds...\n")
    cell_to_keep <- filter_cell_by_threshold(seu)
    cat("Removing", length(setdiff(cell_to_keep, colnames(seu))), "cells...","\n")
    seu <- subset(seu, cells = cell_to_keep, invert = FALSE)
    cat("Running preliminary processing (SCTransform, PCA, UMAP, Neighbors)...\n")
    seu <- SCTransform(seu) %>%
      RunPCA() %>% 
      RunUMAP(dims = 1:30, n.epochs = 1000) %>% 
      FindNeighbors(dims = 1:30) %>%
      FindClusters(resolution = 0.2, algorithm = 4)
    plt_umap <- DimPlot(seu, reduction = 'umap') + theme_bw() + NoLegend()
    print(plt_umap)
    seu$Leiden0_2 <- seu$SCT_snn_res.0.2
    cat("Performing a resolution sweep clustering...\n")
    seu <- get_reso_sweep_cluster(seu)
    cat("Filtering plausible cluster outliers...\n")
    outlier <- get_cluster_outlier(seu)
    cat("Removing:", length(outlier), "cells...", "\n")
    seu <- subset(seu, cells = outlier, invert = TRUE)
    cell_remain <- colnames(seu)
    if (dbfinder == FALSE) {
      seu_raw <- subset(seu_raw, cells = cell_remain)
      cat("Global thresholding passed:", length(cell_to_keep), "cells", "\n")
      cat("Outliers removed:", length(outlier),     "cells", "\n")
      cat("Cells remain    :", length(cell_remain), "cells", "\n")
    }
    if (dbfinder == TRUE)  {
      cat("Running DoubletFinder...\n")
      seu <- run_doubletFinder(seu, pK = dbfinder_pK, 
                               annotation = dbfinder_annotation,
                               expected_prop = dbfinder_expected_prop)
      singlet <- WhichCells(seu, expression = dbFinder_adj == 'Singlet')
      doublet <- WhichCells(seu, expression = dbFinder_adj == 'Doublet')
      seu_raw <- subset(seu_raw, cells = singlet)
      seu_raw <- subset(seu_raw, cells = cell_to_keep)
      rm(seu)
      cat("Sample processed.\n")
      report_seurat_status(seu_raw)
      cat("Global thresholding passed:", length(cell_to_keep), "cells", "\n")
      cat("Outliers removed  :", length(outlier), "cells", "\n")
      cat("Doublets removed  :", length(doublet), "cells", "\n")
      cat("Singlets remaining:", length(singlet), "cells", "\n")
    }
    return(seu_raw)
    cat("==========================================================\n\n")
  })
  names(processed_list) <- names(seurat_list)
  cat("Preprocessing completed.\n")
  return(processed_list)
}
# ens_sm  <- pre_process_ens(list(ens_sm))
# load_sm <- pre_process_ens(load_sm)
################################################################################
# Utility function for rawcount visualization as in the Seurat paper 
# (Choudhary and Satija, 2022: https://doi.org/10.1186/s13059-021-02584-9)
plot_rawcount <- function(seu, seurat_object_version = "v3"){
  if (seurat_object_version == "v3") {
    rawdata <- seu@assays[["RNA"]]@counts
  } else if (seurat_object_version == "v5") {
    rawdata <- seu@assays[["RNA"]]@layers[["counts"]]
  } else {
    cat("Set seurat_object_version to either \"v3\" or \"v5\"\n")
  }
  
  gene_attr <- data.frame(mean = rowMeans(rawdata), 
                          detection_rate = rowMeans(rawdata > 0), 
                          var = apply(rawdata, 1, var))
  gene_attr$log_mean <- log10(gene_attr$mean)
  gene_attr$log_var  <- log10(gene_attr$var)
  rownames(gene_attr)<- rownames(rawdata)
  cell_attr <- data.frame(n_umi = colSums(rawdata), 
                          n_gene = colSums(rawdata > 0))
  rownames(cell_attr) <- colnames(rawdata)
  p1 <- ggplot(gene_attr, aes(log_mean, log_var)) + 
    geom_point(alpha = 0.3, shape = 16, size = 0.2) + 
    geom_density_2d(size = 0.3) + 
    geom_abline(intercept = 0, slope = 1, color = "red") + 
    theme_bw(base_line_size = 0.33) + NoLegend()
  x = seq(from = -3, to = 2, length.out = 1000)
  poisson_model <- data.frame(log_mean = x, 
                              detection_rate = 1 - dpois(0, lambda = 10 ^ x))
  p2 <- ggplot(gene_attr, aes(log_mean, detection_rate)) + 
    geom_point(alpha = 0.3, shape = 16, size = 0.3) + 
    geom_line(data = poisson_model, color = "red") + 
    theme_bw(base_line_size = 0.33) + NoLegend()
  p3 <- ggplot(cell_attr, aes(n_umi, n_gene)) +
    geom_point(alpha = 0.3, size = 0.3, shape = 16) + 
    geom_density_2d(size = 0.33) + 
    theme_bw(base_line_size = 0.33) + NoLegend()
  p4 <- p1 | p2 | p3
  print(p4)
}
################################################################################
# Utility function for gene visualization 
plot_violin <- function(ens, 
                        gene_plot, 
                        col = NULL,
                        group_cell = NULL,
                        x_title = "SCTNormalized Expression", 
                        y_title = "Seurat Ident Clusters",
                        gene_angle = 60){
  require(cols4all)
  # install.packages("remotes")
  # remotes::install_github("mtennekes/cols4all", dependencies = TRUE)
  c2acol <- list(c4a('20',15), c4a('kelly', 22))
  vln <- VlnPlot(ens, features = gene_plot, add.noise = T, fill.by = 'feature',
                 group.by = group_cell, stack = T, alpha = 0.8) 
  v1 <- vln + geom_violin(scale = 'width', 
                          draw_quantiles = c(0.25,0.5,0.75), 
                          color = 'black', size = 0.25, 
                          alpha = 0.25) 
  if (is.null(col)) {
    p <- v1
  }
  else if (col == 'c2a1') {
    p <- v1 +  scale_fill_manual(values = c2acol[[1]])
  }
  else if (col == 'c2a2') {
    p <- v1 +  scale_fill_manual(values = c2acol[[2]])
  }
  else {
    p <- v1 +  scale_fill_manual(values = col)
  }
  plt <- p + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text.x =  element_text(size = 8), 
          axis.text.y =  element_text(size = 8),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          strip.text.x = element_text(size = 8, angle = gene_angle),
          strip.background = element_blank()
    ) + 
    labs(x = x_title)+
    labs(y = y_title) + NoLegend() 
  print(plt)
}
################################################################################
# Utility function for resolution-sweep cluster visualization 
plot_clusttree <- function(seu, reso_min  =  0.1,  reso_max    = 2.0, 
                           reso_increment =  0.1,  seurat_algo = 1,
                           clustree_prop_filter  = 0.2,  
                           clustree_count_filter = 0) {
  sweep_reso <- seq(reso_min, reso_max, by = reso_increment)
  seu_temp   <- seu
  cat("Perform Seurat's clustering, resolution:", sweep_reso,"...\n")
  rev.clustree <- seu_temp@meta.data[grep("^SCT_snn_res.", 
                                          colnames(seu_temp@meta.data))]
  unwanted_metadata <- colnames(rev.clustree)
  for (meta_col in unwanted_metadata) {
    if (meta_col %in% colnames(seu_temp@meta.data)) {
      seu_temp@meta.data[[meta_col]] <- NULL
    }
  }
  for (i in seq_along(sweep_reso)) {
    seu_temp <- FindClusters(seu_temp, resolution = sweep_reso[i], 
                             algorithm = seurat_algo)
  }
  cat("Creating Crustree Plot...\n")
  p1 <- clustree(seu_temp, node_colour = 'grey90',
                 node_text_size = 3, show_axis = T, node_size_range = c(0.5, 6),
                 edge_width = 0.5, prop_filter = clustree_prop_filter,
                 count_filter = clustree_count_filter)
  p2 <- p1 + theme_bw(base_line_size = 0.2, base_size = 8) +
    theme(axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +  
    labs(y = "Cluster Resolutions")
  print(p2)
  rm(seu_temp)
  cat("Complete .\n")
}
# plot_clusttree(ens_sm_lv3, seurat_algo = 4)
################################################################################
#                                    ANALYSES                                  #
################################################################################
# Section 1. Updated Clustering Analysis (for Seurat5 framework)
# 
# An illustrative example of our clustering analysis with some updates to be
# compatible with Seurat5 framework. Here we traeted each library as a separate 
# batch stored in different layers of an Assay5 object. Likewise, doubletFinder
# was applied to each library with pK separately estimated. For our original 
# analysis, please go to section 2 below.
# 1.1) Load dataset.
folders <- list(
  "~/sm022/filtered_feature_bc_matrix",
  "~/sm023/filtered_feature_bc_matrix",
  "~/sm024/filtered_feature_bc_matrix"
)
load_sm <- load_data_10x(folders, seurat_object_version = "v5")
# 1.2) Add metada information for each library.
for (i in seq_along(load_sm)){
  load_sm[[i]]$pre_sorted <- "Baf53bCre-R26Rtomato"
  load_sm[[i]]$region <- "small_intestine"
  load_sm[[i]]$seq <- "10x-v3"
  load_sm[[i]]$age <- "p21"
  load_sm[[i]]$layer <- "submucosal_plexus"
  load_sm[[i]]$orig_barcode <- colnames(load_sm[[i]])
  plot_rawcount(load_sm[[i]], seurat_object_version = "v5")
}
# 1.3) Pre-process data, treating each library as a distinct batch.
load_sm <- pre_process_ens(load_sm, dbfinder_expected_prop = 0.075)
# 1.4) Merge data accorging to Seurat5 framework(one library/layer).
ens_sm <- merge(
  load_sm[[1]],
  y = c(load_sm[[2]], load_sm[[3]]),
  add.cell.ids = c("sm022_", "sm023_", "sm024_"),
  project = "ENSP24"
)
# Save QC data in the working directory
saveRDS(ens_sm, file = "~/ens_sm.rds")
# 1.5) Perform iterative analysis to remove non-enteric neurons.
# Level 1 Clustering:
ens_sm_lv1 <- analyze_sctseurat(ens_sm, mask_gene = 'sex')
ens_sm_lv1 <- PrepSCTFindMarkers(ens_sm_lv1)
ens_sm_lv1_markers <- FindAllMarkers(ens_sm_lv1)
# Save marker table in the working directory
write.table(ens_sm_lv1_markers, file = "ens_sm_lv1_markers.csv", sep = ",", 
            col.names = NA, qmethod = "double")
# Vizualization of some markers
plot_top2_marker(ens_sm_lv1, marker = ens_sm_lv1_markers) + NoLegend()
gene_plot = c('Elavl4','Tubb3','Sox10','Plp1','Cdh1','Tph1','Lgr5','Pyy','H2-Aa',
              'Pecam1','Myh11','Acta2','Pdgfra','Col1a1','Dcn','Lyve1')
plot_violin(ens_sm_lv1, gene_plot = gene_plot, col = 'c2a2')
# Define non-enteric neuron clusters to be removed in the next level
cell_removed_lv1 <- WhichCells(ens_sm_lv1, 
                               idents = c(21, 24, 27, #epithelial cells
                                          30, #immune cells
                                          28, #fibroblasts, smooth muscle, 
                                              #endothelial cells
                                          31  #lymphatic epithelial cells
                                          ))
# Level 2 Clustering:
ens_sm_lv2 <- analyze_sctseurat(ens_sm,
                                mask_gene = 'sex',
                                remove_cell = cell_removed_lv1)
ens_sm_lv2 <- PrepSCTFindMarkers(ens_sm_lv2)
ens_sm_lv2_markers <- FindAllMarkers(ens_sm_lv2)
write.table(ens_sm_lv2_markers, file = "ens_sm_lv2_markers.csv", sep = ",", 
            col.names = NA, qmethod = "double")
plot_top2_marker(ens_sm_lv2, ens_sm_lv2_markers) + NoLegend()
gene_plot = c('Elavl4','Tubb3','Gfap','Sox10','Plp1','Penk','Nos1','Ndufa4l2',
              'Cck','Ntng1','Tac1','Calca','Nmu','Cbln2')
plot_violin(ens_sm_lv2, gene_plot = gene_plot, col = 'c2a2')
cell_removed_lv2 <- WhichCells(ens_sm_lv2, 
                               idents = c(21, 25,     #glia or glia-neurons
                                          22,23,24,28 #myenteric neurons
                               ))
# Level 3 Clustering:
ens_sm_lv3 <- analyze_sctseurat(ens_sm, 
                                mask_gene = 'sex', 
                                cluster_reso = 0.8,
                                remove_cell = c(cell_removed_lv2, 
                                                cell_removed_lv1),
                                neighbors_dims = 1:50
                                )
ens_sm_lv3 <- PrepSCTFindMarkers(ens_sm_lv3)
ens_sm_lv3_markers <- FindAllMarkers(ens_sm_lv3)
write.table(ens_sm_lv3_markers, file = "ens_sm_lv3_markers.csv", sep = ",", 
            col.names = NA, qmethod = "double")
plot_top2_marker(ens_sm_lv3, ens_sm_lv3_markers) + NoLegend()
gene_plot = c('nCount_RNA','nFeature_RNA','percent.mt','Cbln2','Adgrg6','Vip',
              'Npy','Th','Dbh','Cd24a','Sst','Pcdh7','Dmkn', 'Paip2b')
plot_violin(ens_sm_lv3, gene_plot = gene_plot, col = 'c2a2', gene_angle = 90)
# At this level, clusters should contains only enteric neurons.
ens_sm_lv3 <- FindClusters(ens_sm_lv3, algorithm = 4, resolution = 0.4)
DimPlot(ens_sm_lv3, group.by = 'SCT_snn_res.0.4') + theme_bw() + NoLegend()
# Save analyzed data in the working directory
saveRDS(ens_sm_lv3, file = "~/ens_sm_lv3.rds")
################################################################################
################################################################################
# Section 2. Reproducing Clustering Analysis for P24 dataset
# 
# Setting seurat_object_version = "v3" in  load_data_10x() would effectively
# create an Assay (v3) Seurat object which were used to produced our results, 
# For this we merged sm022, sm023 and sm023 to a single Seurat object and 
# computed an SCT model for this merged object. Similary pK for estimated for 
# the merged object before running doubletFinder, for which we used pK = 0.09 
# and expected_prop = 0.075. Please note that we prodeced our results with
# Seurat version 4.1.3, and with paramSweep_v3() and doubletFinder_v3(). Still,
# we found no significant disparity in outcomes when used current versions (
# Seurat version 5.1.0, doubletFinder 2.0.4: 20/12/24).
folders <- list(
  "~/sm022/filtered_feature_bc_matrix",
  "~/sm023/filtered_feature_bc_matrix",
  "~/sm024/filtered_feature_bc_matrix"
)
load_sm <- load_data_10x(folders, seurat_object_version = "v3")
for (i in seq_along(load_sm)){
  load_sm[[i]]$pre_sorted <- "Baf53bCre-R26Rtomato"
  load_sm[[i]]$region <- "small_intestine"
  load_sm[[i]]$seq <- "10x-v3"
  load_sm[[i]]$age <- "p21"
  load_sm[[i]]$layer <- "submucosal_plexus"
  load_sm[[i]]$orig_barcode <- colnames(load_sm[[i]])
  plot_rawcount(load_sm[[i]], seurat_object_version = "v3")
}
ens_sm<- merge(
  load_sm[[1]],
  y = c(load_sm[[2]], load_sm[[3]]),
  add.cell.ids = c("sm022_", "sm023_", "sm024_"),
  project = "ENSP24"
)
# Alternatively data layers (in Assay5) can be joined together:
# ens_sm[["RNA"]] <- JoinLayers(ens_sm[["RNA"]])
ens_sm <- pre_process_ens(list(ens_sm),
                          dbfinder_pK = 0.09,
                          dbfinder_annotation = "SCT_snn_res.0.2",
                          dbfinder_expected_prop = 0.075)
# alternatively one can extract cells passing the QC from the metadata table 
# mdata  <- read.csv("~/NNReviewRound1/sm_metadata.csv", row.names=1)
# ens_sm <- subset(ens_sm, cells = rownames(mdata)) #, Or
# ens_sm <- readRDS("~/sm_p24_lv3.rds") for reanalyzis of processed Seurat object
# Level 1 Clustering
ens_sm_lv1 <- analyze_sctseurat(ens_sm[[1]], 
                                mask_gene = 'sex',
                                cluster_reso = 1.8,
                                cluster_algo = 1,
                                neighbors_dims = 1:50
                                )
# Define non-enteric neuron clusters, based on marker expression.
cell_removed_lv1 <- WhichCells(ens_sm_lv1, 
                               idents = c('24', '27', '17', '22', 
                                          '18', '21', '23')
                               )
# Level 2 Clustering
ens_sm_lv2 <- analyze_sctseurat(ens_sm_lv1, 
                                mask_gene = 'sex',
                                remove_cell = cell_removed_lv1,
                                cluster_reso = 0.8,
                                cluster_algo = 4,
                                neighbors_dims = 1:50
                                )
# Define remaining non-enteric neuron clusters, based on markers.
cell_removed_lv2 <- WhichCells(ens_sm_lv2, 
                               idents = c('13','17','14','19',
                                          '12','18','15','16')
                               )
# Level 3 Clustering
ens_sm_lv3 <- analyze_sctseurat(ens_sm_lv2, 
                                mask_gene = 'sex',
                                remove_cell = cell_removed_lv2,
                                cluster_reso = 0.8,
                                cluster_algo = 4,
                                neighbors_dims = 1:50
                                )
# Join clusters to cell classes based on inspection of marker expression.
# Combine clusters to cell classes, rename to smENC1, smENC2 and smENC3, and add
# to the object metadata as 'cell_class'. 
ens_sm_lv3 <- RenameIdents(
  ens_sm_lv3,
  "8" = "smENC1",
  "5" = "smENC2",
  "4" = "smENC2",
  "3" = "smENC2",
  "6" = "smENC2",
  "2" = "smENC3",
  "1" = "smENC3",
  "8" = "smENC3",
  "7" = "smENC3")
ens_sm_lv3@meta.data$cell.class <- ens_sm_lv3@active.ident
DimPlot(ens_sm_lv3) + theme_bw() + NoLegend()
# Save analyzed data in the working directory
saveRDS(ens_sm_lv3, file = "~/ens_sm_lv3.rds")
# UMAP plot in Figure 1 was performed with the following settings:
# RunUMAP(ens_sm_lv3, seed.use = 1234, dims = 1:20, 
# umap.method = "umap-learn", metric = "cosine", n.neighbors = 45L, 
# min.dist = 0.35, spread = 1, n.epochs = 2000, local.connectivity = 2L)
################################################################################
################################################################################
# Section 3. Reproducing Clustering Analysis for P7 dataset
# Seurat version 4.1.3
folders <- list(
"~/201/filtered_feature_bc_matrix",
"~/202/filtered_feature_bc_matrix"
)
load_p7 <- load_data_10x(folders, seurat_object_version = "v3")
for (i in seq_along(load_p7)){
  load_p7[[i]]$pre_sorted <- "Wint1Cre-R26Rtomato"
  load_p7[[i]]$region <- "small_intestine"
  load_p7[[i]]$seq <- "10x-v2"
  load_p7[[i]]$age <- "p7"
  load_p7[[i]]$layer <- "all"
  load_p7[[i]]$orig_barcode <- colnames(load_p7[[i]])
  plot_rawcount(load_p7[[i]], seurat_object_version = "v3")
}
ens_p7<- merge(
  load_p7[[1]],
  y = c(load_p7[[2]]),
  add.cell.ids = c("p7_201_", "p7_202_"),
  project = "ENSP7"
)
# Filter cells by global thresholds
cell_to_keep <- filter_cell_by_threshold(ens_p7,
                                         nfeature_bounds = c(500, 6500),
                                         ncount_bounds = c(0, 40000),
                                         percent_mt_threshold = 10)
ens_p7 <- subset(ens_p7, cells = cell_to_keep, invert = FALSE)
# Perform analysis with sex and immediate-early genes masked when compute PCA
# Level 1 clustering
ens_p7_lv1 <- analyze_sctseurat(ens_p7,
                                mask_gene = "both", #sex, immediate-early genes
                                cluster_reso = 0.8,
                                cluster_algo = 1,
                                neighbors_dims = 1:50)
# Identify non-ENS clusters to be removed in the next level
cell_removed_p7_lv1 <- WhichCells(ens_p7_lv1,
                                  idents = c('12','18','19','20','21'),
                                  invert = F)
# Level 2 clustering
# Perform re-analysis with non-ENS clusters removed
ens_p7_lv2 <- analyze_sctseurat(ens_p7_lv1,
                                mask_gene = "both", #sex, immediate-early genes
                                remove_cell = cell_removed_p7_lv1,
                                cluster_reso = 1.8,
                                cluster_algo = 1,
                                neighbors_dims = 1:50)
# UMAP plot in Figure 6 was performed with the following settings:
# RunUMAP(ens_p7_lv1, dims = 1:50, n.neighbors = 60L, min.dist = 0.40, 
# n.epochs = 10000, local.connectivity = 2L)
################################################################################
################################################################################
# Section 4. Reproducing label transfaer 
# Seurat version 4.1.3
# For Seurat 5.1.0, setting seurat_object_version = "v3" in  load_data_10x() 
# would effectively create an Assay (v3) Seurat object which were used to 
# produced our results. For this we merged p7_201 and p7_202 to a single Seurat 
# object and computed an SCT model for this merged object.
# ens_reference_map is a merged object of submucosal dataset (current study), 
# myenteric dataset (Morarach et al, 2022) and enteric glia (Zeisel et al. 2018),
# L2 downsampled 1000 cells. These are provided in "ens_reference_datasets.rds".
ens_refmap <- readRDS(file = "~/ens_reference_map.RDS")
ens_query  <- readRDS(file = "~/ens_p7_lv2.RDS")
ens_refmap@active.ident <- ens_refmap$published.type
ens_refmap <- subset(ens_refmap, idents = c("ENC11", "ENC5", "ENC6"))
ens_refmap <- SCTransform(ens_refmap, method = "glmGamPoi") %>%
  RunUMAP(dims = 1:50, n.epochs = 1000)
ens_anchors <- FindTransferAnchors(
  reference = ens_refmap,
  query = ens_query,
  dims = 1:50,
  reference.reduction = "pca",
  query.assay = "SCT",
  reference.assay = "SCT",
  normalization.method = "SCT")
prediction <- TransferData(anchorset = ens_anchors, refdata = 50)
ens_query <- AddMetaData(ens_query, metadata = prediction)
################################################################################
################################################################################
