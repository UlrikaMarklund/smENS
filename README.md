# Enteric Neuron Subtype Analysis Workflow

This repository provides a workflow for analyzing enteric neuron subtypes using single-cell RNA sequencing (scRNA-seq) data processed with the **Seurat** framework. It accompanies the manuscript *Discovering the transcriptomes, connections, and development of submucosal neuron classes in the mouse small intestine* by Li and Morarch et al., 2024. The associated data is deposited in the Gene Expression Omnibus (GEO) database under the identifier **GSE263422**.

The workflow includes key steps such as data loading, preprocessing, clustering, and iterative removal of non-enteric cells. The analysis procedures described in the manuscript are wrapped into a series of custom functions with adjustable parameters for reproducibility and ease of use. The main script, `RSCRIPT.R`, is organized into three parts:

- **Part 1: Libraries** – Load required R packages.  
- **Part 2: Functions** – Define custom functions for preprocessing and analysis.  
- **Part 3: Analyses** – Perform the analyses, organized into four main sections.
- (1.) General Workflow for Cluster Analysis (Compatible with Seurat v5)
- (2.) P21 Submucosal Plexus Cluster Analysis
- (3.) P7 Small Intestine (ENS) Cluster Analysis
- (4.) Label Transfer 

To follow this workflow, first load the Libraries and Functions sections of the `RSCRIPT.R` file. These only need to be loaded once per session. The Analyses section can then be executed as needed.

---
## Analysis Overview
![Workflow Diagram](workflow.png)
#### Main functions
This concise walk through is to illustrate additional steps when working with Seurat v5. *For detailes, refer to the comments and descriptions embedded in the script.*
1. **Data Loading**
   - **`load_10x()`**: Load single-cell RNA sequencing data from 10x Genomics.
2. **Pre-processing**
   - **`filter_cell_by_threshold()`**: Exclude low-quality cells based on thresholds for `nFeature_RNA` and `nCount_RNA`.
   - **`get_cluster_outlier()`**: Exclude extreme outliers, sweeping different resolution values.
   - **`run_dbfinder()`**: Exclude doublets by `DoubletFinder` (with Seurat’s SCT)
3. **Iterative Clustering and removing of non-relevant clusters**
   - **`analyze_sctseurat()`**: Identify and subset non-enteric neuron clusters iteratively.
---
#### 1. Data Loading
#### **1.1 Load Dataset**
Specify the folders containing paths to the CellRanger's output, `filtered_feature_bc_matrix`:
```R
folders <- list(
  "~/sm022/filtered_feature_bc_matrix",
  "~/sm023/filtered_feature_bc_matrix",
  "~/sm024/filtered_feature_bc_matrix"
)
load_sm <- load_data_10x(folders, seurat_object_version = "v5")
```
#### **1.2 Add Metadata Information**
```R
for (i in seq_along(load_sm)) {
  load_sm[[i]]$pre_sorted <- "Baf53bCre-R26Rtomato"
  load_sm[[i]]$region <- "small_intestine"
  load_sm[[i]]$seq <- "10x-v3"
  load_sm[[i]]$age <- "p21"
  load_sm[[i]]$layer <- "submucosal_plexus"
  load_sm[[i]]$orig_barcode <- colnames(load_sm[[i]])
  plot_rawcount(load_sm[[i]], seurat_object_version = "v5")
}
```
---
#### 2. Data Preprocessing
#### **2.1 Preprocess Data**
```R
load_sm <- pre_process_ens(load_sm, dbfinder_expected_prop = 0.075)
```
#### **2.2 Merge Datasets**
```R
ens_sm <- merge(load_sm[[1]], y = c(load_sm[[2]], load_sm[[3]]), add.cell.ids = c("sm022_", "sm023_", "sm024_"), project = "ENSP24")
```
---
#### 3. Iterative Clustering and Cluster Removal
#### **3.1 Level 1 Clustering**
Perform clustering to identify non-enteric neurons:
```R
ens_sm_lv1 <- analyze_sctseurat(ens_sm, mask_gene = 'sex')
ens_sm_lv1 <- PrepSCTFindMarkers(ens_sm_lv1)
ens_sm_lv1_markers <- FindAllMarkers(ens_sm_lv1)
```
Visualize markers:
```R
plot_top2_marker(ens_sm_lv1, marker = ens_sm_lv1_markers) + NoLegend()
```
Identify non-enteric clusters:
```R
cell_removed_lv1 <- WhichCells(ens_sm_lv1, idents = c(21, 24, 27, 30, 28, 31))
```
#### **3.2 Level 2 Clustering**
Refine the analysis:
```R
ens_sm_lv2 <- analyze_sctseurat(ens_sm, mask_gene = 'sex', remove_cell = cell_removed_lv1)
ens_sm_lv2 <- PrepSCTFindMarkers(ens_sm_lv2)
ens_sm_lv2_markers <- FindAllMarkers(ens_sm_lv2)
```
Further identify non-enteric clusters:
```R
cell_removed_lv2 <- WhichCells(ens_sm_lv2, idents = c(21, 25, 22, 23, 24, 28))
```
#### **3.3 Level 3 Clustering**
```R
ens_sm_lv3 <- analyze_sctseurat(ens_sm, mask_gene = 'sex', cluster_reso = 0.8, remove_cell = c(cell_removed_lv2, cell_removed_lv1), neighbors_dims = 1:50)
ens_sm_lv3 <- PrepSCTFindMarkers(ens_sm_lv3)
ens_sm_lv3_markers <- FindAllMarkers(ens_sm_lv3)
```
---
#### 4. Cluster Annotation
#### **4.1 Define Cell Classes**
Combine similar clusters into cell classes:
```R
ens_sm_lv3 <- RenameIdents(
  ens_sm_lv3,
  "8" = "smENC1",
  "5" = "smENC2",
  "4" = "smENC2",
  "3" = "smENC2",
  "6" = "smENC2",
  "2" = "smENC3",
  "1" = "smENC3",
  "7" = "smENC3"
)
ens_sm_lv3@meta.data$cell_class <- ens_sm_lv3@active.ident
DimPlot(ens_sm_lv3) + theme_bw() + NoLegend()
```
---
## Library Requirements
- **R 4.4.1**
  
This script works with:
- **cols4all:** v0.7-1
- **reticulate:** v1.39.0
- **patchwork:** v1.3.0
- **clustree:** v0.5.1
- **ggraph:** v2.2.1
- **ggplot2:** v3.5.1
- **DoubletFinder:** v2.0.4
- **dplyr:** v1.1.4
- **Seurat:** v5.1.0
- **SeuratObject:** v5.0.2
- **sp:** v2.1-4
- **stringr:** v1.5.1
### Reproducibility Requirements
However, to reproduce the results, specific package versions are required:
- **Seurat:** v4.1.3
- **DoubletFinder:** v2.0.3 (for functions `paramSweep_v3()` and `doubletFinder_v3()`)
---
## Citation
```
Discovering the transcriptomes, connections, and development of submucosal neuron classes in the mouse small intestine.
Wei Li1,#, Khomgrit Morarach1,#, Ziwei Liu1, Yanan Chen1, Sanghita Banerjee1, Fernando López-Redondo1,2, Ashley L. Harb3,4, Elham Jalalvand1, Anastassia Mikhailova1, David R. Linden3,4, Ulrika Marklund1,*.
```
---
