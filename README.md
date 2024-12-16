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
*For detailes, refer to the comments and descriptions embedded in the script.*
1. **Data Loading**
   - **`load_10x()`**: Load single-cell RNA sequencing data from 10x Genomics.
2. **Pre-processing**
   - **`filter_cell_by_threshold()`**: Exclude low-quality cells based on thresholds for `nFeature_RNA` and `nCount_RNA`.
   - **`get_cluster_outlier()`**: Exclude extreme outliers, sweeping different resolution values.
   - **`run_dbfinder()`**: Exclude doublets by `DoubletFinder` (with Seurat’s SCT)
3. **Iterative Clustering and removing of non-relevant clusters**
   - **`analyze_sctseurat()`**: Identify and subset non-enteric neuron clusters iteratively.
---
#### **Example of Dataset loading**
Specify the folders containing paths to the CellRanger's output, `filtered_feature_bc_matrix`:
```R
folders <- list(
  "~/sm022/filtered_feature_bc_matrix",
  "~/sm023/filtered_feature_bc_matrix",
  "~/sm024/filtered_feature_bc_matrix"
)
load_sm <- load_data_10x(folders, seurat_object_version = "v5")
```
---
## Library Requirements
- **R 4.4.1**
  
This script works with the following package versions:
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
To reproduce the results, specific package versions are required:
- **Seurat:** v4.1.3
- **DoubletFinder:** v2.0.3 (for functions `paramSweep_v3()` and `doubletFinder_v3()`)
---
## Citation
```
Discovering the transcriptomes, connections, and development of submucosal neuron classes in the mouse small intestine.
Wei Li1,#, Khomgrit Morarach1,#, Ziwei Liu1, Yanan Chen1, Sanghita Banerjee1, Fernando López-Redondo1,2, Ashley L. Harb3,4, Elham Jalalvand1, Anastassia Mikhailova1, David R. Linden3,4, Ulrika Marklund1,*.
```
---
