# Reproducible R Analysis Capsule
This capsule includes the R code and environment used for the manuscript, The BET Degrader ZBC260 Inhibits Breast Cancer Progression by Reprograming Macrophages in the Tumor Microenvironment, submitted to Cancer Research.

## Data Files
1.	Myeloid_cells_finals.rds – This file contains metadata generated from myeloid cell subclustering analysis. It was created after performing subclustering on the "myeloid cells" subset extracted from the total cells Seurat object ("obj"). The metadata includes detailed cluster information specific to myeloid cells and has been used for differential gene expression analysis within this subset. Additionally, this metadata will be utilized for transcription factor analysis focused on myeloid cell populations. While this metadata is suitable for studying cell-to-cell interactions within the myeloid subset, analyses involving interactions between myeloid cells and other cell types should use the total cells Seurat object, which includes all major cell identities such as tumor cells.

2.	Myeloid_cells_processed.rds – Description (Not needed)-please delete this metadata file.

## Scripts (Should be arranged sequence wise like total cells, myeloid subcluster) differential cells chat and TF (please check figure order)
1.	Cell_chat_analysis.R - This script analyzes interactions between myeloid cells and tumor cells. To perform this analysis accurately, both myeloid and tumor cells are first subsetted from the total cells Seurat object ("obj").
2.	Differential gene expression.R - This script performs differential gene expression analysis specifically on myeloid cells. It uses the myeloid cell metadata typically derived from the myeloid cell subset of the total cells Seurat object. The goal is to identify genes that are differentially expressed between distinct myeloid cell clusters or conditions, providing insights into functional differences within the myeloid population.
3.	Myeloid_cell_subcluster.R - This script performs subclustering of myeloid cells extracted from the total cells Seurat object “obj”. It refines the clustering of the myeloid cell subset to identify finer subpopulations within the major myeloid cell group. The output includes updated cluster assignments and metadata, which can be used for downstream analyses such as differential gene expression, transcription factor analysis, and cell-cell interaction studies focused on myeloid cells.
4.	T.F._analysis.R - This script performs transcription factor (TF) analysis specifically on myeloid cells. It examines how TF activity changes in response to ZBC260 treatment, using myeloid cell metadata and expression data to identify treatment-associated regulatory shifts within the myeloid cell population.
5.	Total_cell_population.R - This script processes and analyzes the total cell population using the complete Seurat object ("obj"). It includes identification and annotation of major cell types, including myeloid cells, tumor cells, and other immune or stromal populations.

## How to Run
Each analysis is meant to be run independently. Within each script it lists the appropriate packages that need to be installed for the analysis. The full list of packages for the project is listed below. 

## Required Packages

Core Analysis & Single-Cell Tools
Seurat – Comprehensive toolkit for single-cell RNA-seq data analysis and visualization.
SeuratObject – Underlying object class infrastructure for Seurat workflows.
SingleCellExperiment – Standardized S4 class for storing single-cell data with metadata.
SCENIC – Reconstructs gene regulatory networks and identifies active transcription factors in single cells.
AUCell – Calculates gene set activity in individual cells based on expression rankings.
RcisTarget – Identifies transcription factor binding motifs enriched in gene sets.
GENIE3 – Infers gene regulatory networks from expression data using tree-based methods.
CellChat – Infers, analyzes, and visualizes intercellular communication networks from scRNA-seq data.
scCATCH – Automatically annotates cell types in scRNA-seq data using known markers.
glmgampoi – Fast and scalable GLM-based method for analyzing overdispersed count data in single-cell experiments.
nebulosa – Enhances feature plots using kernel density estimation in single-cell datasets.

Visualization
ggplot2 – Grammar of graphics-based plotting system for creating complex visualizations.
ggpubr – Simplifies publication-ready plots using ggplot2 extensions.
ggrepel – Improves label readability by preventing text overlap in ggplot2 plots.
ggridges – Creates ridge plots (joy plots) for visualizing distributions.
ggalluvial – Creates alluvial/sankey diagrams to show categorical flow.
pheatmap – Generates customizable heatmaps with clustering options.
circlize – Creates circular plots, including chord diagrams and genomic plots.
patchwork – Combines multiple ggplot2 plots into a single layout.
networkD3 – Generates interactive network diagrams using D3.js.

Data Manipulation & General Analysis
tidyverse – Meta-package of core packages (dplyr, ggplot2, tidyr, etc.) for data science.
tidygraph – Provides a tidy API for graph data manipulation and analysis using the grammar of the tidyverse.
dplyr – Efficient data manipulation with verbs like filter, mutate, and summarize.
tibble – Lightweight, modern alternative to data frames with better printing and usability.
readr – Fast functions for reading rectangular text data (.csv, .tsv, etc.).
readxl – Reads Excel files (.xls, .xlsx) without external dependencies.
writexl – Writes data frames to Excel files (.xlsx) natively.
reshape2 – Reshapes data between wide and long formats.
Matrix – Provides dense and sparse matrix classes and methods.
stringr – Consistent functions for string manipulation.
cluster – Implements clustering algorithms and evaluation metrics.
parallel – Base R package for parallel computing support.
future – Asynchronous processing for parallel execution.
future.apply – Parallel versions of apply() functions using the future backend.
doParallel – Backend for the foreach package supporting parallel execution.
foreach – Provides looping constructs that support parallelism.
openxlsx – Read, write, and format Excel files without Java.
hiver – Efficient string comparison and matching using hierarchical index vectors.

Enrichment & Bioinformatics
clusterProfiler – Functional enrichment analysis for GO, KEGG, and more.
enrichR – Accesses multiple gene enrichment databases via the Enrichr API.
org.mm.eg.db – Genome-wide annotation package for mouse (Mus musculus).
dose – Disease ontology and semantic similarity analysis for enrichment tools.

Developer Tools
BiocManager – Installs and manages Bioconductor packages.
devtools – Streamlines R package development and GitHub installation.

## Notes
Please add any notes that may be needed. If nothing can be deleted. 
