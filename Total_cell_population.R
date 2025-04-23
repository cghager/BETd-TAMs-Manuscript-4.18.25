# 1. Load required libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(tidygraph)
library(tibble)
library(openxlsx)
library(SeuratObject)
library(devtools)
library(enrichR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(writexl)
library(tidyverse)
library(stringr)
library(BiocManager)
library(glmGamPoi)
library(cluster)
library(parallel)
library(DOSE)
library(Nebulosa)
library(future)
library(ggrepel)
library(ggpubr)
library(CellChat)

# 2. Set up file directory variables
tables_dir <- file.path(getwd(), "results", "tables")
figures_dir <- file.path(getwd(), "results", "figures")

# Ensure directories exist
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)



# 3. Load the samples (D2A1 samples, last 6 samples)
obj1 <- CreateSeuratObject(counts = Read10X(data.dir = "R:/hemonc-burness/Burness_Lab/Individual Data Folders/Deeksha Sharma/7.R DATA and Scripts/Burness Data/D2A1_CT1/11529-CH-7_sample_filtered_feature_bc_matrix"))
obj2 <- CreateSeuratObject(counts = Read10X(data.dir = "R:/hemonc-burness/Burness_Lab/Individual Data Folders/Deeksha Sharma/7.R DATA and Scripts/Burness Data/D2A1_CT1/11529-CH-8_sample_filtered_feature_bc_matrix"))
obj3 <- CreateSeuratObject(counts = Read10X(data.dir = "R:/hemonc-burness/Burness_Lab/Individual Data Folders/Deeksha Sharma/7.R DATA and Scripts/Burness Data/D2A1_CT1/11529-CH-9_sample_filtered_feature_bc_matrix"))
obj4 <- CreateSeuratObject(counts = Read10X(data.dir = "R:/hemonc-burness/Burness_Lab/Individual Data Folders/Deeksha Sharma/7.R DATA and Scripts/Burness Data/D2A1_CT1/11529-CH-10_sample_filtered_feature_bc_matrix"))
obj5 <- CreateSeuratObject(counts = Read10X(data.dir = "R:/hemonc-burness/Burness_Lab/Individual Data Folders/Deeksha Sharma/7.R DATA and Scripts/Burness Data/D2A1_CT1/11529-CH-11_sample_filtered_feature_bc_matrix"))
obj6 <- CreateSeuratObject(counts = Read10X(data.dir = "R:/hemonc-burness/Burness_Lab/Individual Data Folders/Deeksha Sharma/7.R DATA and Scripts/Burness Data/D2A1_CT1/11529-CH-12_sample_filtered_feature_bc_matrix"))

# 3.1 Define sample names
sample_names <- c("obj1", "obj2", "obj3", "obj4", "obj5", "obj6")

# 3.2 Merge Seurat objects
merged_samples <- merge(x = obj1, y = list(obj2, obj3, obj4, obj5, obj6), add.cell.ids = sample_names)

# 3.3 Check the merged object
print(merged_samples)

# 3.4 Save the merged data
saveRDS(object = merged_samples, file = paste0(getwd(), "/processed-data/", "all_1-merged-object.RDS"))

# 3.5 Join layers after merging
obj <- JoinLayers(merged_samples)

# 3.6 Define new metadata for control and treatment condition

# 3.7 Remove individual objects from the environment
rm(list = c("obj1", "obj2", "obj3", "obj4", "obj5", "obj6"))

# 3.8 Save the merged data
saveRDS(object = obj, file = paste0(getwd(), "/processed-data/", "all_1-joint-object.RDS"))

# 4. Add metadata (sample, cell line, and treatment)
library(stringr)
sample.by.barcode <- word(rownames(obj[[]]), start = 1, end = 1, sep = fixed("_"))
obj$sample <- sample.by.barcode
cell.line <- rep(NA, length(sample.by.barcode))
cell.line[sample.by.barcode %in% c("11529-CH-7", "11529-CH-8", "11529-CH-9", "11529-CH-10", "11529-CH-11", "11529-CH-12")] <- "D2A1"
obj$cell_line <- cell.line
treatment <- rep(NA, length(sample.by.barcode))
treatment[sample.by.barcode %in% c("11529-CH-7", "11529-CH-8", "11529-CH-9")] <- "Vehicle"
treatment[sample.by.barcode %in% c("11529-CH-10", "11529-CH-11", "11529-CH-12")] <- "ZBC260"
obj$treatment <- treatment

# 4.9 Save the merged object with added metadata
saveRDS(object = obj, file = paste0(getwd(), "/processed-data/", "D2A1_1-merged-object.RDS"))

# 5. QC filtering
obj <- readRDS(file = paste0(getwd(), "/processed-data/", "D2A1_1-merged-object.RDS"))
ncells.before.qc.filter <- dim(obj)[2]
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# 5.5 Visualize QC metrics
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# 5.9 Save QC results
write.csv(qc_results, file = paste0(tables_dir, "1-qc_filtering_results.csv"), row.names = FALSE)
ggsave(filename = paste0(figures_dir, "1-QCmetrics-FeatureScatter.png"))

# 6. Perform filtering
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 100000 & percent.mt < 10)
ncells.after.qc.filter <- dim(obj)[2]
ncells.filtered <- ncells.before.qc.filter - ncells.after.qc.filter
percent.filtered <- ncells.filtered / ncells.before.qc.filter * 100

# 6.3 Create a data frame for QC results
qc_results <- data.frame(Metric = c("Number of cells before QC", "Number of cells after QC", "Number of cells filtered", "Percent filtered"),
                         Value = c(ncells.before.qc.filter, ncells.after.qc.filter, ncells.filtered, round(percent.filtered, 2)))
write.csv(qc_results, file = paste0(tables_dir, "1-qc_filtering_results.csv"), row.names = FALSE)

# 6.5 Save RDS file after filtration
saveRDS(object = obj, file = paste0(getwd(), "/processed-data/", "d2a1_2-QCfiltered-object.RDS"))

# 7. SCTransform
obj <- readRDS(file = paste0(getwd(), "/processed-data/", "d2a1_2-QCfiltered-object.RDS"))
obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = TRUE)

# Save after SCTransform
saveRDS(object = obj, file = paste0(getwd(), "/processed-data/", "d2a1_3-SCtransform-object.RDS"))

# 8. Clustering and UMAP
obj <- readRDS(file = paste0(getwd(), "/processed-data/", "d2a1_3-SCtransform-object.RDS"))
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, k.param = 20)
obj <- FindClusters(obj, resolution = 1)

# Visualize clusters
DimPlot(obj, label = TRUE)

# 9. Find markers for all clusters
all_markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save markers to Excel
write.xlsx(all_markers, "all_markers_Myeloid.xlsx")

# 9.2 Create a DotPlot for the top markers
top3_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DotPlot(obj, features = top3_markers$gene) + ggtitle("Top 3 Markers per Cluster") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 9.3 Biased Feature Plot (e.g., specific gene expression)
FeaturePlot(object = obj, features = c("Krt18", "Csfr1r", "Eng", "Col3a1", "Cd3e", "Igkc"))

# 9.4 Rename clusters based on biological interpretation
new.cluster.ids <- c("Tumor cells", "Tumor cells", "Tumor cells", "Myeloid cells", "Tumor cells", "Tumor cells", 
                     "Tumor cells", "Fibroblast", "Tumor cells", "Myeloid cells", "Endothelial cells", "Fibroblast", 
                     "Lymphoid cells", "Myeloid cells", "Myeloid cells", "Fibroblast", "Tumor cells", "Lymphoid cells", "Myeloid cells")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

# Plot UMAP with renamed clusters
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.5)

# 10. Split clusters by treatment
obj@meta.data$orig.ident <- sub("_.*", "", rownames(obj@meta.data))
treatment <- rep(NA, nrow(obj@meta.data))
treatment[obj@meta.data$orig.ident %in% c("obj1", "obj2", "obj3")] <- "Vehicle"
treatment[obj@meta.data$orig.ident %in% c("obj4", "obj5", "obj6")] <- "ZBC260"
obj@meta.data$treatment <- treatment

# Save Seurat object with treatment information
saveRDS(obj, file = "D2A1_obj_with_treatment.rds")

# Split UMAP visualization by treatment
DimPlot(obj, reduction = "umap", split.by = "treatment") + ggtitle("UMAP: Control vs BETd Conditions") + theme_minimal()
