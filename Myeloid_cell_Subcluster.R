## ================================================
## LIBRARIES
## ================================================
# Core
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(tibble)
library(tidyverse)
library(stringr)
library(parallel)

# Visualization
library(pheatmap)
library(ggridges)

# Single-cell & Bioinformatics
library(scCATCH)
library(Nebulosa)
library(future)
library(cluster)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(glmGamPoi)

# Enrichment
library(enrichR)

# I/O
library(openxlsx)
library(writexl)
library(readr)
library(readxl)

# Dev
library(devtools)
library(BiocManager)

## ================================================
## LOAD AND SUBSET SEURAT OBJECT
## ================================================
# Subset Seurat object to Myeloid cells
myeloid_cells <- subset(obj, idents = "Myeloid cells")
dim(myeloid_cells)
DimPlot(myeloid_cells, reduction = "umap") + ggtitle("UMAP of Myeloid Cells")

## Check available treatment conditions
unique(myeloid_cells@meta.data$treatment)

## ================================================
## NORMALIZATION & SAVING
## ================================================
# SCTransform normalization
myeloid_cells <- SCTransform(myeloid_cells, vars.to.regress = "percent.mt", verbose = TRUE)

# Save RDS after SCTransform
saveRDS(myeloid_cells, file = paste0(getwd(), "/processed-data/d2a1_myeloid_cells_SCtransform-myeloid_myeloid_object.RDS"))

# Load the RDS object (for future sessions)
myeloid_cells <- readRDS(file = paste0(getwd(), "/processed-data/d2a1_myeloid_cells_SCtransform-myeloid_myeloid_object.RDS"))

## ================================================
## DIMENSION REDUCTION & CLUSTERING
## ================================================
myeloid_cells <- RunPCA(myeloid_cells, verbose = FALSE)
myeloid_cells <- RunUMAP(myeloid_cells, dims = 1:30, verbose = FALSE)
myeloid_cells <- FindNeighbors(myeloid_cells, dims = 1:30, verbose = FALSE)
myeloid_cells <- FindClusters(myeloid_cells, resolution = 0.3, verbose = FALSE)

# Plot clusters
DimPlot(myeloid_cells, pt.size = 1)
DimPlot(myeloid_cells, split.by = "treatment", pt.size = 1)

## ================================================
## FIND CLUSTER MARKERS
## ================================================
all_markers <- FindAllMarkers(
  object = myeloid_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Get top 10 markers per cluster
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

head(top10_markers, 20)

# Heatmap of top 10 markers
DoHeatmap(myeloid_cells, features = top10_markers$gene) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 10)
  )

## ================================================
## AVERAGE EXPRESSION HEATMAP (SCALED)
## ================================================
avg_expr <- AverageExpression(myeloid_cells, features = top10_markers$gene, assays = "SCT")$SCT
scaled_expr <- t(scale(t(avg_expr)))

pheatmap(
  scaled_expr,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Top 10 Markers Per Cluster (SCT Data)"
)

## ================================================
## CELL COUNT & PROPORTION ANALYSIS
## ================================================
# Add cluster IDs to metadata
myeloid_cells@meta.data$new.cluster.ids <- Idents(myeloid_cells)
metadata <- myeloid_cells@meta.data

# Count cells per treatment and cluster
cell_counts <- metadata %>%
  group_by(new.cluster.ids, treatment) %>%
  summarise(cell_count = n()) %>%
  ungroup()

print(cell_counts)

# Save counts
write.xlsx(cell_counts, "cell_counts.xlsx", rowNames = FALSE)
saveRDS(cell_counts, "cell_counts.rds")

# Load back if needed
cell_counts <- read.xlsx("cell_counts.xlsx")
cell_counts <- readRDS("cell_counts.rds")
print(cell_counts)

# Calculate proportions
cell_proportions <- cell_counts %>%
  group_by(treatment) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup()

## ================================================
## VISUALIZATION: COUNTS & PROPORTIONS
## ================================================
# Stacked barplot - cell counts
ggplot(cell_counts, aes(x = treatment, y = cell_count, fill = new.cluster.ids)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Treatment",
    y = "Count of Cells",
    fill = "Cluster",
    title = "Count of Cells in Each Cluster by Treatment"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14)
  )

# Stacked barplot - proportions
ggplot(cell_proportions, aes(x = treatment, y = proportion, fill = new.cluster.ids)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Treatment",
    y = "Proportion of Cells",
    fill = "Cluster",
    title = "Proportion of Cells in Each Cluster by Treatment"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14)
  )

# Dodged barplot - per cluster
ggplot(cell_counts, aes(x = new.cluster.ids, y = cell_count, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    x = "Cluster",
    y = "Cell Count",
    fill = "Treatment",
    title = "Cell Count Distribution by Cluster and Treatment"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
