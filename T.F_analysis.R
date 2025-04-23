# 1. Load Required Libraries

# Install and load necessary packages only once
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AUCell", "RcisTarget", "GENIE3", "zoo", "mixtools", "rbokeh", 
                       "DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne", 
                       "doMC", "doRNG", "devtools"))
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SingleCellExperiment)
library(doParallel)
library(foreach)
library(Matrix)
library(writexl)

# Optionally install from GitHub (if necessary)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC") 

# Check installed versions
packageVersion("SCENIC")
# 2. Subset Myeloid Cells
# Check available identities (clusters)
unique(Idents(obj))

# Get Myeloid cell names
myeloid_cells <- WhichCells(obj, ident = "Myeloid cells")

# Remove prefixes to match exprMatrix column names
prefix_map <- c(
  "obj1_" = "D2A1_S7_", "obj2_" = "D2A1_S8_", "obj3_" = "D2A1_S9_", 
  "obj4_" = "D2A1_S10_", "obj5_" = "D2A1_S11_", "obj6_" = "D2A1_S12_"
)

# Apply the prefix map to Myeloid cells
myeloid_cells_fixed <- myeloid_cells
for (prefix in names(prefix_map)) {
  myeloid_cells_fixed <- gsub(prefix, prefix_map[prefix], myeloid_cells_fixed)
}

# Check intersection between fixed cell names and exprMatrix columns
exprMatrix_barcodes <- sub("^D2A1_S[0-9]+_", "", colnames(exprMatrix))
matched_cells <- colnames(exprMatrix)[exprMatrix_barcodes %in% sub("^D2A1_S[0-9]+_", "", myeloid_cells_fixed)]

# Subset the expression matrix for myeloid cells
exprMatrix_myeloid <- exprMatrix[, matched_cells]

# Save the subset matrix
saveRDS(exprMatrix_myeloid, file = "exprMatrix_myeloid.rds")
write.csv(as.matrix(exprMatrix_myeloid), "exprMatrix_myeloid.csv")

# 3. Filter Top 2000 Variable Genes
exprMatrix_myeloid <- readRDS("exprMatrix_myeloid.rds")
gene_variances <- apply(exprMatrix_myeloid, 1, var)
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:2000]
exprMatrix_myeloid_top2000 <- exprMatrix_myeloid[top_genes, ]

# Save filtered matrix
saveRDS(exprMatrix_myeloid_top2000, file = "exprMatrix_myeloid_top2000.rds")
write.csv(as.matrix(exprMatrix_myeloid_top2000), "exprMatrix_myeloid_top2000.csv")

# 4. Filter for Genes Kept after SCENIC Preprocessing
genesKept <- readRDS("genesKept_myeloid.rds")
common_genes <- intersect(genesKept, rownames(exprMatrix_myeloid))
exprMat_filtered <- exprMatrix_myeloid[common_genes, ]
saveRDS(exprMat_filtered, file = "exprMat_filtered.rds")
write.csv(as.matrix(exprMat_filtered), "exprMat_filtered.csv")

# 5. Correlation Analysis
exprMat_filtered_log <- log2(exprMat_filtered + 1)

# Save log-transformed data
saveRDS(exprMat_filtered_log, file = "exprMat_filtered_log.rds")
write.csv(as.matrix(exprMat_filtered_log), "exprMat_filtered_log.csv")

# 6. GENIE3 Analysis with Parallelization
exprMat_filtered_log <- readRDS("exprMat_filtered_log.rds")
numCores <- detectCores() - 2  # Use available cores efficiently
cl <- makeCluster(numCores)
registerDoParallel(cl)

scenicOptions@settings$nCores <- numCores
runGenie3(exprMat_filtered_log, scenicOptions)

stopCluster(cl)
saveRDS(scenicOptions, file = "scenicOptions_Genie3.rds")

# 7. Build and Score the GRN
if (file.exists("scenicOptions_Genie3.rds")) {
  scenicOptions <- readRDS("scenicOptions_Genie3.rds")
}

scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123  # Set seed for reproducibility
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$splitBy <- "chunk"

saveRDS(scenicOptions, file = "scenicOptions_optimized.rds")

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file = "scenicOptions_step1.rds")

scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod = "top5perTarget")
saveRDS(scenicOptions, file = "scenicOptions_step2.rds")

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
saveRDS(scenicOptions, file = "scenicOptions_step3.rds")

# 8. Visualization and Results Check
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
auc_matrix <- as.matrix(getAUC(regulonAUC))

# Plot the top 20 most variable regulons
top_regulons <- names(tail(sort(apply(auc_matrix, 1, var)), 20))
pheatmap(auc_matrix[top_regulons, ], main = "Top 20 Regulon Activities")

# 9. Subset AUC Matrix by Conditions (Vehicle vs ZBC260)
dataset_ids <- sub("^(D2A1_S[0-9]+)_.*", "\\1", colnames(auc_matrix))
vehicle_cells <- colnames(auc_matrix)[dataset_ids %in% c("D2A1_S7", "D2A1_S8", "D2A1_S9")]
zbc260_cells <- colnames(auc_matrix)[dataset_ids %in% c("D2A1_S10", "D2A1_S11", "D2A1_S12")]

# Save condition-specific AUC matrices
auc_vehicle <- auc_matrix[, vehicle_cells, drop=FALSE]
auc_zbc260 <- auc_matrix[, zbc260_cells, drop=FALSE]
saveRDS(auc_vehicle, file = "auc_vehicle.rds")
saveRDS(auc_zbc260, file = "auc_zbc260.rds")

# Combine AUC data and plot
auc_combined <- data.frame(
  Regulon = rep(rownames(auc_vehicle), 2),
  AUC = c(rowMeans(auc_vehicle), rowMeans(auc_zbc260)),
  Condition = rep(c("Vehicle", "ZBC260"), each = nrow(auc_vehicle))
)

top_regulons <- unique(c(top10_vehicle, top10_zbc260))
auc_top <- auc_combined[auc_combined$Regulon %in% top_regulons, ]

# Violin plot
ggplot(auc_top, aes(x = Regulon, y = AUC, fill = Condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top") +
  labs(title = "Top 10 Regulon Activity (Vehicle vs ZBC260)",
       x = "Regulon", y = "AUC Score") +
  scale_fill_manual(values = c("Vehicle" = "#E69F00", "ZBC260" = "#0072B2")) 