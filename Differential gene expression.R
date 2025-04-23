
# Assuming you ran differential expression analysis with Seurat or a similar package
cluster0.markers <- FindMarkers(myeloid_cells, ident.1 = "0")

# Check the structure of cluster0.markers again
str(cluster0.markers)

# Filter for significant upregulated genes (adjusted p-value < 0.05 and avg_log2FC > 1)
significant_upregulated_genes <- cluster0.markers[
  cluster0.markers$p_val_adj < 0.05 & cluster0.markers$avg_log2FC > 0.1, 
]
# Save the significant upregulated genes as RDS file
saveRDS(significant_upregulated_genes, "significant_upregulated_genes_c0.rds")

# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_upregulated_genes)

# Filter for significant upregulated genes (adjusted p-value < 0.05 and avg_log2FC > 1)
significant_downregulated_genes <- cluster0.markers[
  cluster0.markers$p_val_adj < 0.05 & cluster0.markers$avg_log2FC <0.1, 
]

# Check how many genes are selected
nrow(significant_downregulated_genes)

# Preview the first few genes
head(significant_downregulated_genes)


# Extract the gene symbols (row names) for the significant upregulated genes
upregulated_genes <- rownames(significant_upregulated_genes)

# Define the Enrichr databases to query
dbs <- c("GO_Biological_Process_2021", 
         "KEGG_2021_Mouse", 
         "WikiPathways_2021_Mouse")

# Perform enrichment analysis using Enrichr for the upregulated genes
enrich_results <- enrichr(upregulated_genes, dbs)

# Check column names for KEGG, GO, and WikiPathways results
print(colnames(enrich_results[["KEGG_2021_Mouse"]]))
print(colnames(enrich_results[["GO_Biological_Process_2021"]]))
print(colnames(enrich_results[["WikiPathways_2021_Mouse"]]))

# Save all results as RDS files
saveRDS(enrich_results[["KEGG_2021_Mouse"]], "KEGG_2021_Mouse_c0_upregulated.rds")
saveRDS(enrich_results[["GO_Biological_Process_2021"]], "GO_Biological_Process_2021_c0_upregulated.rds")
saveRDS(enrich_results[["WikiPathways_2021_Mouse"]], "WikiPathways_2021_Mouse_c0_upregulated.rds")

# Save all results as Excel files using writexl package
library(writexl)

# Combine the data into a list of data frames for easy saving
results_list <- list(
  KEGG_2021_Mouse = enrich_results[["KEGG_2021_Mouse"]],
  GO_Biological_Process_2021 = enrich_results[["GO_Biological_Process_2021"]],
  WikiPathways_2021_Mouse = enrich_results[["WikiPathways_2021_Mouse"]]
)

# Write to Excel file
write_xlsx(results_list, "Enrichment_Results_c0_upregulated.xlsx")

# Visualize top enriched pathways using barplot (based on p-value)
barplot(-log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10]), 
        names.arg = enrich_results[["KEGG_2021_Mouse"]]$Term[1:10], 
        las = 2, col = "blue", main = "Top 10 KEGG Pathways (Upregulated)", ylab = "-log10(p-value)")

# Adjust margin size to prevent text cutoff
par(mar = c(5, 12, 3, 2))  # Increase left margin for pathway names

# Extract top 10 pathways and their adjusted p-values
top_terms <- enrich_results[["KEGG_2021_Mouse"]]$Term[1:10]
top_pvals <- -log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10])

# Define the Cluster 0 color (coral color)
cluster0_color <- "#F8766D"  # Soft coral color

# Adjust margins for better visibility
par(mar = c(5, 15, 3, 2))  # Increase left margin for pathway names

# KEGG Pathways Bar Plot
barplot(top_pvals, 
        names.arg = top_terms, 
        horiz = TRUE, 
        las = 1,  # Keep y-axis labels readable
        col = cluster0_color,  # Use Cluster 0 color (coral)
        main = "Top 10 KEGG Pathways (Upregulated)", 
        xlab = "-log10(p-value)", 
        cex.lab = 1.2,  # Increase axis label size
        cex.names = 0.8)  # Adjust text size if necessary

# Add an x-axis
axis(1)  

# Add a vertical axis line at x = 0
abline(v = 0, col = "black", lwd = 2)



#########-----------CLUSTER1##############################################################################################


# Assuming you ran differential expression analysis with Seurat or a similar package
cluster1.markers <- FindMarkers(myeloid_cells, ident.1 = "1")  # Change to Cluster 1

# Check the structure of cluster1.markers
str(cluster1.markers)

# Filter for significant upregulated genes (adjusted p-value < 0.05 and avg_log2FC > 0.1)
significant_upregulated_genes <- cluster1.markers[
  cluster1.markers$p_val_adj < 0.05 & cluster1.markers$avg_log2FC > 0.1, 
]
# Save the significant upregulated genes as RDS file
saveRDS(significant_upregulated_genes, "significant_upregulated_genes_c1.rds")


# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_upregulated_genes)

# Filter for significant downregulated genes (adjusted p-value < 0.05 and avg_log2FC < 0.1)
significant_downregulated_genes <- cluster1.markers[
  cluster1.markers$p_val_adj < 0.05 & cluster1.markers$avg_log2FC < 0.1, 
]

# Check how many genes are selected
nrow(significant_downregulated_genes)

# Preview the first few genes
head(significant_downregulated_genes)

# Extract the gene symbols (row names) for the significant upregulated genes
upregulated_genes <- rownames(significant_upregulated_genes)

# Define the Enrichr databases to query
dbs <- c("GO_Biological_Process_2021", 
         "KEGG_2021_Mouse", 
         "WikiPathways_2021_Mouse")

# Perform enrichment analysis using Enrichr for the upregulated genes
enrich_results <- enrichr(upregulated_genes, dbs)

# Check column names for KEGG, GO, and WikiPathways results
print(colnames(enrich_results[["KEGG_2021_Mouse"]]))
print(colnames(enrich_results[["GO_Biological_Process_2021"]]))
print(colnames(enrich_results[["WikiPathways_2021_Mouse"]]))

# Save all results as RDS files
saveRDS(enrich_results[["KEGG_2021_Mouse"]], "KEGG_2021_Mouse_c1_upregulated.rds")
saveRDS(enrich_results[["GO_Biological_Process_2021"]], "GO_Biological_Process_2021_c1_upregulated.rds")
saveRDS(enrich_results[["WikiPathways_2021_Mouse"]], "WikiPathways_2021_Mouse_c1_upregulated.rds")

# Save all results as Excel files using writexl package
library(writexl)

# Combine the data into a list of data frames for easy saving
results_list <- list(
  KEGG_2021_Mouse = enrich_results[["KEGG_2021_Mouse"]],
  GO_Biological_Process_2021 = enrich_results[["GO_Biological_Process_2021"]],
  WikiPathways_2021_Mouse = enrich_results[["WikiPathways_2021_Mouse"]]
)

# Write to Excel file
write_xlsx(results_list, "Enrichment_Results_c1_upregulated.xlsx")

# Visualize top enriched pathways using barplot (based on p-value)
barplot(-log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10]), 
        names.arg = enrich_results[["KEGG_2021_Mouse"]]$Term[1:10], 
        las = 2, col = "blue", main = "Top 10 KEGG Pathways (Cluster 1 - Upregulated)", ylab = "-log10(p-value)")

# Adjust margin size to prevent text cutoff
par(mar = c(5, 12, 3, 2))  # Increase left margin for pathway names

# Extract top 10 pathways and their adjusted p-values
top_terms <- enrich_results[["KEGG_2021_Mouse"]]$Term[1:10]
top_pvals <- -log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10])

# Define the Cluster 1 color (Olive Green color)
cluster1_color <- "#aba300"  # Olive Green

# Adjust margins for better visibility
par(mar = c(5, 15, 3, 2))  # Increase left margin for pathway names

# KEGG Pathways Bar Plot
barplot(top_pvals, 
        names.arg = top_terms, 
        horiz = TRUE, 
        las = 1,  # Keep y-axis labels readable
        col = cluster1_color,  # Use Cluster 1 color (Olive Green)
        main = "Top 10 KEGG Pathways (Upregulated) - Cluster 1", 
        xlab = "-log10(p-value)", 
        cex.lab = 1.2,  # Increase axis label size
        cex.names = 0.8)  # Adjust text size if necessary

# Add an x-axis
axis(1)  

# Add a vertical axis line at x = 0
abline(v = 0, col = "black", lwd = 2)






#########---CLUSTER2-----

# Assuming you ran differential expression analysis with Seurat or a similar package
cluster2.markers <- FindMarkers(myeloid_cells, ident.1 = "2")  # Change to Cluster 2

# Check the structure of cluster2.markers
str(cluster2.markers)

# Filter for significant upregulated genes (adjusted p-value < 0.05 and avg_log2FC > 0.1)
significant_upregulated_genes <- cluster2.markers[
  cluster2.markers$p_val_adj < 0.05 & cluster2.markers$avg_log2FC > 0.1, 
]

# Save the significant upregulated genes as RDS file
saveRDS(significant_upregulated_genes, "significant_upregulated_genes_c2.rds")

# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_upregulated_genes)

# Filter for significant downregulated genes (adjusted p-value < 0.05 and avg_log2FC < 0.1)
significant_downregulated_genes <- cluster2.markers[
  cluster2.markers$p_val_adj < 0.05 & cluster2.markers$avg_log2FC < 0.1, 
]

# Save the significant downregulated genes as RDS file
saveRDS(significant_downregulated_genes, "significant_downregulated_genes_c2.rds")

# Check how many genes are selected
nrow(significant_downregulated_genes)
# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_downregulated_genes)

# Extract the gene symbols (row names) for the significant upregulated genes
upregulated_genes <- rownames(significant_upregulated_genes)

# Define the Enrichr databases to query
dbs <- c("GO_Biological_Process_2021", 
         "KEGG_2021_Mouse", 
         "WikiPathways_2021_Mouse")

# Perform enrichment analysis using Enrichr for the upregulated genes
enrich_results <- enrichr(upregulated_genes, dbs)

# Check column names for KEGG, GO, and WikiPathways results
print(colnames(enrich_results[["KEGG_2021_Mouse"]]))
print(colnames(enrich_results[["GO_Biological_Process_2021"]]))
print(colnames(enrich_results[["WikiPathways_2021_Mouse"]]))

# Save all results as RDS files
saveRDS(enrich_results[["KEGG_2021_Mouse"]], "KEGG_2021_Mouse_c2_upregulated.rds")
saveRDS(enrich_results[["GO_Biological_Process_2021"]], "GO_Biological_Process_2021_c2_upregulated.rds")
saveRDS(enrich_results[["WikiPathways_2021_Mouse"]], "WikiPathways_2021_Mouse_c2_upregulated.rds")

# Save all results as Excel files using writexl package
library(writexl)

# Combine the data into a list of data frames for easy saving
results_list <- list(
  KEGG_2021_Mouse = enrich_results[["KEGG_2021_Mouse"]],
  GO_Biological_Process_2021 = enrich_results[["GO_Biological_Process_2021"]],
  WikiPathways_2021_Mouse = enrich_results[["WikiPathways_2021_Mouse"]]
)

# Write to Excel file
write_xlsx(results_list, "Enrichment_Results_c2_upregulated.xlsx")

# Visualize top enriched pathways using barplot (based on p-value)
barplot(-log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10]), 
        names.arg = enrich_results[["KEGG_2021_Mouse"]]$Term[1:10], 
        las = 2, col = "blue", main = "Top 10 KEGG Pathways (Cluster 2 - Upregulated)", ylab = "-log10(p-value)")

# Adjust margin size to prevent text cutoff
par(mar = c(5, 12, 3, 2))  # Increase left margin for pathway names

# Extract top 10 pathways and their adjusted p-values
top_terms <- enrich_results[["KEGG_2021_Mouse"]]$Term[1:10]
top_pvals <- -log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10])

# Define the Cluster 2 color (Teal color)
cluster2_color <- "#7CAE00"  # Teal color for Cluster 2

# Adjust margins for better visibility
par(mar = c(5, 15, 3, 2))  # Increase left margin for pathway names

# KEGG Pathways Bar Plot
barplot(top_pvals, 
        names.arg = top_terms, 
        horiz = TRUE, 
        las = 1,  # Keep y-axis labels readable
        col = cluster2_color,  # Use Cluster 2 color (Teal)
        main = "Top 10 KEGG Pathways (Upregulated) - Cluster 2", 
        xlab = "-log10(p-value)", 
        cex.lab = 1.2,  # Increase axis label size
        cex.names = 0.8)  # Adjust text size if necessary

# Add an x-axis
axis(1)  

# Add a vertical axis line at x = 0
abline(v = 0, col = "black", lwd = 2)

# Save the raw markers as RDS file for future use
saveRDS(cluster2.markers, "cluster2_markers_raw.rds")












################Cluster 3---------------------



# Assuming you ran differential expression analysis with Seurat or a similar package
cluster3.markers <- FindMarkers(myeloid_cells, ident.1 = "3")  # Change to Cluster 3

# Check the structure of cluster3.markers
str(cluster3.markers)

# Filter for significant upregulated genes (adjusted p-value < 0.05 and avg_log2FC > 0.1)
significant_upregulated_genes <- cluster3.markers[
  cluster3.markers$p_val_adj < 0.05 & cluster3.markers$avg_log2FC > 0.1, 
]

# Save the significant upregulated genes as RDS file
saveRDS(significant_upregulated_genes, "significant_upregulated_genes_c3.rds")

# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_upregulated_genes)

# Filter for significant downregulated genes (adjusted p-value < 0.05 and avg_log2FC < 0.1)
significant_downregulated_genes <- cluster3.markers[
  cluster3.markers$p_val_adj < 0.05 & cluster3.markers$avg_log2FC < 0.1, 
]

# Save the significant downregulated genes as RDS file
saveRDS(significant_downregulated_genes, "significant_downregulated_genes_c3.rds")

# Check how many genes are selected
nrow(significant_downregulated_genes)
# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_downregulated_genes)

# Extract the gene symbols (row names) for the significant upregulated genes
upregulated_genes <- rownames(significant_upregulated_genes)

# Define the Enrichr databases to query
dbs <- c("GO_Biological_Process_2021", 
         "KEGG_2021_Mouse", 
         "WikiPathways_2021_Mouse")

# Perform enrichment analysis using Enrichr for the upregulated genes
enrich_results <- enrichr(upregulated_genes, dbs)

# Check column names for KEGG, GO, and WikiPathways results
print(colnames(enrich_results[["KEGG_2021_Mouse"]]))
print(colnames(enrich_results[["GO_Biological_Process_2021"]]))
print(colnames(enrich_results[["WikiPathways_2021_Mouse"]]))

# Save all results as RDS files
saveRDS(enrich_results[["KEGG_2021_Mouse"]], "KEGG_2021_Mouse_c3_upregulated.rds")
saveRDS(enrich_results[["GO_Biological_Process_2021"]], "GO_Biological_Process_2021_c3_upregulated.rds")
saveRDS(enrich_results[["WikiPathways_2021_Mouse"]], "WikiPathways_2021_Mouse_c3_upregulated.rds")

# Save all results as Excel files using writexl package
library(writexl)

# Combine the data into a list of data frames for easy saving
results_list <- list(
  KEGG_2021_Mouse = enrich_results[["KEGG_2021_Mouse"]],
  GO_Biological_Process_2021 = enrich_results[["GO_Biological_Process_2021"]],
  WikiPathways_2021_Mouse = enrich_results[["WikiPathways_2021_Mouse"]]
)

# Write to Excel file
write_xlsx(results_list, "Enrichment_Results_c3_upregulated.xlsx")

# Visualize top enriched pathways using barplot (based on p-value)
barplot(-log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10]), 
        names.arg = enrich_results[["KEGG_2021_Mouse"]]$Term[1:10], 
        las = 2, col = "blue", main = "Top 10 KEGG Pathways (Cluster 3 - Upregulated)", ylab = "-log10(p-value)")

# Adjust margin size to prevent text cutoff
par(mar = c(5, 12, 3, 2))  # Increase left margin for pathway names

# Extract top 10 pathways and their adjusted p-values
top_terms <- enrich_results[["KEGG_2021_Mouse"]]$Term[1:10]
top_pvals <- -log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10])

# Define the Cluster 3 color (Light Blue color)
cluster3_color <- "#00C19A"  # Light Blue color for Cluster 3

# Adjust margins for better visibility
par(mar = c(5, 15, 3, 2))  # Increase left margin for pathway names

# KEGG Pathways Bar Plot
barplot(top_pvals, 
        names.arg = top_terms, 
        horiz = TRUE, 
        las = 1,  # Keep y-axis labels readable
        col = cluster3_color,  # Use Cluster 3 color (Light Blue)
        main = "Top 10 KEGG Pathways (Upregulated) - Cluster 3", 
        xlab = "-log10(p-value)", 
        cex.lab = 1.2,  # Increase axis label size
        cex.names = 0.8)  # Adjust text size if necessary

# Add an x-axis
axis(1)  

# Add a vertical axis line at x = 0
abline(v = 0, col = "black", lwd = 2)

# Save the raw markers as RDS file for future use
saveRDS(cluster3.markers, "cluster3_markers_raw.rds")







#############------------------Cluster 4--------------##############################################



# Assuming you ran differential expression analysis with Seurat or a similar package
cluster4.markers <- FindMarkers(myeloid_cells, ident.1 = "4")  # Change to Cluster 4

# Check the structure of cluster4.markers
str(cluster4.markers)

# Filter for significant upregulated genes (adjusted p-value < 0.05 and avg_log2FC > 0.1)
significant_upregulated_genes <- cluster4.markers[
  cluster4.markers$p_val_adj < 0.05 & cluster4.markers$avg_log2FC > 0.1, 
]

# Save the significant upregulated genes as RDS file
saveRDS(significant_upregulated_genes, "significant_upregulated_genes_c4.rds")

# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_upregulated_genes)

# Filter for significant downregulated genes (adjusted p-value < 0.05 and avg_log2FC < 0.1)
significant_downregulated_genes <- cluster4.markers[
  cluster4.markers$p_val_adj < 0.05 & cluster4.markers$avg_log2FC < 0.1, 
]

# Save the significant downregulated genes as RDS file
saveRDS(significant_downregulated_genes, "significant_downregulated_genes_c4.rds")

# Check how many genes are selected
nrow(significant_downregulated_genes)
# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_downregulated_genes)

# Extract the gene symbols (row names) for the significant upregulated genes
upregulated_genes <- rownames(significant_upregulated_genes)

# Define the Enrichr databases to query
dbs <- c("GO_Biological_Process_2021", 
         "KEGG_2021_Mouse", 
         "WikiPathways_2021_Mouse")

# Perform enrichment analysis using Enrichr for the upregulated genes
enrich_results <- enrichr(upregulated_genes, dbs)

# Check column names for KEGG, GO, and WikiPathways results
print(colnames(enrich_results[["KEGG_2021_Mouse"]]))
print(colnames(enrich_results[["GO_Biological_Process_2021"]]))
print(colnames(enrich_results[["WikiPathways_2021_Mouse"]]))

# Save all results as RDS files
saveRDS(enrich_results[["KEGG_2021_Mouse"]], "KEGG_2021_Mouse_c4_upregulated.rds")
saveRDS(enrich_results[["GO_Biological_Process_2021"]], "GO_Biological_Process_2021_c4_upregulated.rds")
saveRDS(enrich_results[["WikiPathways_2021_Mouse"]], "WikiPathways_2021_Mouse_c4_upregulated.rds")

# Save all results as Excel files using writexl package
library(writexl)

# Combine the data into a list of data frames for easy saving
results_list <- list(
  KEGG_2021_Mouse = enrich_results[["KEGG_2021_Mouse"]],
  GO_Biological_Process_2021 = enrich_results[["GO_Biological_Process_2021"]],
  WikiPathways_2021_Mouse = enrich_results[["WikiPathways_2021_Mouse"]]
)

# Write to Excel file
write_xlsx(results_list, "Enrichment_Results_c4_upregulated.xlsx")

# Visualize top enriched pathways using barplot (based on p-value)
barplot(-log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10]), 
        names.arg = enrich_results[["KEGG_2021_Mouse"]]$Term[1:10], 
        las = 2, col = "blue", main = "Top 10 KEGG Pathways (Cluster 4 - Upregulated)", ylab = "-log10(p-value)")

# Adjust margin size to prevent text cutoff
par(mar = c(5, 12, 3, 2))  # Increase left margin for pathway names

# Extract top 10 pathways and their adjusted p-values
top_terms <- enrich_results[["KEGG_2021_Mouse"]]$Term[1:10]
top_pvals <- -log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10])

# Define the Cluster 4 color (Light Purple color)
cluster4_color <- "#00A9FF"  # Light Purple color for Cluster 4

# Adjust margins for better visibility
par(mar = c(5, 15, 3, 2))  # Increase left margin for pathway names

# KEGG Pathways Bar Plot
barplot(top_pvals, 
        names.arg = top_terms, 
        horiz = TRUE, 
        las = 1,  # Keep y-axis labels readable
        col = cluster4_color,  # Use Cluster 4 color (Light Purple)
        main = "Top 10 KEGG Pathways (Upregulated) - Cluster 4", 
        xlab = "-log10(p-value)", 
        cex.lab = 1.2,  # Increase axis label size
        cex.names = 0.8)  # Adjust text

# Add an x-axis
axis(1)  

# Add a vertical axis line at x = 0
abline(v = 0, col = "black", lwd = 2)

# Save the raw markers as RDS file for future use
saveRDS(cluster4.markers, "cluster4_markers_raw.rds")

################------------CLUSTER5---------------


# Assuming you ran differential expression analysis with Seurat or a similar package
cluster5.markers <- FindMarkers(myeloid_cells, ident.1 = "5")  # Change to Cluster 5

# Check the structure of cluster5.markers
str(cluster5.markers)

# Filter for significant upregulated genes (adjusted p-value < 0.05 and avg_log2FC > 0.1)
significant_upregulated_genes <- cluster5.markers[
  cluster5.markers$p_val_adj < 0.05 & cluster5.markers$avg_log2FC > 0.1, 
]

# Save the significant upregulated genes as RDS file
saveRDS(significant_upregulated_genes, "significant_upregulated_genes_c5.rds")

# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_upregulated_genes)

# Filter for significant downregulated genes (adjusted p-value < 0.05 and avg_log2FC < 0.1)
significant_downregulated_genes <- cluster5.markers[
  cluster5.markers$p_val_adj < 0.05 & cluster5.markers$avg_log2FC < 0.1, 
]

# Save the significant downregulated genes as RDS file
saveRDS(significant_downregulated_genes, "significant_downregulated_genes_c5.rds")

# Check how many genes are selected
nrow(significant_downregulated_genes)

# Preview the first few genes
head(significant_downregulated_genes)

# Extract the gene symbols (row names) for the significant upregulated genes
upregulated_genes <- rownames(significant_upregulated_genes)

# Define the Enrichr databases to query
dbs <- c("GO_Biological_Process_2021", 
         "KEGG_2021_Mouse", 
         "WikiPathways_2021_Mouse")

# Perform enrichment analysis using Enrichr for the upregulated genes
enrich_results <- enrichr(upregulated_genes, dbs)

# Check column names for KEGG, GO, and WikiPathways results
print(colnames(enrich_results[["KEGG_2021_Mouse"]]))
print(colnames(enrich_results[["GO_Biological_Process_2021"]]))
print(colnames(enrich_results[["WikiPathways_2021_Mouse"]]))

# Save all results as RDS files
saveRDS(enrich_results[["KEGG_2021_Mouse"]], "KEGG_2021_Mouse_c5_upregulated.rds")
saveRDS(enrich_results[["GO_Biological_Process_2021"]], "GO_Biological_Process_2021_c5_upregulated.rds")
saveRDS(enrich_results[["WikiPathways_2021_Mouse"]], "WikiPathways_2021_Mouse_c5_upregulated.rds")

# Save all results as Excel files using writexl package
library(writexl)

# Combine the data into a list of data frames for easy saving
results_list <- list(
  KEGG_2021_Mouse = enrich_results[["KEGG_2021_Mouse"]],
  GO_Biological_Process_2021 = enrich_results[["GO_Biological_Process_2021"]],
  WikiPathways_2021_Mouse = enrich_results[["WikiPathways_2021_Mouse"]]
)

# Write to Excel file
write_xlsx(results_list, "Enrichment_Results_c5_upregulated.xlsx")

# Visualize top enriched pathways using barplot (based on p-value)
barplot(-log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10]), 
        names.arg = enrich_results[["KEGG_2021_Mouse"]]$Term[1:10], 
        las = 2, col = "blue", main = "Top 10 KEGG Pathways (Cluster 5 - Upregulated)", ylab = "-log10(p-value)")

# Adjust margin size to prevent text cutoff
par(mar = c(5, 12, 3, 2))  # Increase left margin for pathway names

# Extract top 10 pathways and their adjusted p-values
top_terms <- enrich_results[["KEGG_2021_Mouse"]]$Term[1:10]
top_pvals <- -log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10])

# Define the Cluster 5 color (Light Purple color)
cluster5_color <- "#8494FF"  # 

# Adjust margins for better visibility
par(mar = c(5, 15, 3, 2))  # Increase left margin for pathway names

# KEGG Pathways Bar Plot
barplot(top_pvals, 
        names.arg = top_terms, 
        horiz = TRUE, 
        las = 1,  # Keep y-axis labels readable
        col = cluster5_color,  # Use Cluster 5 color (Light Purple)
        main = "Top 10 KEGG Pathways (Upregulated) - Cluster 5", 
        xlab = "-log10(p-value)", 
        cex.lab = 1.2,  # Increase axis label size
        cex.names = 0.8)  # Adjust text
# Add an x-axis
axis(1)  

# Add a vertical axis line at x = 0
abline(v = 0, col = "black", lwd = 2)

# Save the raw markers as RDS file for future use
saveRDS(cluster5.markers, "cluster5_markers_raw.rds")

##############-------------CLUSTER 6##############################################


# Assuming you ran differential expression analysis with Seurat or a similar package
cluster6.markers <- FindMarkers(myeloid_cells, ident.1 = "6")  # Change to Cluster 6

# Check the structure of cluster6.markers
str(cluster6.markers)

# Filter for significant upregulated genes (adjusted p-value < 0.05 and avg_log2FC > 0.1)
significant_upregulated_genes <- cluster6.markers[
  cluster6.markers$p_val_adj < 0.05 & cluster6.markers$avg_log2FC > 0.1, 
]

# Save the significant upregulated genes as RDS file
saveRDS(significant_upregulated_genes, "significant_upregulated_genes_c6.rds")

# Check how many genes are selected
nrow(significant_upregulated_genes)

# Preview the first few genes
head(significant_upregulated_genes)

# Filter for significant downregulated genes (adjusted p-value < 0.05 and avg_log2FC < 0.1)
significant_downregulated_genes <- cluster6.markers[
  cluster6.markers$p_val_adj < 0.05 & cluster6.markers$avg_log2FC < 0.1, 
]

# Save the significant downregulated genes as RDS file
saveRDS(significant_downregulated_genes, "significant_downregulated_genes_c6.rds")

# Check how many genes are selected
nrow(significant_downregulated_genes)

# Preview the first few genes
head(significant_downregulated_genes)

# Extract the gene symbols (row names) for the significant upregulated genes
upregulated_genes <- rownames(significant_upregulated_genes)

# Define the Enrichr databases to query
dbs <- c("GO_Biological_Process_2021", 
         "KEGG_2021_Mouse", 
         "WikiPathways_2021_Mouse")

# Perform enrichment analysis using Enrichr for the upregulated genes
enrich_results <- enrichr(upregulated_genes, dbs)

# Check column names for KEGG, GO, and WikiPathways results
print(colnames(enrich_results[["KEGG_2021_Mouse"]]))
print(colnames(enrich_results[["GO_Biological_Process_2021"]]))
print(colnames(enrich_results[["WikiPathways_2021_Mouse"]]))

# Save all results as RDS files
saveRDS(enrich_results[["KEGG_2021_Mouse"]], "KEGG_2021_Mouse_c6_upregulated.rds")
saveRDS(enrich_results[["GO_Biological_Process_2021"]], "GO_Biological_Process_2021_c6_upregulated.rds")
saveRDS(enrich_results[["WikiPathways_2021_Mouse"]], "WikiPathways_2021_Mouse_c6_upregulated.rds")

# Save all results as Excel files using writexl package
library(writexl)

# Combine the data into a list of data frames for easy saving
results_list <- list(
  KEGG_2021_Mouse = enrich_results[["KEGG_2021_Mouse"]],
  GO_Biological_Process_2021 = enrich_results[["GO_Biological_Process_2021"]],
  WikiPathways_2021_Mouse = enrich_results[["WikiPathways_2021_Mouse"]]
)

# Write to Excel file
write_xlsx(results_list, "Enrichment_Results_c6_upregulated.xlsx")

# Visualize top enriched pathways using barplot (based on p-value)
barplot(-log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10]), 
        names.arg = enrich_results[["KEGG_2021_Mouse"]]$Term[1:10], 
        las = 2, col = "blue", main = "Top 10 KEGG Pathways (Cluster 6 - Upregulated)", ylab = "-log10(p-value)")

# Adjust margin size to prevent text cutoff
par(mar = c(5, 12, 3, 2))  # Increase left margin for pathway names

# Extract top 10 pathways and their adjusted p-values
top_terms <- enrich_results[["KEGG_2021_Mouse"]]$Term[1:10]
top_pvals <- -log10(enrich_results[["KEGG_2021_Mouse"]]$Adjusted.P.value[1:10])

# Define the Cluster 6 color (Light Purple color)
cluster6_color <- "#FF68A1"  # Light Purple color for Cluster 6

# Adjust margins for better visibility
par(mar = c(5, 15, 3, 2))  # Increase left margin for pathway names

# KEGG Pathways Bar Plot
barplot(top_pvals, 
        names.arg = top_terms, 
        horiz = TRUE, 
        las = 1,  # Keep y-axis labels readable
        col = cluster6_color,  # Use Cluster 6 color (Light Purple)
        main = "Top 10 KEGG Pathways (Upregulated) - Cluster 6", 
        xlab = "-log10(p-value)", 
        cex.lab = 1.2,  # Increase axis label size
        cex.names = 0.8)  # Adjust text



# Add an x-axis
axis(1)  

# Add a vertical axis line at x = 0
abline(v = 0, col = "black", lwd = 2)

# Save the raw markers as RDS file for future use
saveRDS(cluster6.markers, "cluster6_markers_raw.rds")

############Effect of treatment

# Load required libraries
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(openxlsx)

###-------------------#
### Differential Expression Analysis
###-------------------#

# Perform differential expression analysis between ZBC260 and Vehicle treatments
diff_expr_myeloid_cells <- FindMarkers(
  object = myeloid_cells, 
  ident.1 = "ZBC260", 
  ident.2 = "Vehicle", 
  group.by = "treatment"
)

# View the top results
head(diff_expr_myeloid_cells)

# Save results as RDS
saveRDS(diff_expr_myeloid_cells, file = "diff_expr_results_myeloid_cells1.rds")

# Load results if needed
diff_expr_myeloid_cells <- readRDS("diff_expr_results_myeloid_cells1.rds")

###-------------------#
### Filter Significant Genes
###-------------------#

# Filter for significant upregulated genes (adj p-value < 0.05 & log2FC > 0.1)
significant_upregulated_genes <- diff_expr_myeloid_cells[
  diff_expr_myeloid_cells$p_val_adj < 0.05 & diff_expr_myeloid_cells$avg_log2FC > 0.1, 
]

# Filter for significant downregulated genes (adj p-value < 0.05 & log2FC < -0.1)
significant_downregulated_genes <- diff_expr_myeloid_cells[
  diff_expr_myeloid_cells$p_val_adj < 0.05 & diff_expr_myeloid_cells$avg_log2FC < -0.1, 
]

# Extract gene symbols for enrichment analysis
upregulated_genes <- rownames(significant_upregulated_genes)
downregulated_genes <- rownames(significant_downregulated_genes)

###-------------------#
### Pathway Enrichment Analysis (KEGG & GO)
###-------------------#

# Helper: Map gene symbols to Entrez IDs
map_genes_to_entrez <- function(gene_symbols) {
  mapped <- bitr(gene_symbols, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Mm.eg.db)
  return(mapped$ENTREZID)
}

# Initialize results list
pathway_results <- list()

# Example: run enrichment for upregulated genes (can repeat for downregulated)
cluster_id <- "Upregulated"  # or any label you choose
cluster_genes <- upregulated_genes

# Map gene symbols to Entrez IDs
cluster_genes_entrez <- map_genes_to_entrez(cluster_genes)
cluster_genes_entrez <- cluster_genes_entrez[!is.na(cluster_genes_entrez)]

# KEGG enrichment analysis
kegg_enrichment <- enrichKEGG(
  gene = cluster_genes_entrez,
  organism = "mmu", 
  pvalueCutoff = 0.05
)
pathway_results[[cluster_id]] <- kegg_enrichment

# GO enrichment analysis (Biological Process)
go_enrichment <- enrichGO(
  gene = cluster_genes_entrez,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05
)
pathway_results[[paste0(cluster_id, "_GO")]] <- go_enrichment

## Visualization
###-------------------#

# Visualize KEGG results
for (cluster_id in names(pathway_results)) {
  if (inherits(pathway_results[[cluster_id]], "enrichResult")) {
    message("Cluster ", cluster_id, " KEGG Pathways:")
    print(dotplot(pathway_results[[cluster_id]], showCategory = 10))
  }
}

# Visualize GO results
for (cluster_id in names(pathway_results)) {
  go_key <- paste0(cluster_id, "_GO")
  if (inherits(pathway_results[[go_key]], "enrichResult")) {
    message("Cluster ", cluster_id, " GO Terms:")
    print(dotplot(pathway_results[[go_key]], showCategory = 10))
  }
}

###-------------------#
### Save Enrichment Results
###-------------------#

# Save KEGG results
for (cluster_id in names(pathway_results)) {
  if (inherits(pathway_results[[cluster_id]], "enrichResult")) {
    write.xlsx(as.data.frame(pathway_results[[cluster_id]]), 
               file = paste0("KEGG_Enrichment_Cluster_", cluster_id, "_Mouse.xlsx"))
  }
}

# Save GO results
for (cluster_id in names(pathway_results)) {
  go_key <- paste0(cluster_id, "_GO")
  if (inherits(pathway_results[[go_key]], "enrichResult")) {
    write.xlsx(as.data.frame(pathway_results[[go_key]]), 
               file = paste0("GO_Enrichment_Cluster_", cluster_id, "_Mouse.xlsx"))
  }
}


###Feature plots of biased genes

FeaturePlot(myeloid_cells, features = "C1qa") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Adgre1") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "C1qc") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Serpinh1") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )


FeaturePlot(myeloid_cells, features = "Cd74") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Mrc1") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )


FeaturePlot(myeloid_cells, features = "Arg1") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )


FeaturePlot(myeloid_cells, features = "H2-Ab1") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Itgax") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Itgam") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )




FeaturePlot(myeloid_cells, features = "Cd14") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Ctsk") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Cd68") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Gatm") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Gatm") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Cpd") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Cxcl9") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )
FeaturePlot(myeloid_cells, features = "Fcgr4") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
  )

FeaturePlot(myeloid_cells, features = "Fgl2") +
  theme_minimal() +
  theme(
    axis.line = element_line(), 
    axis.text = element_text(), 
    axis.title = element_text(),
    panel.grid = element_blank()  # Removes grid lines
