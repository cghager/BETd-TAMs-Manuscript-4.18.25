# 1. Load Required Libraries
########################################

# Core packages
library(dplyr)
library(ggplot2)
library(tibble)
library(stringr)
library(openxlsx)
library(writexl)
library(tidyverse)
library(future)
library(future.apply)
library(parallel)
library(ggrepel)
library(ggpubr)
library(reshape2)

# Seurat & Analysis Tools
library(Seurat)
library(SeuratObject)
library(glmGamPoi)
library(cluster)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(Nebulosa)

# Visualization & Networks
library(patchwork)
library(pheatmap)
library(circlize)
library(networkD3)
library(HiveR)
library(ggalluvial)

# Cell-Cell Communication
library(CellChat)

# NicheNet (for signaling inference)
if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("saeyslab/nichenetr")

########################################
# 2. Subset and Prepare Seurat Objects
########################################

# Subset Myeloid and Tumor cells
myeloid_cells <- subset(obj, idents = "Myeloid cells")
tumor_cells <- subset(obj, idents = "Tumor cells")

# Add metadata for cell type
myeloid_cells$cell_type <- "Myeloid"
tumor_cells$cell_type <- "Tumor"

# Merge objects
combined_cells <- merge(myeloid_cells, y = tumor_cells, add.cell.ids = c("Myeloid", "Tumor"))

# Save and reload combined object
saveRDS(combined_cells, file = "combined_myeloid_tumor.rds")
combined_cells <- readRDS("combined_myeloid_tumor.rds")

# Split by treatment
vehicle_cells <- subset(combined_cells, subset = treatment == "Vehicle")
zbc260_cells <- subset(combined_cells, subset = treatment == "ZBC260")

# Save Vehicle cells for CellChat analysis
saveRDS(vehicle_cells, file = "Vehicle_combined_chat_obj.rds")
Vehicle_obj <- readRDS("Vehicle_combined_chat_obj.rds")


# 3. Initialize CellChat Analysis
########################################

gc()
options(future.globals.maxSize = 5 * 1024^3)
future::plan("multisession", workers = 2)

# Create CellChat object
cellchat <- createCellChat(object = Vehicle_obj, group.by = "ident", assay = "SCT")
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

saveRDS(cellchat, "cellchat_combined_Vehicle.rds")
cellchat <- readRDS("cellchat_combined_Vehicle.rds")

# 4. Communication Probability & Network
########################################

cellchat <- computeCommunProb(cellchat, nboot = 5)
saveRDS(cellchat, "cellchat_combined_prob_Vehicle.rds")
cellchat <- readRDS("cellchat_combined_prob_Vehicle.rds")

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

# Circle plots
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction Count")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction Strength")

########################################
# 5. Network Visualization & Data Extraction
########################################

# Interaction table
df.net <- subsetCommunication(cellchat)

# Bubble plot
netVisual_bubble(cellchat, sources.use = 1:4, targets.use = 1:4, remove.isolate = TRUE) +
  theme(legend.text = element_text(face = "bold"))

# 3D to 2D Data Conversion for Top Pathways
top_pathways <- names(sort(apply(cellchat@net$prob, 3, sum), decreasing = TRUE))[1:20]
df <- melt(cellchat@net$prob[, , top_pathways])
colnames(df) <- c("Sender", "Receiver", "Pathway", "Probability")

# Heatmap
ggplot(df, aes(x = Sender, y = Receiver, fill = Probability)) +
  geom_tile() +
  facet_wrap(~Pathway) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()

# Bubble plot
ggplot(df, aes(x = Sender, y = Receiver, size = Probability, color = Pathway)) +
  geom_point(alpha = 0.8) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

########################################
# 6. Ligand-Receptor Interaction Plots
########################################

top_ligands <- names(sort(apply(cellchat@net$prob, 3, sum), decreasing = TRUE)[1:20])
df_ligands <- cellchat@net$prob[,,top_ligands] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sender") %>%
  pivot_longer(cols = -Sender, names_to = "Receiver", values_to = "Interaction") %>%
  filter(!is.na(Interaction) & Interaction > 0)

# Dot plot
ggplot(df_ligands, aes(x = Sender, y = Receiver, size = Interaction, color = Interaction)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  theme_minimal()

########################################
# 7. Chord Diagram
########################################

df_chord <- df_ligands %>% select(Sender, Receiver, Interaction)
circos.clear()
chordDiagram(df_chord, 
             col = colorRamp2(c(0, max(df_chord$Interaction)), c("blue", "red")),
             transparency = 0.5,
             annotationTrack = "grid")

# Add sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_name <- get.cell.meta.data("sector.index")
  if (!is.null(sector_name) && sector_name != "") {
    circos.text(mean(get.cell.meta.data("xlim")), get.cell.meta.data("ylim")[2] - 1.5,
                labels = sector_name, facing = "clockwise", niceFacing = TRUE, cex = 0.8)
  }
}, bg.border = NA)

########################################
# 8. Sankey Diagram (Alluvial)
########################################

ggplot(df, aes(axis1 = Sender, axis2 = Receiver, y = Probability, fill = Pathway)) +
  geom_alluvium() +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  labs(title = "Sankey Diagram of Top 5 Pathways")



##########Repeat same for ZBC260
cellchat_ZBC260 <- subsetData(cellchat_ZBC260)
cellchat_ZBC260 <- identifyOverExpressedGenes(cellchat_ZBC260)
cellchat_ZBC260 <- identifyOverExpressedInteractions(cellchat_ZBC260)
cellchat_ZBC260 <- projectData(cellchat_ZBC260, PPI.mouse)

cellchat_ZBC260 <- computeCommunProb(cellchat_ZBC260, nboot = 5)
cellchat_ZBC260 <- aggregateNet(cellchat_ZBC260)
saveRDS(cellchat_ZBC260, file = "cellchat_combined_prob_ZBC260.rds")




#Graphing
groupSize <- as.numeric(table(cellchat_ZBC260@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_ZBC260@net$count, vertex.weight = groupSize, title.name = "Number of interactions")
netVisual_circle(cellchat_ZBC260@net$weight, vertex.weight = groupSize, title.name = "Interaction weights/strength")



p <- netVisual_bubble(cellchat_ZBC260, sources.use = 1:4, targets.use = 1:4, remove.isolate = TRUE)
p + theme(legend.text = element_text(face = "bold"))



#TOP-LIGAND RECEPTO INTERACTIONS


df_ligands <- reshape2::melt(cellchat@net$prob[,,top_ligands])
colnames(df_ligands) <- c("Sender", "Receiver", "Pathway", "Interaction")
df_ligands <- df_ligands %>% filter(Interaction > 0)


ggplot(df_ligands, aes(x = Sender, y = Receiver, fill = Interaction)) +
  geom_tile() +
  facet_wrap(~Pathway) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()
ggplot(df_ligands, aes(x = Sender, y = Receiver, size = Interaction, color = Interaction)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  theme_minimal()

chordDiagram(df_chord, transparency = 0.5, col = colorRamp2(c(0, max(df_chord$Interaction)), c("blue", "red")))
