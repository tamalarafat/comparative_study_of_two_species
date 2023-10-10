# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_2/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/On_coefficient/Analysis_objects/Seurat_object_with_CC/seurat_object_of_K_50.RData")

Idents(integrated.data) <- "RNA_snn_res.0.2"

# Cluster 17 cells
cluster_cells = WhichCells(integrated.data, idents = "12")

# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmf"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

inmf_mat = inmf_mat[cluster_cells, ]

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H$GEP <- apply(df_H, 1, function(x) parse_number(colnames(df_H)[which.max(x)]))

df_H$cell_id = rownames(df_H)
df_H$cell_id = factor(df_H$cell_id, levels = df_H$cell_id)

df_H <-
  df_H %>% dplyr::select(-GEP) %>% melt(
    id.vars = "cell_id",
    variable.name = "GEPs",
    value.name = "Expression"
  )

axis_text = as.character(seq(1, 50, 2))

axis_text_face = rep("plain", length(axis_text))

axis_text_face[which(axis_text == "7")] = "bold"

p <-
  ggplot(df_H, aes(x = cell_id, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion of\nGEP usage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("Cluster 12 cells") +
  scale_y_discrete(breaks = str_c("GEP_", seq(1, 50, 2)), labels = axis_text) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black", face = axis_text_face),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, colour = "black")) + NoLegend()

ggsave(filename = "Figure_2F_RCO_cluster.png", plot = p, width = 12, height = 18, dpi = 400)

######
# 2nd heatmap without sorting
######
# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmf"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H = df_H[cluster_cells, ]

df_H$cell_id = rownames(df_H)
df_H$cell_id = factor(df_H$cell_id, levels = df_H$cell_id)

df_H <-
  df_H %>% melt(
    id.vars = "cell_id",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(df_H, aes(x = cell_id, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion of\nGEP usage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("Cluster 12 cells") +
  scale_y_discrete(breaks = str_c("GEP_", seq(1, 50, 2)), labels = axis_text) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black", face = axis_text_face),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, colour = "black")) + NoLegend()

ggsave(filename = "Figure_2F_RCO_cluster_without_sorting.png", plot = p, width = 12, height = 18, dpi = 400)


######
# 3rd heatmap - rco cells sorted
######

# Check the RCO expressing cells
Cells_with_expression_detection = names(GetAssayData(integrated.data, assay = "RNA", slot = "counts")["AT5G67651", GetAssayData(integrated.data, assay = "RNA", slot = "counts")["AT5G67651", ] != 0])

writeLines(Cells_with_expression_detection, "RCO_expressing_cells.txt")

cluster_expressing_cells = intersect(Cells_with_expression_detection, cluster_cells)

writeLines(cluster_expressing_cells, "RCO_expressing_cells_cluster_12.txt")

non_expressing_cells = setdiff(cluster_cells, cluster_expressing_cells)

cells_ordered = c(cluster_expressing_cells, non_expressing_cells)

# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmf"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H = df_H[cells_ordered, ]

df_H$cell_id = rownames(df_H)
df_H$cell_id = factor(df_H$cell_id, levels = df_H$cell_id)

df_H <-
  df_H %>% melt(
    id.vars = "cell_id",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(df_H, aes(x = cell_id, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion of\nGEP usage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("Cluster 12 cells") +
  scale_y_discrete(breaks = str_c("GEP_", seq(1, 50, 2)), labels = axis_text) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black", face = axis_text_face),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, colour = "black")) + NoLegend()

ggsave(filename = "Figure_2F_RCO_cluster_RCO_cells_sorted.png", plot = p, width = 12, height = 18, dpi = 400)


######
# 4th heatmap - only rco expressing cells
######
# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H = df_H[cluster_expressing_cells, ]

df_H$cell_id = rownames(df_H)
df_H$cell_id = factor(df_H$cell_id, levels = df_H$cell_id)

df_H <-
  df_H %>% melt(
    id.vars = "cell_id",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(df_H, aes(x = cell_id, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion of\nGEP usage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("Cluster 12 - RCO expressing cells") +
  scale_y_discrete(breaks = str_c("GEP_", seq(1, 50, 2)), labels = axis_text) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black", face = axis_text_face),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, colour = "black")) + NoLegend()

ggsave(filename = "Figure_2F_RCO_cluster_only_RCO_cells.png", plot = p, width = 12, height = 18, dpi = 400)


######
# 3rd heatmap - rco cells sorted
######

# Check the RCO expressing cells
Cells_with_expression_detection = names(GetAssayData(integrated.data, assay = "RNA", slot = "counts")["AT5G67651", GetAssayData(integrated.data, assay = "RNA", slot = "counts")["AT5G67651", ] != 0])

cluster_expressing_cells = intersect(Cells_with_expression_detection, cluster_cells)

rest_expressing_cells = setdiff(Cells_with_expression_detection, cluster_expressing_cells)

cells_ordered = c(cluster_expressing_cells, rest_expressing_cells)

# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmf"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H = df_H[cells_ordered, ]

df_H$cell_id = rownames(df_H)
df_H$cell_id = factor(df_H$cell_id, levels = df_H$cell_id)

df_H <-
  df_H %>% melt(
    id.vars = "cell_id",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(df_H, aes(x = cell_id, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion of\nGEP usage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("All RCO expressing cells") +
  scale_y_discrete(breaks = str_c("GEP_", seq(1, 50, 2)), labels = axis_text) +
  geom_vline(xintercept = length(cluster_expressing_cells)) + 
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black", face = axis_text_face),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, colour = "black")) + NoLegend()

ggsave(filename = "Figure_2F_all_RCO_cells_sorted.png", plot = p, width = 12, height = 18, dpi = 400)
