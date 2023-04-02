source("/home/yasir/Documents/comparative_study_of_two_species/Functions/Load_libraries.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/load_rdata.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/order_factorized_matrix.R")

# # Load all the functions stored in scripts from the folder housing the scripts
# scripts_list <- list.files("/home/yasir/Documents/comparative_study_of_two_species/Functions", pattern = "*.R$", full.names = TRUE) 
# sapply(scripts_list, source, .GlobalEnv)

# Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A")

load("/home/yasir/Documents/comparative_study_of_two_species/Strategy_2/seurat_object_of_K_50.RData")

p <-  FeaturePlot(integrated.data, features = "AT1G19790", reduction = "umap", dims = c(1, 2), pt.size = 1, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[8])) + 
  ggtitle(str_c("AT1G19790", " - " , "SRS7")) + 
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.ticks.length = unit(.30, "cm"), 
    axis.text = element_text(size = 18, face = "bold"),
    title = element_text(size = 24, face = "bold"),
    legend.key.size = unit(2,"line"),
    legend.key = element_rect(size = 20),
    legend.text = element_text(size = 18, face = "bold"))
ggsave(filename = "Figure_1G_1.png", plot = p, width = 14, height = 14, dpi = 300)

p <-  FeaturePlot(integrated.data, features = "AT1G62360", reduction = "umap", dims = c(1, 2), pt.size = 1, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[8])) + 
  ggtitle(str_c("AT1G62360", " - " , "STM")) + 
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.ticks.length = unit(.30, "cm"), 
    axis.text = element_text(size = 18, face = "bold"),
    title = element_text(size = 24, face = "bold"),
    legend.key.size = unit(2,"line"),
    legend.key = element_rect(size = 20),
    legend.text = element_text(size = 18, face = "bold"))
ggsave(filename = "Figure_1G_2.png", plot = p, width = 14, height = 14, dpi = 300)

