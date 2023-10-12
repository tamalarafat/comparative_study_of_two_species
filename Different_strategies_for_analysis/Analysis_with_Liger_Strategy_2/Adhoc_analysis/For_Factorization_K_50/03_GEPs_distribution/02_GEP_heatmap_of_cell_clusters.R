# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_marker_identification", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_matrix_manipulation", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_5/On_coefficient/Seurat_objects/seurat_object_of_K_50.RData")

Idents(integrated.data) <- "RNA_snn_res.0.2"

cluster.ids = levels(Idents(integrated.data))

for (i in c(1:length(cluster.ids))){
  heatmap_cells_coefficient(seuratObject = integrated.data, store_dir = getwd(), store_folder = "Heatmap_of_clusters", cell_ids = WhichCells(integrated.data, idents = cluster.ids[i]), figureName = str_c("Cluster_", cluster.ids[i]))
}

