# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_marker_identification", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_matrix_manipulation", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# Storing directory
res_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_7/Analyses_on_K_50/On_coefficient"

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_7/Analyses_on_K_50/On_coefficient/Seurat_objects/seurat_object_of_K_50.RData")

Idents(integrated.data) <- "RNA_snn_res.0.2"

between_groups_differentially_expressed_genes(seuratObject = integrated.data, between_group_variable = "Species", store_dir = res_dir, DEGtest = "wilcox")
