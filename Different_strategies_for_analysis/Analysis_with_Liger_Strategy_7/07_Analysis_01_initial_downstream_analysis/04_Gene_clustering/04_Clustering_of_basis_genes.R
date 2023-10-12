# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_marker_identification", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_matrix_manipulation", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

###############################
###
# Convert the liger object to a seurat object to perform downstream analysis
###

storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_7/Analyses_on_K_50/On_coefficient"

cluster_basis_genes(seurat_object_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_7/Analyses_on_K_50/On_coefficient/Seurat_objects/seurat_object_with_all_factorized_K.RData", 
                    reduction_name_pattern = "^inmf",
                    store_data = TRUE,
                    store_dir = storing_dir,
                    store_folder = "Basis_objects")
