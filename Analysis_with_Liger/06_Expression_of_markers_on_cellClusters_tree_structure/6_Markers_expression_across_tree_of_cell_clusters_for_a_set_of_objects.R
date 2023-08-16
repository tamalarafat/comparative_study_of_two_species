# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)


# Load the known cell type markers file
known_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Thaliana_known_cell_type_markers.csv")

# Storing directory
storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/On_coefficient"

clus_tree_series_of_seurat_objects(seurat_object_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/On_coefficient/Seurat_objects", 
                                   file_name_pattern = "seurat_object_of_", 
                                   marker_file = known_markers, 
                                   query_pattern = "RNA_snn_res.", 
                                   gene_ID_column = "AT_ID", 
                                   gene_name_column = "Gene_Name", 
                                   store_dir = storing_dir)

