projects_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/"

# projects_dir = "~/Documents/Projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_1 <- list.files(paste0(projects_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE)
sapply(scripts_list_1, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(projects_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE)
sapply(scripts_list_2, source, .GlobalEnv)

# Dir - containing the DEG files
DEG_dir = paste0(projects_dir, "comparative_study_of_two_species/Different_strategies_for_analysis/Analysis_with_Liger_Strategy_7/07_Analysis_01_initial_downstream_analysis/Analysis_outputs/Differentially_expressed_genes/Conserved_DEGs_grouped_by_Species/Conserved_markers_DEtest_wilcox")

# Storing directory
storing_dir = paste0(projects_dir, "comparative_study_of_two_species/Different_strategies_for_analysis/Analysis_with_Liger_Strategy_7/07_Analysis_01_initial_downstream_analysis/Analysis_outputs")

specific_markers_list = specific_conserved_marker_finder(DEG_file = DEG_dir, 
                                                         max_pct2_detection = 0.1, 
                                                         include_pct_diff = TRUE, 
                                                         pct_diff = 0.3, 
                                                         store_outputs = TRUE, 
                                                         store_dir = storing_dir, 
                                                         store_folder = "Cluster_specific_conserved_markers")
