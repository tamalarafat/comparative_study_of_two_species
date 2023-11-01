project_dir = "/home/yasir/Documents/Projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# Markers from cluster DEGs
# Dir - containing the DEG files
DEG_dir = paste0(project_dir, "comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs/Differentially_expressed_genes/DEGs_for_a_species_plus_cluster/Species_and_cluster_DEtest_wilcox")

# Storing directory
storing_dir = paste0(project_dir, "comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs")

temp_files = str_c(list.files(DEG_dir, pattern = "Cluster", full.names = TRUE))

if (!dir.exists(str_c(storing_dir, "/", "Species_and_cluster_specific_markers"))){
  dir.create(str_c(storing_dir, "/", "Species_and_cluster_specific_markers"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Create a folder to store the markers file
temp_dir = str_c(storing_dir, "/", "Species_and_cluster_specific_markers", "/")


for (i in c(1:length(temp_files))){
    cluster_deg = loadRData(temp_files[i])
    
    specific_markers = specific_marker_finder(DEG_file_dir = cluster_deg, max_pct2_detection = 0.1, pct_diff = 0.3, include_pct_diff = TRUE, store_outputs = FALSE)
    
    write.csv(specific_markers, file = str_c(temp_dir, sub(pattern = "\\.[^.]+$", replacement = "", x = basename(temp_files[i])), "_specific_markers.csv"), row.names = TRUE)
  }


