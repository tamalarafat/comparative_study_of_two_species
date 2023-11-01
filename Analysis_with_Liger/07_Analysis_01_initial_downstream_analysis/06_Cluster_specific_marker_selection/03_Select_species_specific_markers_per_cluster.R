project_dir = "/home/yasir/Documents/Projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# Markers from cluster DEGs
# Dir - containing the DEG files
DEG_dir = paste0(project_dir, "comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs/Differentially_expressed_genes/DEGs_between_species_of_a_cluster/DEG_between_Species_DEtest_wilcox")

# Storing directory
storing_dir = paste0(project_dir, "comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs")

temp_files = str_sort(str_c(list.files(DEG_dir, full.names = TRUE), "/"), numeric = TRUE)

if (!dir.exists(str_c(storing_dir, "/", "Specific_markers_between_species_per_cluster"))){
  dir.create(str_c(storing_dir, "/", "Specific_markers_between_species_per_cluster"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Create a folder to store the markers file
temp_holder = str_c(storing_dir, "/", "Specific_markers_between_species_per_cluster", "/")


for (i in c(1:length(temp_files))){
  
  if (!dir.exists(str_c(temp_holder, "/", basename(temp_files[i])))){
    dir.create(str_c(temp_holder, "/", basename(temp_files[i])), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(temp_holder, "/", basename(temp_files[i]), "/")
  
  temp_file_path = list.files(temp_files[i], full.names = TRUE, pattern = "Cluster")
  
  if (length(temp_file_path) != 0) {
    for (j in c(1:length(temp_file_path))) {
      cluster_deg = loadRData(temp_file_path[j])
      
      specific_markers = specific_marker_finder(DEG_file_dir = cluster_deg, max_pct2_detection = 0.1, pct_diff = 0.3, include_pct_diff = TRUE, store_outputs = FALSE)
      
      write.csv(specific_markers, file = str_c(temp_dir, sub(pattern = "\\.[^.]+$", replacement = "", x = basename(temp_file_path[j])), "_specific_markers.csv"), row.names = TRUE)
    }
  }
  
  else {
    next
  }
}

