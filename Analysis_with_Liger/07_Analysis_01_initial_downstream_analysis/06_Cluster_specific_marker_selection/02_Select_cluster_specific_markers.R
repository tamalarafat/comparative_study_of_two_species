# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Projects_Yasir/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Markers from cluster DEGs
# Dir - containing the DEG files
DEG_dir = "/home/yasir/Documents/Projects_Yasir/comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

# Storing directory
storing_dir = "/home/yasir/Documents/Projects_Yasir/comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs"

specific_markers_list = specific_marker_finder(DEG_file_dir = DEG_dir, 
                                               pct_detection = 0.1, 
                                               pct_diff = 0.3, 
                                               store_dir = storing_dir, 
                                               store.marker.file = TRUE)
