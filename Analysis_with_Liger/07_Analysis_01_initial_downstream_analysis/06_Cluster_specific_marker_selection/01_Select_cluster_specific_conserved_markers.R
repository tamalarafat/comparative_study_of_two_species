# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Projects_Yasir/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Markers from cluster DEGs
# Dir - containing the DEG files
DEG_dir = "/home/yasir/Documents/Projects_Yasir/comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs/Differentially_expressed_genes/Conserved_DEGs_grouped_by_Species/Conserved_markers_DEtest_wilcox"

DEG_files = str_sort(list.files(DEG_dir, pattern = "Cluster_"), numeric = TRUE)

Cluster_DEGs = lapply(str_c(DEG_dir, "/", DEG_files), loadRData)

temp_names = names(Cluster_DEGs)

# Storing directory
storing_dir = "/home/yasir/Documents/Projects_Yasir/comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs"

specific_markers_list = specific_conserved_marker_finder(DEG_file = Cluster_DEGs, 
                                                         max.pct2.detection = 0.1, 
                                                         include.pct.difference = TRUE, 
                                                         pct.difference = 0.3, 
                                                         store.marker.file = TRUE, 
                                                         store_dir = storing_dir, 
                                                         store_folder = "Cluster_specific_conserved_markers")
