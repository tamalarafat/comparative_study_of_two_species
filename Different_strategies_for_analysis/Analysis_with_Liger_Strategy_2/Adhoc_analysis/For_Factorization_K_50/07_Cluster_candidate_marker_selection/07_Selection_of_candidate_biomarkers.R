# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_marker_identification", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_matrix_manipulation", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# Load the Seurat object
load("/netscratch/dep_tsiantis/common/Yasir/Comparative_species_analysis/Liger_analysis/Liger_analysis_without_cell_cycle_effect/Analysis_objects/seurat_object_of_K_50_without_CC.RData")

# Directory path to store the outputs
storing_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/comparative_study_of_two_species/Analysis_with_Liger/08_Analysis_02_controlling_for_cell_cycle_effect/Analysis_outputs"

# Path to the DEG files - differentially genes using the "Findmarkers" function (Seurat) per cluster
DEG_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/comparative_study_of_two_species/Analysis_with_Liger/08_Analysis_02_controlling_for_cell_cycle_effect/Analysis_outputs/Differentially_expressed_genes/Conserved_marker_grouped_by_Species/Conserved_markers_DEtest_wilcox"

# Load the TFs list 
at_TFs = read.delim("/home/ytamal2/Documents/2023/PhD_projects_Yasir/comparative_study_of_two_species/Analysis_with_Liger/Analysis_objects/Annotation_files/TFs_list/Ath_TF_list.txt", header = 1)
at_TFs = unique(at_TFs$Gene_ID)

# Cluster IDs of the active clustering solution
cluster_IDs = str_sort(levels(integrated.data), numeric = TRUE)

# Empty table to store the degs table
deg_table_list = list()

# From the deg files directory get the file names and load the degs tables
if (is.character(DEG_dir)) {
  
  # get the file names
  file_names = str_sort(list.files(DEG_dir, pattern = "Cluster"), numeric = TRUE)
  
  for (i in c(1:length(file_names))){
    
    marker_file = loadRData(str_c(DEG_dir, "/", file_names[i]))
    
    deg_table_list[[i]] <- marker_file
    
    if (grepl("[[:digit:]]", file_names[i]) != TRUE){
      names(deg_table_list)[i] <- sub("(.*)\\.RData$", "\\1", file_names[i])
    }
    
    else {
      names(deg_table_list)[i] <- str_c("cluster_", parse_number(file_names[i]))
    }
  }
}

# For each cluster ID define the GEP IDs from which the candidate markers will be selected
GEP_IDs = list(cluster_0 = c(20, 35), 
               cluster_1 = 35, 
               cluster_2 = 23, 
               cluster_3 = 43, 
               cluster_4 = c(20, 35), 
               cluster_5 = 25, 
               cluster_6 = c(24, 43), 
               cluster_7 = 35, 
               cluster_8 = 19, 
               cluster_9 = 2, 
               cluster_10 = 11, 
               cluster_11 = 35, 
               cluster_12 = 7, 
               cluster_13 = 17, 
               cluster_14 = 1, 
               cluster_15 = 46, 
               cluster_16 = 15, 
               cluster_17 = 5, 
               cluster_18 = 9, 
               cluster_19 = 33)


# Creating an empty list to which GEP IDs and DEGs table for each cluster will be stored. 
# This list will be served as query input to identify candidate biomarkers using the GEP IDs and degs for each cluster
query_list = list()

for (i in c(1:length(cluster_IDs))) {
  query_list[[i]] = list(cluster_ID = cluster_IDs[i], # First element is the cluster ID
                        GEP_IDs = GEP_IDs[[grep(pattern = str_c("^", cluster_IDs[i], "$", sep = ""), parse_number(names(GEP_IDs)))]], # Second element is the GEP IDs
                        deg_table = deg_table_list[[grep(pattern = str_c("^", cluster_IDs[i], "$", sep = ""), parse_number(names(deg_table_list)))]]) # 3rd element is the degs table
  
  names(query_list)[i] = str_c("cluster_", cluster_IDs[i])
}


# Now, iterate over the list item to identify candidate biomarkers for each cluster using the query inputs
for (i in c(1:length(query_list))){
  query_inputs = query_list[[i]]
  
  cluster_candidates = candidate_markers_GEP_and_DEGs(seurat_object = integrated.data, 
                                                      GEP_IDs = query_inputs$GEP_IDs, 
                                                      cluster_ID = query_inputs$cluster_ID, 
                                                      DEG_file = query_inputs$deg_table, 
                                                      find_candidates = 10, 
                                                      reduction_name = "inmf", 
                                                      max_pct2_detection = 0.1, 
                                                      pct_diff = 0.3, 
                                                      include_pct_diff = TRUE, 
                                                      specify_gene_ID = at_TFs, 
                                                      incorporate_column_name = "TFs", 
                                                      store_outputs = TRUE, 
                                                      store_dir = storing_dir)
}

