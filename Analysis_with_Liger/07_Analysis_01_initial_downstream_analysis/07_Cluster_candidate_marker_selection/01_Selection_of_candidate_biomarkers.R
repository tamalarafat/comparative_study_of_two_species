# projects_dir = "~/Documents/2023/PhD_projects_Yasir/"

projects_dir = "~/Documents/Projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files(paste0(projects_dir, "scExplorer/Functions"), pattern = "*.R$", full.names = TRUE)

# Split the vector into two chunks
num_chunks <- 4
chunks <- split(scripts_list, 1:num_chunks)

# Loop through each chunk and source the scripts
for (i in c(1:num_chunks)) {
  lapply(chunks[[i]], source, .GlobalEnv)
}

source(paste0(projects_dir, "scExplorer/Functions/specific_marker_finder.R"))

# sapply(scripts_list, source, .GlobalEnv)

# # Load the Seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/On_coefficient/Analysis_objects/Seurat_object_with_CC/seurat_object_of_K_50.RData")

# Directory path to store the outputs
storing_dir = paste0(projects_dir, "comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs")

# Path to the DEG files - differentially genes using the "Findmarkers" function (Seurat) per cluster
DEG_dir = paste0(projects_dir, "comparative_study_of_two_species/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs/Differentially_expressed_genes/Conserved_DEGs_grouped_by_Species/Conserved_markers_DEtest_wilcox")

# Load the TFs list 
at_TFs = read.delim(paste0(projects_dir, "comparative_study_of_two_species/Analysis_with_Liger/Analysis_objects/Annotation_files/TFs_list/Ath_TF_list.txt"), header = 1)
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

