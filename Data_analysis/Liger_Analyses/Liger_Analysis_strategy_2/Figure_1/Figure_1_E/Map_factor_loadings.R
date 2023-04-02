source("/home/yasir/Documents/comparative_study_of_two_species/Functions/Load_libraries.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/load_rdata.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/map_factor_loadings.R")

load("/home/yasir/Documents/comparative_study_of_two_species/Strategy_2/seurat_object_of_K_50.RData")

map_factor_loadings(seuratObject = integrated.data, store_dir = getwd(), factor_to_plot = "inmf", reduction_name = "umap")
