# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/On_coefficient/Seurat_objects/seurat_object_of_K_50.RData")

############################ Set the resolution parameter ############################
Idents(integrated.data) <- integrated.data$RNA_snn_res.0.2
Idents(integrated.data) <- factor(Idents(integrated.data), levels = seq(0, length(levels(Idents(integrated.data))) - 1))
paste("What is the active ident?", paste(levels(integrated.data@active.ident), collapse = ", "))
###########################################################################

# Directory - Storing results
res_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/Figure_1/Select_UMAP/sinmf"


# Lets get the metadata file
md = integrated.data@meta.data


# UMAP - 2 components reduction

n_neighbors = c(20, 30, 40, 50)
n_seed = c(1:50)

for (i in c(1:length(n_neighbors))){
  for (j in c(1:length(n_seed))){
    integrated.data <- RunUMAP(integrated.data, reduction = "sinmf", dims = 1:50, n.components = 2, n.neighbors = n_neighbors[i], seed.use = n_seed[j])
    
    RDimension_plot(seuratObject = integrated.data, store_dir = res_dir, store_folder = "Seed_and_NN",dimension_reduction_name = "umap", figure_name_pref = str_c("_nn", n_neighbors[i], "_seed", n_seed[j]))
    }
}


