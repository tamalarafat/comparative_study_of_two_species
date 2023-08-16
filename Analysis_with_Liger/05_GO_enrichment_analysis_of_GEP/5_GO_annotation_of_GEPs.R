# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)


###
# Perform GO annotation of the genes in each GEP for different factorization
###

Thaliana_genes =  read.delim("/netscratch/dep_tsiantis/grp_laurent/tamal/2021/Final_Analyses/Input_Files/ATgenes.txt", header = FALSE, col.names = "genes")
Thaliana_genes <- as.character(Thaliana_genes$genes)

storing_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2"

GO_annotation_of_GEPs(gene_clusters_file = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/On_coefficient/Basis_objects/Basis.RData", 
                      ATgenes = Thaliana_genes, 
                      store_dir = storing_dir, 
                      store_folder = "GO_annotation_of_GEPs")

