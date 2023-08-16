source("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions/Load_libraries.R")
source("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions/load_rdata.R")
source("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions/optimize_to_new_K.R")

###
# Optimize to a new K from previously factorized K (Higher than the new K) 
###

storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2"

optimize_to_new_k(optimized_liger_object = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/1_Integration_script/integrated_wt_species_liger.RData", 
                  choice_of_new_K = c(10:49),
                  store_dir = storing_dir, 
                  store_folder = "Liger_objects")
