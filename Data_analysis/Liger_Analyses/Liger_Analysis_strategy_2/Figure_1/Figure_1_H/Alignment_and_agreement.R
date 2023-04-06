# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/Liger_objects/Liger_object_K_50.RData")

a = calcAlignment(Liger_object)
b = calcAgreement(Liger_object, ndims = ncol(Liger_object@H.norm))

# Creating an empty list to store the grouping results
df_list <- list()


df_list[[1]] = data.frame("Alignment" = a, "Agreement" = b)
names(df_list)[1] <- "K_50"

df_list[[2]] = data.frame("Alignment" = c, "Agreement" = d)
names(df_list)[2] <- "K_49"

temp <- do.call(rbind.data.frame, df_list)

# Liger objects directory
object_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/Liger_objects"

list_objects = str_sort(list.files(path = object_dir, pattern = "Liger_object_K"), numeric = TRUE)

file_name = str_c(object_dir, list_objects[1])

ncol(Liger_object@H.norm)

