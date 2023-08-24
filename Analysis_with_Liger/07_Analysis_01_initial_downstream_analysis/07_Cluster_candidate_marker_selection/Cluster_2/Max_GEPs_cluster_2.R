# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_2/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_2/Analysis_objects/Seurat_object_with_CC/seurat_object_of_K_50.RData")

Idents(integrated.data) <- "RNA_snn_res.0.2"

target_ID = "2"

# Load the gene description file
ATH = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_2/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$AT_ID

cluster_cells = WhichCells(integrated.data, idents = target_ID)

# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmf"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

inmf_mat = inmf_mat[cluster_cells, ]

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H$GEP <- apply(df_H, 1, function(x) parse_number(colnames(df_H)[which.max(x)]))

table(df_H$GEP)

cluster_GEP_usage = data.frame(table(df_H$GEP))
colnames(cluster_GEP_usage) = c("GEP", "cell_count")

save(cluster_GEP_usage, file = str_c("cluster_", target_ID, "_GEP_usage.RData"))


