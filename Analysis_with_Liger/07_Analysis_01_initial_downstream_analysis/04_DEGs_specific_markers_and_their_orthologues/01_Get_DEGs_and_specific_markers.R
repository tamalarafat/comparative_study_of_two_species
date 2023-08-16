# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_2/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_2/Analysis_objects/Seurat_object_with_CC/seurat_object_of_K_50.RData")

Idents(integrated.data) <- "RNA_snn_res.0.2"

# Load the gene description file
ATH = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_2/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$AT_ID

# Dir - containing the DEG files
DEG_dir = "/home/yasir/Documents/Thesis_PhD/Chapter_2/Liger_analysis_with_CC/Differentially_expressed_genes/Conserved_marker_grouped_by_Species/Conserved_markers_DEtest_wilcox"

markers_list = get_conserved_DEGs_details_thaliana(DEG_file_dir = DEG_dir, gene_description_file = ATH, return_markers_list = TRUE)

save(markers_list, file = "DEGs_and_Markers.RData")
