# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2023/Beginning_of_a_compendium/a_scRNAseq_analysis_pipeline/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load the Orthologues table - Contains orthologues ID of Arabidopsis thaliana and Cardamine hirsuta.
# The orthologues table has been generated with the TAIR10 Arabidopsis thaliana gene IDs
###

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

# Loading data -  loading of sparse data matrices provided by 10X genomics using Seurat's "Read10X" function
COL_data_5E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_5th_ALL_2/filtered_feature_bc_matrix/")

# Convert the gene ids in the count matrix to desired orthologues gene ids

# The function "prepare_ortho_data" takes a sparse matrix (loaded with the "Read10X" function) as input.
# parameters :: ortho_data - the orthologues table 
# parameters :: ortho_column_name_of_gene_ids - The name of the column in the orthologues table that contains the gene IDs that we want to convert.
# parameters :: ortho_column_name_to_assign - The name of the column in the orthologues table that contains the gene IDs that we want to assign.

# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(COL_data_5E)

# extracting the Cardamine IDs that are present in orthologues table 
thaliana_ortho_genes = as.character(ortho_table$A.thaliana.TAIR10)

# not all the thaliana ids are present in the ortho data - 51 thaliana genes are missing in the thaliana data
thaliana_ortho_genes = intersect(thaliana_genes, thaliana_ortho_genes)

# Let's subset the data with the ortho genes
Ortho_COL_REP_4 <- COL_data_5E[thaliana_ortho_genes, ]
# save(Ortho_COL_REP_4, file = "Arabidopsis_REP_4_Ortho_data.RData")
