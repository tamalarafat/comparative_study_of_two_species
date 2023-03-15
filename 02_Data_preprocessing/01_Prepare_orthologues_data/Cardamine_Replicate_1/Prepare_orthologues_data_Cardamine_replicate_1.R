# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2023/Beginning_of_a_compendium/a_scRNAseq_analysis_pipeline/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load the Orthologues table - Contains orthologues ID of Arabidopsis thaliana and Cardamine hirsuta.
# The orthologues table has been generated with the TAIR10 Arabidopsis thaliana gene IDs
###

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

# Loading data -  loading of sparse data matrices provided by 10X genomics using Seurat's "Read10X" function
OX_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_OX_RNA_1ST_2_Newest/filtered_feature_bc_matrix/")

# Convert the gene ids in the count matrix to desired orthologues gene ids

# The function "prepare_ortho_data" takes a sparse matrix (loaded with the "Read10X" function) as input.
# parameters :: ortho_data - the orthologues table 
# parameters :: ortho_column_name_of_gene_ids - The name of the column in the orthologues table that contains the gene IDs that we want to convert.
# parameters :: ortho_column_name_to_assign - The name of the column in the orthologues table that contains the gene IDs that we want to assign.

Ortho_OX_REP_1 = prepare_ortho_data(input_data = OX_data_1E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

# save(Ortho_OX_REP_1, file = "Cardamine_REP_1_Ortho_data.RData")