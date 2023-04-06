# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load data tables
###

###
# Protoplasting induced genes
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Ox_Co0_leaf_protoplast_v12_final_August_2022_ortho.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

###
# Orthologues table
###

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

###
# WT C. hirsuta
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
OX_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_OX_RNA_1ST_2_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 2nd Experiment (leaf 6 and 7)
OX_data_2E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_2nd_ALL_2_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 3rd Experiment (leaf 5 and 6)
OX_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_3rd_ALL_3000_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 7th Experiment (leaf 6 and 7)
OX_data_7E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_7th_ALL_2_Newest/filtered_feature_bc_matrix/")

# Convert the gene ids in the data table to ortho gene ids
OX_DF_1 = prepare_ortho_data(input_data = OX_data_1E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_2 = prepare_ortho_data(input_data = OX_data_2E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_3 = prepare_ortho_data(input_data = OX_data_3E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

OX_DF_7 = prepare_ortho_data(input_data = OX_data_7E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

###
# WT A. thaliana
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
COL_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_1ST_2/filtered_feature_bc_matrix/")

# Load data - WT OX 2nd Experiment (leaf 6 and 7)
COL_data_2E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col_RNA_2nd_ALL_2/filtered_feature_bc_matrix/")

# Load data - WT OX 3rd Experiment (leaf 5 and 6)
COL_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_3rd_ALL/filtered_feature_bc_matrix/")

# Load data - WT OX 7th Experiment (leaf 6 and 7)
COL_data_5E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_5th_ALL_2/filtered_feature_bc_matrix/")

# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(COL_data_1E)

# extracting the Cardamine IDs that are present in orthologues table 
thaliana_ortho_genes = as.character(ortho_table$A.thaliana.TAIR10)

# not all the thaliana ids are present in the ortho data - 51 thaliana genes are missing in the thaliana data
thaliana_ortho_genes = intersect(thaliana_genes, thaliana_ortho_genes)

# Let's subset the data with the ortho genes
COL_data_1E <- COL_data_1E[thaliana_ortho_genes, ]
COL_data_2E <- COL_data_2E[thaliana_ortho_genes, ]
COL_data_3E <- COL_data_3E[thaliana_ortho_genes, ]
COL_data_5E <- COL_data_5E[thaliana_ortho_genes, ]

# remove the missing genes from the data
OX_DF_1 <- OX_DF_1[thaliana_ortho_genes, ]
OX_DF_2 <- OX_DF_2[thaliana_ortho_genes, ]
OX_DF_3 <- OX_DF_3[thaliana_ortho_genes, ]
OX_DF_7 <- OX_DF_7[thaliana_ortho_genes, ]

# All gene IDs - all datasets
ortho_genes = rownames(COL_data_1E)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(ortho_genes, PP_genes)


##### Remove the protoplasting induced genes
OX_DF_1 <- OX_DF_1[genes_to_keep, ]
OX_DF_2 <- OX_DF_2[genes_to_keep, ]
OX_DF_3 <- OX_DF_3[genes_to_keep, ]
OX_DF_7 <- OX_DF_7[genes_to_keep, ]

COL_data_1E <- COL_data_1E[genes_to_keep, ]
COL_data_2E <- COL_data_2E[genes_to_keep, ]
COL_data_3E <- COL_data_3E[genes_to_keep, ]
COL_data_5E <- COL_data_5E[genes_to_keep, ]

###
# OX - 1 E
###

# First replicate - OX 1E - total cells 6640; filter out genes that are not detected in at least 13 cells
OX_1E <- CreateSeuratObject(counts = OX_DF_1, project = "OX_1E", min.features = 200)

# Add metadata information to the seurat object
OX_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-1", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_1E <- subset(OX_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_1E[["percent.mt"]] <- PercentageFeatureSet(OX_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_1E[["percent.pt"]] <- PercentageFeatureSet(OX_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_1E <- subset(OX_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_1E <- NormalizeData(OX_1E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_1E <- FindVariableFeatures(OX_1E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W1 <- GetAssayData(OX_1E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W1) <- paste("O1", colnames(OX_W1), sep = "_")


###
# OX - 2 E
###

# First replicate - OX 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
OX_2E <- CreateSeuratObject(counts = OX_DF_2, project = "OX_2E", min.features = 200)

# Add metadata information to the seurat object
OX_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-2", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_2E <- subset(OX_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_2E[["percent.mt"]] <- PercentageFeatureSet(OX_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_2E[["percent.pt"]] <- PercentageFeatureSet(OX_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_2E <- subset(OX_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_2E <- NormalizeData(OX_2E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_2E <- FindVariableFeatures(OX_2E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W2 <- GetAssayData(OX_2E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W2) <- paste("O2", colnames(OX_W2), sep = "_")


###
# OX - 3 E
###

# First replicate - OX 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
OX_3E <- CreateSeuratObject(counts = OX_DF_3, project = "OX_3E", min.features = 200)

# Add metadata information to the seurat object
OX_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-3", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_3E <- subset(OX_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_3E[["percent.mt"]] <- PercentageFeatureSet(OX_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_3E[["percent.pt"]] <- PercentageFeatureSet(OX_3E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_3E <- subset(OX_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_3E <- NormalizeData(OX_3E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_3E <- FindVariableFeatures(OX_3E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W3 <- GetAssayData(OX_3E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W3) <- paste("O3", colnames(OX_W3), sep = "_")

###
# OX - 7 E
###

# First replicate - OX 7E - total cells 9090; filter out genes that are not detected in at least 18 cells
OX_7E <- CreateSeuratObject(counts = OX_DF_7, project = "OX_7E", min.features = 200)

# Add metadata information to the seurat object
OX_7E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-7", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_7E <- subset(OX_7E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_7E[["percent.mt"]] <- PercentageFeatureSet(OX_7E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_7E[["percent.pt"]] <- PercentageFeatureSet(OX_7E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_7E <- subset(OX_7E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_7E <- NormalizeData(OX_7E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_7E <- FindVariableFeatures(OX_7E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
OX_W7 <- GetAssayData(OX_7E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W7) <- paste("O7", colnames(OX_W7), sep = "_")


###
# Intersection of highly variable genes - Cardamine hirsuta
###

# Lets find common variable features for the replicates of C. hirsuta
ox_hvgs = Reduce(intersect, list(OX_1E@assays$RNA@var.features, OX_2E@assays$RNA@var.features, OX_3E@assays$RNA@var.features, OX_7E@assays$RNA@var.features))

fileGenerator(ox_hvgs, fileName = "Shared_highly_variable_genes_between_replicates_CH.txt")

###
# COL - 1 E
###

# First replicate - COL 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
COL_1E <- CreateSeuratObject(counts = COL_data_1E, project = "COL_1E", min.features = 200)

# Add metadata information to the seurat object
COL_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-1", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_1E <- subset(COL_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_1E[["percent.mt"]] <- PercentageFeatureSet(COL_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_1E[["percent.pt"]] <- PercentageFeatureSet(COL_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_1E <- subset(COL_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_1E <- NormalizeData(COL_1E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
COL_1E <- FindVariableFeatures(COL_1E, selection.method = "vst", nfeatures = 3000)

# Extract the count table from the seurat object
COL_W1 <- GetAssayData(COL_1E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W1) <- paste("C1", colnames(COL_W1), sep = "_")


###
# COL - 2 E
###

# First replicate - COL 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
COL_2E <- CreateSeuratObject(counts = COL_data_2E, project = "COL_2E", min.features = 200)

# Add metadata information to the seurat object
COL_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-2", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_2E <- subset(COL_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_2E[["percent.mt"]] <- PercentageFeatureSet(COL_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_2E[["percent.pt"]] <- PercentageFeatureSet(COL_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_2E <- subset(COL_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_2E <- NormalizeData(COL_2E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
COL_2E <- FindVariableFeatures(COL_2E, selection.method = "vst", nfeatures = 3000)

# Extract the count table from the seurat object
COL_W2 <- GetAssayData(COL_2E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W2) <- paste("C2", colnames(COL_W2), sep = "_")


###
# COL - 3 E
###

# First replicate - COL 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
COL_3E <- CreateSeuratObject(counts = COL_data_3E, project = "COL_3E", min.features = 200)

# Add metadata information to the seurat object
COL_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-3", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_3E <- subset(COL_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_3E[["percent.mt"]] <- PercentageFeatureSet(COL_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_3E[["percent.pt"]] <- PercentageFeatureSet(COL_3E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_3E <- subset(COL_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_3E <- NormalizeData(COL_3E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
COL_3E <- FindVariableFeatures(COL_3E, selection.method = "vst", nfeatures = 3000)

# Extract the count table from the seurat object
COL_W3 <- GetAssayData(COL_3E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W3) <- paste("C3", colnames(COL_W3), sep = "_")


###
# COL - 5 E
###

# First replicate - COL 5E - total cells 8420; filter out genes that are not detected in at least 18 cells
COL_5E <- CreateSeuratObject(counts = COL_data_5E, project = "COL_5E", min.features = 200)

# Add metadata information to the seurat object
COL_5E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-5", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_5E <- subset(COL_5E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_5E[["percent.mt"]] <- PercentageFeatureSet(COL_5E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_5E[["percent.pt"]] <- PercentageFeatureSet(COL_5E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_5E <- subset(COL_5E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_5E <- NormalizeData(COL_5E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
COL_5E <- FindVariableFeatures(COL_5E, selection.method = "vst", nfeatures = 3000)

# Extract the count table from the seurat object
COL_W5 <- GetAssayData(COL_5E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W5) <- paste("C5", colnames(COL_W5), sep = "_")


###
# Intersection of highly variable genes - Arabidopsis thaliana
###

# Lets find common variable features for the replicates of A. thaliana
col_hvgs = Reduce(intersect, list(COL_1E@assays$RNA@var.features, COL_2E@assays$RNA@var.features, COL_3E@assays$RNA@var.features, COL_5E@assays$RNA@var.features))

fileGenerator(col_hvgs, fileName = "Shared_highly_variable_genes_between_replicates_AT.txt")

###
# Combine the two sets of highly variable genes - Arabidopsis thaliana and Cardamine hirsuta
###

# Combine the HVGs
HVGs_combined = union(ox_hvgs, col_hvgs) # Total = 2294

# fileGenerator(HVGs_combined, "Seurat_3000_HVG_intersect_reps_union_species_without_mincells_2294.txt")

###
# Creating a Liger object, pre-processing, and performing integration
###

# Lets create the liger object
WT_Species <- createLiger(list(WO1 = OX_W1, WO2 = OX_W2, WO3 = OX_W3, WO7 = OX_W7, WC1 = COL_W1, WC2 = COL_W2, WC3 = COL_W3, WC5 = COL_W5), remove.missing = F)

# Add metadata information to the liger object for the replicates in the same way as it was added in the seurat object

# Add replicate information
WT_Species@cell.data$Replicates <- WT_Species@cell.data$dataset
WT_Species@cell.data$Replicates <- factor(WT_Species@cell.data$Replicates, levels = c("WO1", "WO2", "WO3", "WO7", "WC1", "WC2", "WC3", "WC5"), labels = c("WT-OX-1", "WT-OX-2", "WT-OX-3", "WT-OX-7", "WT-COL-1", "WT-COL-2", "WT-COL-3", "WT-COL-5"))

# Add species information
WT_Species@cell.data$Species <- str_sub(WT_Species@cell.data$dataset, 1, nchar(WT_Species@cell.data$dataset) - 1)
WT_Species@cell.data$Species <- factor(WT_Species@cell.data$Species, levels = c("WO", "WC"), labels = c("Hirsuta", "Thaliana"))

# Add genotype information
WT_Species@cell.data$Genotype <- "WT"
WT_Species@cell.data$Genotype <- factor(WT_Species@cell.data$Genotype)

# Add tissue information
WT_Species@cell.data$Tissue <- "Leaf"
WT_Species@cell.data$Tissue <- factor(WT_Species@cell.data$Tissue)

WT_Species@norm.data <- list(WO1 = OX_1E@assays$RNA@data,
                             WO2 = OX_2E@assays$RNA@data, 
                             WO3 = OX_3E@assays$RNA@data, 
                             WO7 = OX_7E@assays$RNA@data, 
                             WC1 = COL_1E@assays$RNA@data, 
                             WC2 = COL_2E@assays$RNA@data, 
                             WC3 = COL_3E@assays$RNA@data, 
                             WC5 = COL_5E@assays$RNA@data)

# Selecting a set of highly variable genes
# WT_Species <- selectGenes(WT_Species, num.genes = 2000, do.plot = FALSE)

WT_Species@var.genes <- HVGs_combined

# Scale the feature count
WT_Species <- scaleNotCenter(WT_Species)

# Check which datasets are we integrating
table(WT_Species@cell.data$dataset)

#### Running suggest K
sg_K <- suggestK(WT_Species, k.test = seq(10, 60, 5), lambda = 5, num.cores = 10)

save(sg_K, file = "suggest_K_wt_species.Rdata")

writeLines(capture.output(sessionInfo()), "Session_info_suggest_K_wt_species.txt")

