# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2023/Beginning_of_a_compendium/a_scRNAseq_analysis_pipeline/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load the Orthologues table - Contains orthologues ID of Arabidopsis thaliana and Cardamine hirsuta.
# The orthologues table has been generated with the TAIR10 Arabidopsis thaliana gene IDs
###

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")


###
# Load the protoplasting induced genes - contain genes whose expression was altered due to the protoplasting procedure
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Ox_Co0_leaf_protoplast_v12_final_August_2022_ortho.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

# Load the orthologues data - saved in the previous step
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Input_Data/Orthologues_data/Arabidopsis_REP_3_Ortho_data.RData")

thaliana_genes = rownames(Ortho_COL_REP_3)

# Remove protoplasting-induced genes from the total set of Arabidopsis data
genes_to_keep = setdiff(thaliana_genes, PP_genes)

# Subsetting the data without the protoplasting induced genes
Ortho_COL_REP_3 <- Ortho_COL_REP_3[genes_to_keep, ]

###
# COL - 1st Experiment
###

# First replicate - COL 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
COL_3E <- CreateSeuratObject(counts = Ortho_COL_REP_3, project = "COL_3E", min.features = 200)

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


generate_data_summary_plot(COL_3E, covariate = "nCount_RNA", figure_name_suffix = "_AT_REP_3", binwidth = 1000)

generate_data_summary_plot(COL_3E, covariate = "nFeature_RNA", figure_name_suffix = "_AT_REP_3", binwidth = 100)

ggsave(plot = VlnPlot(COL_3E, features = "percent.mt") + 
         xlab("") + 
         ylab("percent.mt") + 
         ggtitle("Distribution of percent.mt") + 
         theme_classic() + 
         theme(axis.title.x = element_text(size = 18, face = "bold", colour = "black"), 
               axis.title.y = element_text(size = 18, face = "bold", colour = "black"), 
               axis.ticks.length = unit(.30, "cm"), 
               axis.text = element_text(size = 18, face = "bold", colour = "black"),
               title =  element_text(size = 18, face = "bold", colour = "black"),
               legend.position = "none"), 
       filename = "Mitochondrial_distribution_AT_REP_3.png", width = 12, height = 12, dpi = 300)

ggsave(plot = VlnPlot(COL_3E, features = "percent.pt") + 
         xlab("") + 
         ylab("percent.pt") + 
         ggtitle("Distribution of percent.pt") + 
         theme_classic() + 
         theme(axis.title.x = element_text(size = 18, face = "bold", colour = "black"), 
               axis.title.y = element_text(size = 18, face = "bold", colour = "black"), 
               axis.ticks.length = unit(.30, "cm"), 
               axis.text = element_text(size = 18, face = "bold", colour = "black"),
               title =  element_text(size = 18, face = "bold", colour = "black"),
               legend.position = "none"), 
       filename = "Chloroplast_distribution_AT_REP_3.png", width = 12, height = 12, dpi = 300)
