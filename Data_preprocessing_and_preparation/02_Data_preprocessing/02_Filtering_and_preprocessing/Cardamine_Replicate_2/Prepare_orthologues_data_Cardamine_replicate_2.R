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


# The orthologues table has been generated using the TAIR10 genome, whereas the Arabidopsis scRNAseq datasets have been generated using the Araport 11 genome
# The Araport11 genome has 51 genes more than the TAIR10 genome.
# So, we are going to remove those genes from the Cardamine dataset

# Load the thaliana genes
thaliana_genes = read_lines("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Input_Data/Gene_sets/Arabidopsis_orthologues.txt")


# Load the orthologues data - saved in the previous step
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Input_Data/Orthologues_data/Cardamine_REP_1_Ortho_data.RData")

# Get the gene ids in the cardamine dataset
hirsuta_genes = rownames(Ortho_OX_REP_1)



###
# Filter out genes from the data
###

# not all the thaliana ids are present in the ortho data - 51 thaliana genes are missing in the thaliana data
hirsuta_genes_to_keep = intersect(hirsuta_genes, thaliana_genes)

# Remove protoplasting-induced genes from the total set of Arabidopsis data
genes_to_keep = setdiff(thaliana_genes, PP_genes)

# Subsetting the data without the protoplasting induced genes
Ortho_OX_REP_1 <- Ortho_OX_REP_1[genes_to_keep, ]


###
# OX - 1st Experiment
###

# Cardamine hirsuta (strain OXFORD) - First replicate
OX_1E <- CreateSeuratObject(counts = Ortho_OX_REP_1, project = "OX_1E", min.features = 200)

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

generate_data_summary_plot(OX_1E, covariate = "nCount_RNA", figure_name_suffix = "_CH_REP_1", binwidth = 1000)

generate_data_summary_plot(OX_1E, covariate = "nFeature_RNA", figure_name_suffix = "_CH_REP_1", binwidth = 100)

ggsave(plot = VlnPlot(OX_1E, features = "percent.mt") + 
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
       filename = "Mitochondrial_distribution_CH_REP_1.png", width = 12, height = 12, dpi = 300)

ggsave(plot = VlnPlot(OX_1E, features = "percent.pt") + 
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
       filename = "Chloroplast_distribution_CH_REP_1.png", width = 12, height = 12, dpi = 300)
