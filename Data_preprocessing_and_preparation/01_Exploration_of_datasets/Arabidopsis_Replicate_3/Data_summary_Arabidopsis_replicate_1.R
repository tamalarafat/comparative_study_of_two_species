# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2023/Beginning_of_a_compendium/a_scRNAseq_analysis_pipeline/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Loading data -  loading of sparse data matrices provided by 10X genomics using Seurat's "Read10X" function
COL_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_3rd_ALL/filtered_feature_bc_matrix/")

###
# COL - 3 Experiment
###

# Arabidopsis thaliana (Col-0 or Columbia-0 is an Arabidopsis ecotype) third replicate - COL 3E
COL_3E <- CreateSeuratObject(counts = COL_data_3E, project = "COL_3E", min.features = 200)

# Add metadata information to the seurat object
COL_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-3", "WT", "Leaf")

generate_data_summary_plot(COL_3E, covariate = "nCount_RNA", figure_name_suffix = "_AT_REP_3", binwidth = 1000)

generate_data_summary_plot(COL_3E, covariate = "nFeature_RNA", figure_name_suffix = "_AT_REP_3", binwidth = 100)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_3E[["percent.mt"]] <- PercentageFeatureSet(COL_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_3E[["percent.pt"]] <- PercentageFeatureSet(COL_3E, pattern = "^ATC")

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



