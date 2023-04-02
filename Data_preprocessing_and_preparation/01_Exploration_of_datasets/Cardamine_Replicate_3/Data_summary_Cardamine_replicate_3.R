# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2023/Beginning_of_a_compendium/a_scRNAseq_analysis_pipeline/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Loading data -  loading of sparse data matrices provided by 10X genomics using Seurat's "Read10X" function
OX_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_3rd_ALL_3000_Newest/filtered_feature_bc_matrix/")

###
# OX - 3rd Experiment
###

# First replicate - OX 1E - total cells 6640; filter out genes that are not detected in at least 13 cells
OX_3E <- CreateSeuratObject(counts = OX_data_3E, project = "OX_3E", min.features = 200)

# Add metadata information to the seurat object
OX_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-3", "WT", "Leaf")

generate_data_summary_plot(OX_3E, covariate = "nCount_RNA", figure_name_suffix = "_CH_REP_3", binwidth = 1000)

generate_data_summary_plot(OX_3E, covariate = "nFeature_RNA", figure_name_suffix = "_CH_REP_3", binwidth = 100)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_3E[["percent.mt"]] <- PercentageFeatureSet(OX_3E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_3E[["percent.pt"]] <- PercentageFeatureSet(OX_3E, pattern = "Pt")

ggsave(plot = VlnPlot(OX_3E, features = "percent.mt") + 
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
       filename = "Mitochondrial_distribution_CH_REP_3.png", width = 12, height = 12, dpi = 300)

ggsave(plot = VlnPlot(OX_3E, features = "percent.pt") + 
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
       filename = "Chloroplast_distribution_CH_REP_3.png", width = 12, height = 12, dpi = 300)



