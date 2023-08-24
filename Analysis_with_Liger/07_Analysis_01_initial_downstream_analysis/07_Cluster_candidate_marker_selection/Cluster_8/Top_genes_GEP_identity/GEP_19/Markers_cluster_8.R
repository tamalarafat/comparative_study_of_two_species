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

# GEP genes
GEPgenes = read.delim("GEP_19_genes.txt", header = FALSE, col.names = "gene_ID")

GEP_ID = 19

cluster_ids = levels(Idents(integrated.data))

target_ID = "8"

cluster_GEP_markers = FindMarkers(integrated.data, ident.1 = target_ID, features = GEPgenes$gene_ID, test.use = "wilcox", only.pos = TRUE)
cluster_GEP_markers$gene_ID = rownames(cluster_GEP_markers)

# Get the basis matrix
basis = as.data.frame(integrated.data@reductions[["inmf"]]@feature.loadings)
colnames(basis) = str_c("GEP_", parse_number(colnames(basis)))

# Order the matrix with all geps
basis = order_factorized_matrix(basis)

Basis_subset = basis[GEPgenes$gene_ID, GEP_ID, drop = FALSE]
colnames(Basis_subset)[1] <- "coefficient"
Basis_subset = Basis_subset[order(Basis_subset$coefficient, decreasing = TRUE), , drop = FALSE]
Basis_subset$rank = c(1:nrow(Basis_subset))
Basis_subset$gene_ID= rownames(Basis_subset)

# lets create a description file for the basis genes
description_basis_genes = ATH[rownames(ATH) %in% rownames(Basis_subset), ]

# Description - All DEGs identified for the cl/home/yasir/Documents/Chapter1/Liger_analysis_strategy_1/Results/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Cluster_biomarkers/Cluster_0/Top_genes_GEP_identity/GEP_27/Cluster_0_GEP_27_markers.csvuster
geps_description = merge(Basis_subset, description_basis_genes, by.x = "gene_ID", by.y = "AT_ID", all.x = TRUE)[-c(5, 6)]

geps_description = geps_description[order(geps_description$rank, decreasing = FALSE), ]

geps_description = merge(geps_description, cluster_GEP_markers, by = "gene_ID")

geps_description = geps_description[, c(1, 4, 7, 10, 11, 12, 2, 3, 9, 13, 5, 6, 8)]

geps_description = geps_description[order(geps_description$pct.2, decreasing = FALSE), ]

# if the selection criteria does not work, we picked the top 10 genes
geps_description = geps_description[geps_description$pct.2 < 0.1, ]

if (nrow(geps_description) == 0) {
  print("No genes")} else if (nrow(geps_description) <= 10 & nrow(geps_description) != 0){
    geps_description = geps_description
  } else {
    if (nrow(geps_description[geps_description$avg_log2FC >= 1,]) < 10 & nrow(geps_description[geps_description$avg_log2FC >= 1,]) > 0){
      subset1 = geps_description[geps_description$avg_log2FC >= 1,]
      subset2 = geps_description[geps_description$avg_log2FC < 1,]
      subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
      
      if (nrow(subset2) <= (10 - nrow(subset1))) {
        geps_description = rbind.data.frame(subset1, subset2)
      } else {
        subset2 = head(subset2, (10 - nrow(subset1)))
        geps_description = rbind.data.frame(subset1, subset2)
      }
    }
    else if (nrow(geps_description[geps_description$avg_log2FC >= 1,]) == 0){
      geps_description = geps_description[order(geps_description$avg_log2FC, decreasing = TRUE), ]
      geps_description = geps_description[c(1:10), ]
    }
    else {
      geps_description = geps_description[geps_description$avg_log2FC >= 1, ]
      geps_description = geps_description[c(1:10), ]
    }
  }

rownames(geps_description) <- geps_description$gene_ID

geps_description$source = str_c("GEP_", GEP_ID)

write.csv(geps_description, file = str_c("Cluster_", target_ID, "_GEP_", GEP_ID,"_markers.csv"), row.names = FALSE)

if (nrow(geps_description) < 10){
  
  subset1 = geps_description[ , !colnames(geps_description) %in% "source"]
  
  geps_description = merge(Basis_subset, description_basis_genes, by.x = "gene_ID", by.y = "AT_ID", all.x = TRUE)[-c(5, 6)]
  
  geps_description = geps_description[order(geps_description$rank, decreasing = FALSE), ]
  
  geps_description = merge(geps_description, cluster_GEP_markers, by = "gene_ID")
  
  geps_description = geps_description[, c(1, 4, 7, 10, 11, 12, 2, 3, 9, 13, 5, 6, 8)]
  
  geps_description = geps_description[order(geps_description$pct.2, decreasing = FALSE), ]
  
  subset2 = geps_description[!geps_description$gene_ID %in% rownames(subset1), ]
  
  if (nrow(subset2) <= (10 - nrow(subset1))) {
    geps_description = rbind.data.frame(subset1, subset2)
  } else {
    subset2 = head(subset2, (10 - nrow(subset1)))
    geps_description = rbind.data.frame(subset1, subset2)
  }
  
}


rownames(geps_description) <- geps_description$gene_ID

geps_description$source = str_c("GEP_", GEP_ID)

write.csv(geps_description, file = str_c("Cluster_", target_ID, "_GEP_", GEP_ID,"_markers.csv"), row.names = FALSE)


