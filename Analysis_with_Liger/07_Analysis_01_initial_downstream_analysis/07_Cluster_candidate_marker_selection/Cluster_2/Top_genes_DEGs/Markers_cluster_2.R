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

# Load the markers from GEPs
GEP_markers = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_2/Liger_analysis_with_CC/Cluster_DEGs_Markers_Candidates/Cluster_candidates/Cluster_2/Top_genes_GEP_identity/GEP_23/Cluster_2_GEP_23_markers.csv")

rownames(GEP_markers) = GEP_markers$gene_ID

geps_description = GEP_markers[, -c(4, 7:10)]

cluster_ids = levels(Idents(integrated.data))

target_ID = "2"

cluster_index = 3

# Markers from cluster DEGs
# Dir - containing the DEG files
DEG_dir = "/home/yasir/Documents/Thesis_PhD/Chapter_2/Liger_analysis_with_CC/Differentially_expressed_genes/Conserved_marker_grouped_by_Species/Conserved_markers_DEtest_wilcox"

DEG_files = str_sort(list.files(DEG_dir, pattern = "Cluster"), numeric = TRUE)

cluster_DEG = loadRData(str_c(DEG_dir, "/", DEG_files[cluster_index]))

# Cluster 0 marker
specific_markers_list = loadRData("/home/yasir/Documents/Thesis_PhD/Chapter_2/Liger_analysis_with_CC/Cluster_DEGs_Markers_Candidates/DEGs_and_Markers.RData")[[2]]

names(specific_markers_list)[cluster_index]

specfic_markers = specific_markers_list[[cluster_index]]

Cluster_marker = cluster_DEG[rownames(cluster_DEG) %in% specfic_markers, ]
Cluster_marker$gene_ID <- rownames(Cluster_marker)

# Create ortho genes list and gene description file
description = ATH[rownames(ATH) %in% rownames(Cluster_marker), -c(3, 4)]

# Description - All DEGs identified for the cluster
markers_description = merge(Cluster_marker, description, by.x = "gene_ID", by.y = "AT_ID", all.x = TRUE)

# Extract all the columns in the DEG file containing pct.2 information
temp_all_pct2 = markers_description[ , str_detect(colnames(markers_description), pattern = "pct.2")]

markers_description = markers_description[order(markers_description[, colnames(temp_all_pct2)[1]], markers_description[, colnames(temp_all_pct2)[2]]), ]
rownames(markers_description) = markers_description$gene_ID

# For TFs
markers_description_TF_subset = markers_description[!is.na(markers_description$TF), ]

avg_FC_names = colnames(markers_description)[str_detect(colnames(markers_description), pattern = "avg_log2FC")]

if (nrow(markers_description_TF_subset) == 0) {
  print("No TFs")} else if (nrow(markers_description_TF_subset) <= 10 & nrow(markers_description_TF_subset) != 0){
    markers_description_TF_subset = markers_description_TF_subset
  } else {
    if (nrow(markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]) < 10 & 
        nrow(markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]) > 0){
      subset1 = markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]
      subset2 = markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] < 1 & markers_description_TF_subset[[avg_FC_names[2]]] < 1, ]
      subset2 = subset2[order(subset2[, avg_FC_names[1]], subset2[, avg_FC_names[2]], decreasing = TRUE), ]
      
      if (nrow(subset2) <= (10 - nrow(subset1))) {
        markers_description_TF_subset = rbind.data.frame(subset1, subset2)
      } else {
        subset2 = head(subset2, (10 - nrow(subset1)))
        markers_description_TF_subset = rbind.data.frame(subset1, subset2)
      }
    }
    else if (nrow(markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]) == 0){
      markers_description_TF_subset = markers_description_TF_subset[order(markers_description_TF_subset[, avg_FC_names[1]], markers_description_TF_subset[, avg_FC_names[2]], decreasing = TRUE), ]
      markers_description_TF_subset = markers_description_TF_subset[c(1:10), ]
    }
    else {
      markers_description_TF_subset = markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]
      markers_description_TF_subset = markers_description_TF_subset[c(1:10), ]
    }
  }

rownames(markers_description_TF_subset) <- markers_description_TF_subset$gene_ID

write.csv(markers_description_TF_subset, str_c("Cluster_", target_ID, "_top_10_TFs_all_details.csv"), row.names = FALSE)

pct1_names = colnames(markers_description_TF_subset)[str_detect(colnames(markers_description_TF_subset), pattern = "pct.1")]

markers_description_TF_subset$pct.1 = apply(markers_description_TF_subset[, pct1_names], 1, mean)

pct2_names = colnames(markers_description_TF_subset)[str_detect(colnames(markers_description_TF_subset), pattern = "pct.2")]

markers_description_TF_subset$pct.2 = apply(markers_description_TF_subset[, pct2_names], 1, mean)

markers_description_TF_subset$source = "DEGs"

markers_description_TF_subset = markers_description_TF_subset[, c(1, 14, 17, 19, 20, 15, 16, 18, 21)]

write.csv(markers_description_TF_subset, str_c("Cluster_", target_ID, "_top_10_TFs.csv"), row.names = FALSE)

if (length(intersect(rownames(geps_description), rownames(markers_description_TF_subset))) != 0){
  
  # For TFs
  markers_description_TF_subset = markers_description[!is.na(markers_description$TF), ]
  
  shared_gene = intersect(rownames(geps_description), markers_description_TF_subset$gene_ID)
  
  geps_description[shared_gene, "source"] = str_c(geps_description[shared_gene, "source"], " and DEGs")
  
  markers_description_TF_subset = markers_description_TF_subset[!markers_description_TF_subset$gene_ID %in% shared_gene, ]
  
  if (nrow(markers_description_TF_subset) == 0) {
    print("No TFs")} else if (nrow(markers_description_TF_subset) <= 10 & nrow(markers_description_TF_subset) != 0){
      markers_description_TF_subset = markers_description_TF_subset
    } else {
      if (nrow(markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]) < 10 & 
          nrow(markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]) > 0){
        subset1 = markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]
        subset2 = markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] < 1 & markers_description_TF_subset[[avg_FC_names[2]]] < 1, ]
        subset2 = subset2[order(subset2[, avg_FC_names[1]], subset2[, avg_FC_names[2]], decreasing = TRUE), ]
        
        if (nrow(subset2) <= (10 - nrow(subset1))) {
          markers_description_TF_subset = rbind.data.frame(subset1, subset2)
        } else {
          subset2 = head(subset2, (10 - nrow(subset1)))
          markers_description_TF_subset = rbind.data.frame(subset1, subset2)
        }
      }
      else if (nrow(markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]) == 0){
        markers_description_TF_subset = markers_description_TF_subset[order(markers_description_TF_subset[, avg_FC_names[1]], markers_description_TF_subset[, avg_FC_names[2]], decreasing = TRUE), ]
        markers_description_TF_subset = markers_description_TF_subset[c(1:10), ]
      }
      else {
        markers_description_TF_subset = markers_description_TF_subset[markers_description_TF_subset[[avg_FC_names[1]]] >= 1 & markers_description_TF_subset[[avg_FC_names[2]]] >= 1, ]
        markers_description_TF_subset = markers_description_TF_subset[c(1:10), ]
      }
    }
  
  write.csv(markers_description_TF_subset, str_c("Cluster_", target_ID, "_top_10_TFs_all_details.csv"), row.names = FALSE)
  
  markers_description_TF_subset$pct.1 = apply(markers_description_TF_subset[, pct1_names], 1, mean)
  
  markers_description_TF_subset$pct.2 = apply(markers_description_TF_subset[, pct2_names], 1, mean)
  
  markers_description_TF_subset$source = "DEGs"
  
  markers_description_TF_subset = markers_description_TF_subset[, c(1, 14, 17, 19, 20, 15, 16, 18, 21)]
  
  write.csv(markers_description_TF_subset, str_c("Cluster_", target_ID, "_top_10_TFs.csv"), row.names = FALSE)
  
}

# For genes 
markers_description_genes_subset = markers_description[is.na(markers_description$TF), ]

if (nrow(markers_description_genes_subset) == 0) {
  print("No TFs")} else if (nrow(markers_description_genes_subset) <= 10 & nrow(markers_description_genes_subset) != 0){
    markers_description_genes_subset = markers_description_genes_subset
  } else {
    if (nrow(markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]) < 10 & 
        nrow(markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]) > 0){
      subset1 = markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]
      subset2 = markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] < 1 & markers_description_genes_subset[[avg_FC_names[2]]] < 1, ]
      subset2 = subset2[order(subset2[, avg_FC_names[1]], subset2[, avg_FC_names[2]], decreasing = TRUE), ]
      
      if (nrow(subset2) <= (10 - nrow(subset1))) {
        markers_description_genes_subset = rbind.data.frame(subset1, subset2)
      } else {
        subset2 = head(subset2, (10 - nrow(subset1)))
        markers_description_genes_subset = rbind.data.frame(subset1, subset2)
      }
    }
    else if (nrow(markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]) == 0){
      markers_description_genes_subset = markers_description_genes_subset[order(markers_description_genes_subset[, avg_FC_names[1]], markers_description_genes_subset[, avg_FC_names[2]], decreasing = TRUE), ]
      markers_description_genes_subset = markers_description_genes_subset[c(1:10), ]
    }
    else {
      markers_description_genes_subset = markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]
      markers_description_genes_subset = markers_description_genes_subset[c(1:10), ]
    }
  }

rownames(markers_description_genes_subset) <- markers_description_genes_subset$gene_ID

write.csv(markers_description_genes_subset, str_c("Cluster_", target_ID, "_top_10_genes_all_details.csv"), row.names = FALSE)

markers_description_genes_subset$pct.1 = apply(markers_description_genes_subset[, pct1_names], 1, mean)

markers_description_genes_subset$pct.2 = apply(markers_description_genes_subset[, pct2_names], 1, mean)

markers_description_genes_subset$source = "DEGs"

markers_description_genes_subset = markers_description_genes_subset[, c(1, 14, 17, 19, 20, 15, 16, 18, 21)]

write.csv(markers_description_genes_subset, str_c("Cluster_", target_ID, "_top_10_genes.csv"), row.names = FALSE)


if (length(intersect(rownames(geps_description), rownames(markers_description_genes_subset))) != 0){
  
  markers_description_genes_subset = markers_description[is.na(markers_description$TF), ]
  
  shared_gene = intersect(rownames(geps_description), markers_description_genes_subset$gene_ID)
  
  geps_description[shared_gene, "source"] = str_c(geps_description[shared_gene, "source"], " and DEGs")
  
  markers_description_genes_subset = markers_description_genes_subset[!markers_description_genes_subset$gene_ID %in% shared_gene, ]
  
  if (nrow(markers_description_genes_subset) == 0) {
    print("No TFs")} else if (nrow(markers_description_genes_subset) <= 10 & nrow(markers_description_genes_subset) != 0){
      markers_description_genes_subset = markers_description_genes_subset
    } else {
      if (nrow(markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]) < 10 & 
          nrow(markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]) > 0){
        subset1 = markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]
        subset2 = markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] < 1 & markers_description_genes_subset[[avg_FC_names[2]]] < 1, ]
        subset2 = subset2[order(subset2[, avg_FC_names[1]], subset2[, avg_FC_names[2]], decreasing = TRUE), ]
        
        if (nrow(subset2) <= (10 - nrow(subset1))) {
          markers_description_genes_subset = rbind.data.frame(subset1, subset2)
        } else {
          subset2 = head(subset2, (10 - nrow(subset1)))
          markers_description_genes_subset = rbind.data.frame(subset1, subset2)
        }
      }
      else if (nrow(markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]) == 0){
        markers_description_genes_subset = markers_description_genes_subset[order(markers_description_genes_subset[, avg_FC_names[1]], markers_description_genes_subset[, avg_FC_names[2]], decreasing = TRUE), ]
        markers_description_genes_subset = markers_description_genes_subset[c(1:10), ]
      }
      else {
        markers_description_genes_subset = markers_description_genes_subset[markers_description_genes_subset[[avg_FC_names[1]]] >= 1 & markers_description_genes_subset[[avg_FC_names[2]]] >= 1, ]
        markers_description_genes_subset = markers_description_genes_subset[c(1:10), ]
      }
    }
  
  rownames(markers_description_genes_subset) <- markers_description_genes_subset$gene_ID
  
  write.csv(markers_description_genes_subset, str_c("Cluster_", target_ID, "_top_10_genes_all_details.csv"), row.names = FALSE)
  
  markers_description_genes_subset$pct.1 = apply(markers_description_genes_subset[, pct1_names], 1, mean)
  
  markers_description_genes_subset$pct.2 = apply(markers_description_genes_subset[, pct2_names], 1, mean)
  
  markers_description_genes_subset$source = "DEGs"
  
  markers_description_genes_subset = markers_description_genes_subset[, c(1, 14, 17, 19, 20, 15, 16, 18, 21)]
  
  write.csv(markers_description_genes_subset, str_c("Cluster_", target_ID, "_top_10_genes.csv"), row.names = FALSE)
  
}

cluster_biomarker = rbind(geps_description, markers_description_TF_subset, markers_description_genes_subset)
cluster_biomarker = cluster_biomarker[order(cluster_biomarker$pct.2), ]
cluster_biomarker$rank = c(1:nrow(cluster_biomarker))

write.csv(cluster_biomarker, file = str_c("Cluster_", target_ID, "_bio_markers.csv"))

