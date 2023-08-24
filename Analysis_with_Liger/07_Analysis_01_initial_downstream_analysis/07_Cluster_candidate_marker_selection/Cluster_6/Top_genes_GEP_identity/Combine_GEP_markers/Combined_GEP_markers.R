# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_2/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# load and combine
file_list = list.files(pattern = "Cluster_")

GEP_markers = do.call(rbind.data.frame, lapply(file_list, read.csv))
rownames(GEP_markers) <- GEP_markers$gene_ID

# Order the data
GEP_markers = GEP_markers[order(GEP_markers$pct.2), ]
GEP_markers$rank = c(1:nrow(GEP_markers))
write.csv(GEP_markers, "GEP_markers_cluster_6.csv", row.names = FALSE)
