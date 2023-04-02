source("/home/yasir/Documents/comparative_study_of_two_species/Functions/Load_libraries.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/load_rdata.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/order_factorized_matrix.R")

# # Load all the functions stored in scripts from the folder housing the scripts
# scripts_list <- list.files("/home/yasir/Documents/comparative_study_of_two_species/Functions", pattern = "*.R$", full.names = TRUE) 
# sapply(scripts_list, source, .GlobalEnv)

# Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A")

load("/home/yasir/Documents/comparative_study_of_two_species/Strategy_2/seurat_object_of_K_50.RData")

basis = as.data.frame(integrated.data@reductions[["inmf"]]@feature.loadings)
colnames(basis) = str_c("GEP_", c(1:50))

basis = order_factorized_matrix(basis)

Basis_subset = basis[, 1, drop = FALSE]
Basis_subset$GEP_1 = sort(Basis_subset$GEP_1)
Basis_subset$index = c(1:nrow(Basis_subset))

Basis_label = basis[, 1, drop = FALSE]
Basis_label$GEP_1 = sort(Basis_label$GEP_1)
Basis_label$index = c(1:nrow(Basis_label))
Basis_label = tail(Basis_label, 5)
Basis_label$gene_ID = rownames(Basis_label)

p <- ggplot(Basis_subset, aes(x = index, y = GEP_1)) + geom_point(color = grp_col[3]) +
  ggrepel::geom_text_repel(data = Basis_label, aes(label = gene_ID), nudge_x = -300, size = 8) +
  xlab("Genes") + ylab("GEP 1") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title = element_text(size = 32, face = "bold", color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 28, face = "bold", colour = "black")
  )
ggsave(filename = "Figure_1F_1.png", plot = p, width = 10, height = 16, dpi = 300)


Basis_subset = basis[, 5, drop = FALSE]
Basis_subset$GEP_5 = sort(Basis_subset$GEP_5)
Basis_subset$index = c(1:nrow(Basis_subset))

Basis_label = basis[, 5, drop = FALSE]
Basis_label$GEP_5 = sort(Basis_label$GEP_5)
Basis_label$index = c(1:nrow(Basis_label))
Basis_label = tail(Basis_label, 5)
Basis_label$gene_ID = rownames(Basis_label)

p <- ggplot(Basis_subset, aes(x = index, y = GEP_5)) + geom_point(color = grp_col[3]) +
  ggrepel::geom_text_repel(data = Basis_label, aes(label = gene_ID), nudge_x = -300, size = 8) +
  xlab("Genes") + ylab("GEP 5") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title = element_text(size = 32, face = "bold", color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 28, face = "bold", colour = "black")
  )
ggsave(filename = "Figure_1F_2.png", plot = p, width = 10, height = 16, dpi = 300)

