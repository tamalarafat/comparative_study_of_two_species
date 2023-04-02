source("/home/yasir/Documents/comparative_study_of_two_species/Functions/Load_libraries.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/load_rdata.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/order_factorized_matrix.R")

load("/home/yasir/Documents/comparative_study_of_two_species/Strategy_2/seurat_object_of_K_50.RData")

basis = as.data.frame(integrated.data@reductions[["inmf"]]@feature.loadings)
colnames(basis) = str_c("GEP_", c(1:50))

basis = order_factorized_matrix(basis)

mat_W = data.frame(Genes = rownames(basis), basis)
mat_W$Genes <- factor(mat_W$Genes, levels = mat_W$Genes)
mat_W <-
  mat_W %>% melt(
    id.vars = "Genes",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <- ggplot(mat_W, aes(x = GEPs, y = Genes, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Coefficient\nvalue",
    low = RColorBrewer::brewer.pal(9, "Blues")[1],
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(mat_W$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("GEPs") +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 32, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 28, face = "bold", colour = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 28, face = "bold", colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 28, face = "bold")
  )
ggsave(filename = "Figure_1B.png", plot = p, width = 32, height = 32, dpi = 300)
