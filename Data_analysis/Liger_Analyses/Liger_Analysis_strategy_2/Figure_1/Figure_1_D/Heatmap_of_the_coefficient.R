source("/home/yasir/Documents/comparative_study_of_two_species/Functions/Load_libraries.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/load_rdata.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/order_factorized_matrix.R")

load("/home/yasir/Documents/comparative_study_of_two_species/Strategy_2/seurat_object_of_K_50.RData")

coef = as.data.frame(integrated.data@reductions[["inmf"]]@cell.embeddings)
colnames(coef) = str_c("GEP_", c(1:50))

coef = order_factorized_matrix(coef)

mat_H = data.frame(Cells = rownames(coef), coef)
mat_H$Cells = factor(mat_H$Cells, levels = mat_H$Cells)

mat_H <-
  mat_H %>% melt(
    id.vars = "Cells",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(mat_H, aes(x = Cells, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion\nusage",
    low = RColorBrewer::brewer.pal(9, "Blues")[1],
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(mat_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("Cells") +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 32, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 28, face = "bold", colour = "black"),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 28, face = "bold", colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 28, face = "bold")
  )

ggsave(filename = "Figure_1D.png", plot = p, width = 44, height = 28, dpi = 300)

