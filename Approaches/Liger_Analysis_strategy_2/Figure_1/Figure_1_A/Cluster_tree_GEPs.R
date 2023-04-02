source("/home/yasir/Documents/comparative_study_of_two_species/Functions/Load_libraries.R")
source("/home/yasir/Documents/comparative_study_of_two_species/Functions/load_rdata.R")

# # Load all the functions stored in scripts from the folder housing the scripts
# scripts_list <- list.files("/home/yasir/Documents/comparative_study_of_two_species/Functions", pattern = "*.R$", full.names = TRUE) 
# sapply(scripts_list, source, .GlobalEnv)


# load the gene cluster file that is stored in the basis matrix dataframe
load("/home/yasir/Documents/comparative_study_of_two_species/Strategy_2/Basis_objects/Basis.RData")

Basis_subset = Basis[, str_c("K_", c(10, 20, 30, 40, 50))]

p <- clustree(Basis_subset, prefix = "K_", node_text_size = 7, node_size = 10, prop_filter = 0.2) + 
  labs(colour = "Factorization\nK", edge_colour = "Count", edge_alpha = "Proportion") + 
  theme(legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 22, colour = "black", face = "bold"))

ggsave(filename = "Figure_1A.png", plot = p, width = 24, height = 20, dpi = 300)

