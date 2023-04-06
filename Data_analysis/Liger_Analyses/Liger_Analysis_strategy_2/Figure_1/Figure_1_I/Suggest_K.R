# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A")

###
# Load data tables
###

load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/Suggest_K_and_Lambda/Suggest_K/suggest_K_wt_species.Rdata")

df = sg_K$data
df = df[df$calc == "KL_div", ]
p <- ggplot(df, aes(x = k, y = median_kl)) + geom_line(color = grp_col[3], size = 1.2) + geom_point(color = grp_col[4]) + 
  xlab("K (GEPs)") + 
  ylab("Median KL divergence (across all cells)") + 
  scale_x_continuous(breaks = seq(10, 60, 10)) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title = element_text(size = 32, face = "bold", colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 28, face = "bold", colour = "black"),
    title = element_text(size = 32, face = "bold", colour = "black"))

ggsave(filename = "Figure_1I.png", plot = p, width = 12, height = 12, dpi = 300)

# Sources
# https://github.com/welch-lab/liger/issues/103
# https://github.com/welch-lab/liger/issues/116
# https://www.rdocumentation.org/packages/rliger/versions/0.5.0/topics/suggestK
# https://broadinstitute.github.io/2020_scWorkshop/batch-correction-lab.html

