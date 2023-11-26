library(tidyverse)

theme_set(theme_bw())

seed <- "UP000001805"

lddt_file <- "results/Ncra_opistho/db/db_lddt.tsv"
lddt_proteomes <- read_delim(lddt_file)

ggplot(lddt_proteomes, aes(fct_reorder(proteome, mean_lddt), mean_lddt)) + 
  geom_violin(fill="#35B779") + 
  geom_boxplot(width=0.1, alpha=0, color="black", outlier.alpha = 1) +
  ylim(c(0,100)) +
  theme(axis.text.x = element_text(angle=70, hjust=1))

