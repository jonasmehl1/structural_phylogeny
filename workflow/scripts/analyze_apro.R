suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtree))
# library(ggtern)
library(patchwork)
theme_set(theme_bw())
source("workflow/scripts/functions.R")

labeled_sptree <- read.tree(snakemake@input[["sptree"]])
nodes_clades <- fortify(labeled_sptree) %>% 
  filter(!isTip) %>% 
  mutate(label = fct_reorder(label, y)) %>% 
  select(node, label)

clades <- read_delim(snakemake@input[["groups"]], show_col_types = FALSE)

apro_files <- snakemake@input[["trees"]]
apro_trees <- read.tree(text = sapply(apro_files, readLines))
names(apro_trees) <- sapply(str_split(basename(apro_files), "_"), function(x) paste0(x[2], "_", x[4]))

trees_df <- fortify(apro_trees) 
nodes_df <- trees_df %>% 
  separate(.id, c("targets", "model")) %>% 
  filter(!isTip, label!="") %>% 
  mutate(label = gsub("\\[|\\]|\\'", "", str_replace_all(label, "[a-z]{1,2}[0-9]=", ""))) %>% 
  separate(label, c("pp1", "pp2", "pp3", "f1", "f2", "f3", "q1", "q2", "q3"), ";", convert = TRUE) %>% 
  left_join(nodes_clades) %>% 
  mutate(freq=f1+f2+f3)

pp1_plot <- nodes_df %>% 
    ggplot(aes(q1, label, fill=model, color=model, group=model)) +
    geom_linerange(aes(x=x, xmin=xmin, xmax=xmax, color=model),
                   data = . %>% 
                     group_by(model, label) %>% 
                     summarise(x=mean(q1), xmin=min(q1), xmax=max(q1)),
                   position = position_dodge(width=.7), color="grey70"
                   ) +
    geom_point(aes(shape=targets),
               position = position_dodge(width=.7)) +
  scale_fill_manual(values=palettes_model) +
  scale_color_manual(values = palettes_model_d) +
  scale_shape_manual(values = c(21,22,23)) +
  guides(color = guide_legend(override.aes = list(pch = 21, size=3)),
         shape = guide_legend(override.aes = list(size=3))) +
  theme(legend.position = "bottom")

sptree_plot <- ggtree(labeled_sptree) + 
  geom_tiplab() + 
  geom_nodelab(color="firebrick", 
               nudge_x = -.1, hjust = 1, 
               nudge_y = .2, size=2)+
  scale_x_continuous(expand = c(.2,.2))
  
final_plot <- (sptree_plot | pp1_plot) + plot_layout(widths = c(1,1.5))

ggsave(snakemake@output[[1]], final_plot, width = 12, height = 7)


#apro_files_sptrees <- list.files(opt$apro, pattern="*apro_sptree.nwk", full.names=TRUE)
#apro_sptrees <- read.tree(text = sapply(apro_files_sptrees, readLines))
#names(apro_sptrees) <- paste0(str_split_i(basename(apro_files_sptrees), "_", 2),"_",str_split_i(basename(apro_files), "_", 3))

#fortify(apro_sptrees) %>% 
#  separate(.id, c("method", "alphabet"), remove = FALSE) %>% 
#  ggtree() + 
#  geom_tiplab() +
#  facet_grid(method~alphabet)
