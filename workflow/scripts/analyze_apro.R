# library("optparse")
# 
# option_list <- list(
#   make_option(c("-s","--sptree"),
#               type = "character", default = NULL,
#               help = "species tree with internal node names"
#   ),
#   make_option(c("-m","--meta"),
#               type = "character", default = NULL,
#               help = "taxon information file"
#   ),
#   make_option(c("-a", "--apro"),
#               type = "character", default = NULL,
#               help = "astral pro directory"
#   ),
#   make_option(c("-o", "--output"),
#               type = "character", default = NULL,
#               help = "output filename", metavar = "character"
#   )
# )
# 
# desc <- "Produces informative plot on trees comparison"
# 
# opt_parser <- OptionParser(option_list = option_list, description = desc)
# opt <- parse_args(opt_parser)

# opt <- NULL
# opt$sptree <- "data/sptrees/homo_internal.spTree.nw"
# opt$meta <- "data/meta/Hsap_taxon.csv"
# opt$apro <- "results/Hsap_opistho/reco/"
# opt$output <- ""

library(tidyverse)
library(ggtree)
# library(ggtern)
library(patchwork)
theme_set(theme_bw())
source("workflow/scripts/functions.R")

labeled_sptree <- read.tree(snakemake@input[["sptree"]])
nodes_clades <- fortify(labeled_sptree) %>% 
  filter(!isTip) %>% 
  mutate(label = fct_reorder(label, y)) %>% 
  select(node, label)

clades <- read_delim(snakemake@input[["groups"]])

apro_files <- snakemake@input[["trees"]]
apro_trees <- read.tree(text = sapply(apro_files, readLines))
names(apro_trees) <- paste0(str_split_i(basename(apro_files), "_", 2),"_",str_split_i(basename(apro_files), "_", 3))


trees_df <- fortify(apro_trees) 
nodes_df <- trees_df %>% 
  separate(.id, c("method", "alphabet")) %>% 
  filter(!isTip, label!="") %>% 
  mutate(label = gsub("\\[|\\]|\\'", "", str_replace_all(label, "[a-z]{1,2}[0-9]=", ""))) %>% 
  separate(label, c("pp1", "pp2", "pp3", "f1", "f2", "f3", "q1", "q2", "q3"), ";", convert = TRUE) %>% 
  left_join(nodes_clades) %>% 
  mutate(freq=f1+f2+f3)

pp1_plot <- ggplot(nodes_df, aes(q1, label, color=method, shape=alphabet, size=freq)) + 
  geom_point() + 
  # facet_grid(~alphabet) + 
  scale_color_manual(values=palettes_method) + 
  scale_shape_manual(values = c(15,16,18)) +
  xlim(c(0,1)) +
  theme(legend.position = "bottom")


sptree_plot <- ggtree(labeled_sptree) + 
  geom_tiplab() + 
  geom_nodelab(nudge_x = -.1, hjust = 1, nudge_y = .2)+
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
