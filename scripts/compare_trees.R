library("optparse")

option_list <- list(
  make_option(c("--seq"),
              type = "character", default = NULL,
              help = "sequence trees", dest = "seq"
  ),
  make_option(c("--struct"),
              type = "character", default = NULL,
              help = "structure trees", dest = "struct"
  ),
  make_option(c("--structml"),
              type = "character", default = NULL,
              help = "structure ml trees", dest = "structml"
  ),
  make_option(c("-r", "--ranger"),
              type = "character", default = NULL,
              help = "ranger results", dest = "ranger"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = "output filename", metavar = "character", dest = "outfile"
  )
)

desc <- "Produces informative plot on trees comparison"

opt_parser <- OptionParser(option_list = option_list, description = desc)
opt <- parse_args(opt_parser)

library(tidyverse)
library(ggtree)
library(patchwork)
library(ggpubr)

# opt <- NULL
# opt$ranger <- "test/tcs/prova"
# struct_trees <- "results/Ncra_opistho/trees/UP000001805_fs_foldtrees.txt"
# structml_trees <- "results/Ncra_opistho/trees/UP000001805_fs_mltrees.txt"
# seq_trees <- "results/Ncra_opistho/trees/UP000001805_blast_trees.txt"
# ranger_result <- "results/Ncra_opistho/reco/UP000001805_ranger.tsv"

theme_set(theme_bw())

df_struct <- read_delim(opt$struct, col_names = c("gene", "tree_struct"), show_col_types = FALSE)
# mutate(.id=factor(row_number()))
df_structml <- read_delim(opt$structml, col_names = c("gene", "tree_structml"), show_col_types = FALSE) 
df_seq <- read_delim(opt$seq, col_names = c("gene", "tree_seq"), show_col_types = FALSE)
# mutate(.id=factor(row_number()))

df <- full_join(df_struct, df_seq) %>% 
  full_join(df_structml) %>% 
  mutate(overlap=NA, n_seqs=NA, n_structs=NA, jaccard=NA, over_b=NA, over_f=NA, 
         RF_seq_struct=NA, RF_seq_structml=NA, RF_struct_structml=NA)

for (row in 1:nrow(df)){
  if (!is.na(df[row,]$tree_struct)) {
    tst <- read.tree(text=df[row,]$tree_struct)
    tst_labs <- tst$tip.label
  }
  if (!is.na(df[row,]$tree_structml)) {
    tstml <- read.tree(text=df[row,]$tree_structml)
    tstml_labs <- tstml$tip.label
  }
  if (!is.na(df[row,]$tree_seq)) {
    tse <- read.tree(text=df[row,]$tree_seq)
    tse_labs <- tse$tip.label
  }
  
  if (!is.na(df[row,]$tree_struct) & !is.na(df[row,]$tree_seq)) {
    overlap <- length(intersect(tst_labs, tse_labs))
    jaccard <- overlap/(length(tse_labs)+length(tst_labs)-overlap)
    df[row, ]$overlap <- overlap
    df[row, ]$n_seqs <- length(tse_labs)
    df[row, ]$n_structs <- length(tst_labs)
    df[row, ]$jaccard <- jaccard
    df[row, ]$over_b <- overlap/length(tse_labs)
    df[row, ]$over_f <- overlap/length(tst_labs)
    if (overlap > 3){
      df[row, ]$RF_seq_struct <- TreeDist::TreeDistance(tst,tse)
    }
  }
  if (!is.na(df[row,]$tree_structml) & !is.na(df[row,]$tree_seq)) {
    if (overlap > 3){
      df[row, ]$RF_seq_structml <- TreeDist::TreeDistance(tse,tstml)
    }
  }
  if (!is.na(df[row,]$tree_structml) & !is.na(df[row,]$tree_struct)) {
    if (overlap > 3){
      df[row, ]$RF_struct_structml <- TreeDist::TreeDistance(tst,tstml)
    }
  }
}

rf_plot <- select(df, gene, jaccard, contains("RF")) %>% 
  pivot_longer(cols = contains("RF")) %>% 
  ggplot(aes(name, value)) + 
  geom_violin(aes(fill=name)) + 
  geom_boxplot(alpha=.7, color="grey20", width=.1) +
  geom_line(aes(group=gene), linetype=2, alpha=.1) +
  theme(legend.position = "none") + 
  labs(y="RF Similarity")
# facet_grid(~name)
ranger <- read_delim(opt$ranger, show_col_types = FALSE)

my_comparisons <- as.list(as.data.frame((combn(unique(ranger$method), 2))))


ranger_plot <- ggboxplot(ranger, x = "method", y = "dups",
          color = "method", palette = "jco")+
  stat_compare_means(comparisons = my_comparisons) + 
  theme(legend.position = "none")

final_plot <- rf_plot | ranger_plot
ggsave(opt$outfile, final_plot, width = 7, height = 5)

# trees_struct <- read.tree(text=df_struct$tree)
# names(trees_struct) <- df_struct$gene
# 
# trees_seq <- read.tree(text=df_seq$tree)
# names(trees_seq) <- df_seq$gene
# 
# trees_struct_df <- fortify(trees_struct) %>% 
#   left_join(select(df_struct, -tree), by=c(".id"="gene"))
# trees_seq_df <- fortify(trees_seq) %>% 
#   left_join(select(df_seq, -tree), by=c(".id"="gene"))
