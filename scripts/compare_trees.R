library("optparse")

option_list <- list(
  make_option(c("-s","--scores"),
              type = "character", default = NULL,
              help = "scores"
  ),
  make_option(c("-t", "--trees"),
              type = "character", default = NULL,
              help = "tree file"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = "output filename", metavar = "character"
  )
)

desc <- "Produces informative plot on trees comparison"

opt_parser <- OptionParser(option_list = option_list, description = desc)
opt <- parse_args(opt_parser)

library(tidyverse)
library(ggtree)
library(patchwork)
library(ggpubr)
source("scripts/functions.R")

# opt <- NULL
# opt$scores <- "results/Hsap_opistho/reco/UP000005640_scores.tsv"
# opt$trees <- "results/Hsap_opistho/trees/UP000005640_unrooted_trees.txt"

theme_set(theme_bw())

df_trees <- read_delim(opt$trees, col_names = c("gene", "target", "alphabet", "tree"))

ts <- read.tree(text = df_trees$tree)
names(ts) <- paste(df_trees$gene, df_trees$target, df_trees$alphabet, sep = "_")

df_bs <- fortify(ts) %>%
  filter(!isTip) %>% 
  mutate(support=as.numeric(label)) %>% 
  separate(.id, c("gene", "target", "alphabet")) %>% 
  filter(alphabet!="fident")

plot_bs <- ggplot(df_bs, aes(support, target, fill=alphabet)) +
  geom_violin(position = position_dodge(width = 1)) +
  ggdist::stat_pointinterval(position = position_dodge(width = 1), color="black") +
  scale_fill_manual(values = palettes_alphabet) + 
  theme_classic() + 
  labs(y="") +
  theme(legend.position = "bottom")


rf_df <- NULL

# for this you should filter for trees with at least 10 (or 50%) common?

for (gene in unique(df_trees$gene)){
  print(gene)
  single_gene <- df_trees[df_trees$gene==gene, ]
  t <- read.tree(text = single_gene$tree)
  
  # do trees have at least 4 in common?
  labs <- sapply(t, function(x) x$tip.label)
  names(labs) <- letters[1:9]
  min_int <- min(crossprod(table(stack(labs))))

  if (min_int < 5){
    print("Less then 5 sequences")
    next
  }
  m <- TreeDist::TreeDistance(t, t)
  rownames(m) <- paste(single_gene$alphabet, single_gene$target, sep="_")
  colnames(m) <- paste(single_gene$alphabet, single_gene$target, sep="_")
  out_df <- as_tibble(m, rownames = "method") %>% 
    pivot_longer(!method, names_to = "method2") %>% 
    filter(method>=method2) %>% 
    mutate(gene=gene)
  rf_df <- rbind(rf_df, out_df)
}


rf_text <- group_by(rf_df, method, method2) %>% 
  summarise(value=median(value, na.rm=T)) %>% 
  separate(method2, c("pipeline","targets")) %>% 
  separate(method, c("pipeline2", "targets2"))

# filter(rf_df, method!=method2) %>% 
rf_plot <- rf_df %>% 
  separate(method, c("pipeline","targets")) %>% 
  separate(method2, c("pipeline2", "targets2")) %>% 
  ggplot(aes(value)) +
  geom_rug(linewidth=.1) +
  geom_tile(aes(x=0.5, y=0.5, fill=value), data = rf_text) +
  geom_text(aes(x=0.5, y=.5, label=round(value, 2)), color="white", data = rf_text) +
  geom_density(aes(y = after_stat(scaled)), data= . %>% filter(paste(targets, pipeline)!=paste(targets2, pipeline2))) + 
  scale_fill_distiller(palette="YlGnBu", direction = 1, limits = c(0,1)) + 
  scale_x_continuous(breaks = c(0,.5,1)) +
  coord_fixed(xlim = c(0,1), expand=0) +
  # facet_grid(method~method2, switch = "both") +
  ggh4x::facet_nested(pipeline+targets~pipeline2+targets2, nest_line = element_line(linetype = 1)) +
  labs(x="Clustering Information distance") +
  theme_classic() +
  theme(legend.position = "none", 
        panel.background = element_rect(fill='transparent'),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # strip.text = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())


# Reconciliation ranger+TCS scores

scores <- read_delim(opt$scores, show_col_types = FALSE)

# Ranger results

ranger_plot <- mutate(scores, n_events = (dups+losses)/n_tips) %>% 
  ggplot(aes(n_events, targets, fill=alphabet)) + 
  geom_boxplot() +
  scale_fill_manual(values = palettes_alphabet) + 
  labs(y="", x="(D+L)/# Tips") +
  guides(fill = guide_legend(nrow = 3)) + 
  theme(legend.position = "bottom")

TCS_df <- scores %>% 
  select(1:3,8) %>% 
  pivot_wider(names_from = alphabet, values_from = score)

aa_3di <- ggplot(TCS_df, aes(aa, `3Di`, color=targets)) + 
  geom_point() + 
  geom_abline(intercept = 1) +
  coord_fixed() +
  labs(title="TCS score") +
  guides(color = guide_legend(nrow = 3)) + 
  scale_color_manual(values = palettes_method) 
aa_fident <- ggplot(TCS_df, aes(aa, fident, color=targets)) + 
  geom_point() +
  geom_abline(intercept = 1) +
  coord_fixed() + 
  guides(color = guide_legend(nrow = 3)) + 
  scale_color_manual(values = palettes_method) 
fident_3di <- ggplot(TCS_df, aes(fident, `3Di`, color=targets)) + 
  geom_point() + 
  geom_abline(intercept = 1) +
  coord_fixed() +
  guides(color = guide_legend(nrow = 3)) + 
  scale_color_manual(values = palettes_method) 

plot_tcs <- (aa_3di / aa_fident / fident_3di) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

plot_varr2t <- ggplot(scores, aes(variance_r2t, color=alphabet)) +
  geom_density() + 
  facet_grid(targets~.) +
  scale_x_continuous(limits = c(0,.25)) + 
  scale_color_manual(values = palettes_alphabet) +
  theme(legend.position = 'none')


final_plot <- (plot_bs | rf_plot | plot_tcs | ranger_plot | plot_varr2t) + plot_layout(widths = c(1, 2.6, 1, 1, 1))
ggsave(opt$output, final_plot, width = 20, height = 10)

# ggsave("test/rf/rfplot.png", rf_plot, width = 7, height = 7, dpi=300, bg = "white")
# 
# df <- full_join(df_struct, df_seq) %>% 
#   full_join(df_structml) %>% 
#   mutate(overlap=NA, n_seqs=NA, n_structs=NA, jaccard=NA, over_b=NA, over_f=NA, 
#          RF_seq_struct=NA, RF_seq_structml=NA, RF_struct_structml=NA)
# 
# for (row in 1:nrow(df)){
#   if (!is.na(df[row,]$tree_struct)) {
#     tst <- read.tree(text=df[row,]$tree_struct)
#     tst_labs <- tst$tip.label
#   }
#   if (!is.na(df[row,]$tree_structml)) {
#     tstml <- read.tree(text=df[row,]$tree_structml)
#     tstml_labs <- tstml$tip.label
#   }
#   if (!is.na(df[row,]$tree_seq)) {
#     tse <- read.tree(text=df[row,]$tree_seq)
#     tse_labs <- tse$tip.label
#   }
#   
#   if (!is.na(df[row,]$tree_struct) & !is.na(df[row,]$tree_seq)) {
#     overlap <- length(intersect(tst_labs, tse_labs))
#     jaccard <- overlap/(length(tse_labs)+length(tst_labs)-overlap)
#     df[row, ]$overlap <- overlap
#     df[row, ]$n_seqs <- length(tse_labs)
#     df[row, ]$n_structs <- length(tst_labs)
#     df[row, ]$jaccard <- jaccard
#     df[row, ]$over_b <- overlap/length(tse_labs)
#     df[row, ]$over_f <- overlap/length(tst_labs)
#     if (overlap > 3){
#       df[row, ]$RF_seq_struct <- TreeDist::TreeDistance(tst,tse)
#     }
#   }
#   if (!is.na(df[row,]$tree_structml) & !is.na(df[row,]$tree_seq)) {
#     if (overlap > 3){
#       df[row, ]$RF_seq_structml <- TreeDist::TreeDistance(tse,tstml)
#     }
#   }
#   if (!is.na(df[row,]$tree_structml) & !is.na(df[row,]$tree_struct)) {
#     if (overlap > 3){
#       df[row, ]$RF_struct_structml <- TreeDist::TreeDistance(tst,tstml)
#     }
#   }
# }
# 
# rf_plot <- select(df, gene, jaccard, contains("RF")) %>% 
#   pivot_longer(cols = contains("RF")) %>%
#   mutate(name = gsub("_", " vs ", gsub("RF_", "", name))) %>% 
#   ggplot(aes(name, value)) + 
#   geom_violin(aes(fill=name)) + 
#   geom_boxplot(alpha=.7, color="grey20", width=.1) +
#   geom_line(aes(group=gene), linetype=2, alpha=.1) +
#   scale_y_continuous(limits = c(0,1)) +
#   theme(legend.position = "none") + 
#   labs(y="RF Similarity")

# facet_grid(~name)
# ranger <- read_delim(opt$ranger, show_col_types = FALSE)
# 
# my_comparisons <- as.list(as.data.frame((combn(unique(ranger$method), 2))))
# 
# ranger_plot <- ggboxplot(ranger, x = "method", y = "dups",
#           color = "method", palette = "jco")+
#   stat_compare_means(comparisons = my_comparisons) + 
#   theme(legend.position = "none")
# 
# final_plot <- rf_plot | ranger_plot
# ggsave(opt$outfile, final_plot, width = 7, height = 5)

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
