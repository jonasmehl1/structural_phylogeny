library("optparse")

option_list <- list(
  make_option(c("-b", "--blast"),
              type = "character", default = NULL,
              help = "blast file", dest = "blast"
  ),
  make_option(c("-f", "--fs"),
              type = "character", default = NULL,
              help = "foldseek file", dest = "fs"
  ),
  make_option(c("--fb"),
              type = "character", default = NULL,
              help = "foldseek brh", dest = "fs_brh"
  ),
  make_option(c("--bb"),
              type = "character", default = NULL,
              help = "blast brh", dest = "blast_brh"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = "output filename", metavar = "character", dest = "outfile"
  )
)

desc <- "Produces informative plot on foldseek and blast results comparison"

opt_parser <- OptionParser(option_list = option_list, description = desc)
opt <- parse_args(opt_parser)

library(tidyverse)
library(patchwork)
source("scripts/functions.R")
theme_set(theme_bw())

# opt <- NULL
# opt$fs <- "results/Ncra_opistho/homology/UP000001805_fs.tsv"
# opt$blast <- "results/Ncra_opistho/homology/UP000001805_blast.tsv"
# opt$outfile <- "results/Ncra_opistho/plots/plot_intersection.png"
# opt$fs_brh <- "results/Ncra_opistho/homology/UP000001805_fs_brh.tsv"
# opt$blast_brh <- "results/Ncra_opistho/homology/UP000001805_blast_brh.tsv"

fs <- read_delim(opt$fs, col_names = fs_columns, show_col_types = FALSE)
blast <- read_delim(opt$blast, col_names = blast_columns, show_col_types = FALSE)

seeds <- union(unique(fs$query), unique(blast$query))

evals <- c(Inf, 1e-2, 1e-3, 1e-5)
min_covs <- c(0, 30, 60)

fs_brhs <- read_delim(opt$fs_brh, col_names = c("query", "target"), show_col_types = FALSE) %>% 
  group_by(query) %>% 
  summarise(fs_brhs=list(target))

blast_brhs <- read_delim(opt$blast_brh, col_names = c("query", "target"), show_col_types = FALSE) %>% 
  group_by(query) %>% 
  summarise(blast_brhs=list(target))

brhs <- full_join(fs_brhs, blast_brhs, by="query")

df_overlap <- tibble()

for (cov in min_covs){
  for (eval1 in evals){
    for (eval2 in evals){
      
      print(paste(cov, "-", eval1, "-", eval2))
      
      a <- filter(blast, query %in% seeds, evalue<eval1, qcov >= cov) %>% 
        group_by(query) %>% 
        summarise(blast_hits=list(target), n_seqs=n())
      
      b <- filter(fs, query %in% seeds, evalue<eval2, qcov >= cov/100) %>% 
        group_by(query) %>% 
        summarise(fs_hits=list(target), n_structs=n())
      
      to_add <- full_join(a, b, by="query") %>% 
        full_join(brhs, by="query") %>% 
        rowwise() %>% 
        mutate(eval_blast = factor(eval1), 
               eval_fs = factor(eval2), 
               cov = cov,
               fs_unique = list(setdiff(fs_hits, blast_hits)), 
               blast_unique = list(setdiff(blast_hits, fs_hits)), 
               common_hits = list(intersect(blast_hits, fs_hits)),
               fs_brhs = list(intersect(fs_hits, fs_brhs)),
               n_fs_brhs = length(fs_brhs), 
               blast_brhs = list(intersect(blast_hits, blast_brhs)),
               n_blast_brhs = length(blast_brhs), 
               overlap = length(common_hits),
               jaccard = overlap/(n_seqs+n_structs-overlap),
               overlap_blast = overlap/n_seqs,
               overlap_fs = overlap/n_structs,
               blast_fs_brh = length(intersect(blast_hits, fs_brhs))/n_fs_brhs,
               fs_blast_brh = length(intersect(fs_hits, blast_brhs))/n_blast_brhs,
               fs_unique_brh = length(intersect(fs_unique, fs_brhs))/n_fs_brhs,
               common_fs_brh = 1-fs_unique_brh,
               blast_unique_brh = length(intersect(blast_unique, blast_brhs))/n_blast_brhs,
               common_blast_brh = 1-blast_unique_brh,
        ) %>% 
        select(-blast_hits, -fs_hits, -fs_brhs, -blast_brhs, 
               -fs_unique, -blast_unique, -common_hits, -overlap)
      
      df_overlap <- bind_rows(df_overlap, to_add)
    }
  }
}


plot_intersection <- df_overlap %>% 
  filter(n_seqs > 1 | n_structs > 1) %>% 
  select(eval_blast, eval_fs, cov, jaccard, overlap_blast, overlap_fs) %>% 
  pivot_longer(cols = !c(eval_blast, eval_fs, cov), names_to = "intersection") %>%
  mutate(cov=paste("cov:", cov), intersection = factor(intersection, levels = c("overlap_blast","jaccard","overlap_fs"))) %>% 
  group_by(eval_blast, eval_fs, intersection, cov) %>% 
  summarise(mean_overlap=mean(value, na.rm = T), sd_overlap=sd(value, na.rm = T)) %>% 
  ggplot(aes(eval_blast, eval_fs)) + 
  geom_tile(aes(fill=mean_overlap), color="black") +
  # geom_text(aes(label=round(mean_overlap, 2)), color="white") +
  # geom_text(aes(label=paste(round(mean_overlap, 2), "±", round(sd_overlap, 2))), color="white") +
  scale_fill_viridis_c(limits = c(0,1), option = "B") +
  facet_grid(cov~intersection) +
  labs(x="Blast E-value", y="Foldseek E-value") + 
  coord_cartesian(expand = 0) +
  theme_minimal() + 
  theme(legend.position = "none")

plot_intersection

plot_brh <- df_overlap %>% 
  filter(n_seqs > 1 | n_structs > 1, cov==30, eval_fs==1e-3, eval_blast==1e-3) %>% 
  select(query, blast_fs_brh, fs_blast_brh) %>% 
  pivot_longer(!query, values_to="overlap") %>% 
  # mutate(name = factor(name, levels=c("blast_fs_brh", "fs_blast_brh"))) %>% 
  ggplot(aes(name, overlap)) + 
  geom_violin(fill="#1B7F79") + 
  theme(axis.text.x = element_text(angle=70, hjust=1), axis.title.x=element_blank())

plot_brh_id <- df_overlap %>% 
  filter(n_seqs > 1 | n_structs > 1, cov==30, eval_fs==1e-3, eval_blast==1e-3) %>% 
  select(query, common_blast_brh, common_fs_brh) %>% 
  pivot_longer(!query, values_to="overlap") %>% 
  # mutate(name = factor(name, levels=c("blast_fs_brh", "fs_blast_brh"))) %>% 
  ggplot(aes(name, overlap)) + 
  geom_violin(fill="#FF4858") + 
  theme(axis.text.x = element_text(angle=70, hjust=1), axis.title.x=element_blank())


final_plot <- (plot_intersection / (plot_brh | plot_brh_id)) + plot_layout(heights = c(3,1))

ggsave(opt$outfile, final_plot, width = 6, height = 8, dpi = 300)


# n_spec <- full_join(count(group_by(blast, target)), count(group_by(fs, target)), by="target", suffix = c("_blast", "_fs"))
# 
# pivot_longer(n_spec, !c(target)) %>% 
#   ggplot(aes(name, value)) + 
#   geom_boxplot()
# 
# 
# ggplot(df_overlap, aes(n_seqs, common_blast_brh)) + 
#   stat_bin_2d() + 
#   ggh4x::facet_nested(rows = vars(cov,eval_blast,eval_fs)) +
#   scale_fill_gradientn(colours = wespal)

# IN CASE YOULL WANT TO EXPLORE CORRELATIONS BETWEEN FS AND BLAST MEASURES
# a <- full_join(select(blast, query, target, pident, length, evalue, bitscore, qcov),
#           select(fs, query, target, pident, length, evalue, bitscore, qcov, prob, lddt, alntmscore),
#           by = c("query", "target"), suffix = c("_blast", "_fs")) %>% 
#   filter(!is.na(pident_blast), !is.na(pident_fs))
# ggplot(a, aes(pident_blast, pident_fs)) +
#   stat_bin_2d(bins = 60) +
#   geom_smooth() + 
#   scale_fill_gradientn(colours = wespal)
# ggplot(a, aes(bitscore_blast, bitscore_fs)) +
#   stat_bin_2d(bins = 60) +
#   geom_smooth() + 
#   scale_fill_gradientn(colours = wespal)
# ggplot(a, aes(evalue_blast, evalue_fs)) +
#   stat_bin_2d(bins = 60) +
#   geom_smooth() + 
#   scale_x_log10() +
#   scale_y_log10() +
#   scale_fill_gradientn(colours = wespal)
# ggplot(a, aes(qcov_blast, qcov_fs)) +
#   stat_bin_2d(bins = 60) +
#   geom_smooth() + 
#   scale_fill_gradientn(colours = wespal)
# 
# a %>% 
#   pivot_longer(cols = c(prob, lddt, alntmscore)) %>% 
#   ggplot(aes(evalue_blast, value)) + 
#   stat_bin_2d(bins = 60) +
#   geom_smooth() + 
#   scale_x_log10() +
#   scale_fill_gradientn(colours = wespal) + 
#   facet_grid(~name)