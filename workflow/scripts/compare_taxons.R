library("optparse")

option_list <- list(
  make_option(c("-i", "--input"),
              type = "character", default = NULL,
              help = "input table file", dest = "input"
  ),
  make_option(c("-m", "--meta"),
              type = "character", default = NULL,
              help = "meta file", dest = "meta"
  ),
  # make_option(c("-r", "--ref"),
  #             type = "character", default = NULL,
  #             help = "reference proteome seeds", dest = "ref"
  # ),
  make_option(c("-t", "--taxidmap"),
              type = "character", default = NULL,
              help = "taxidmap file", dest = "taxidmap"
  ),
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

desc <- "Produces informative plot on foldseek and blast results taxonomic distribution"

opt_parser <- OptionParser(option_list = option_list, description = desc)
opts <- parse_args(opt_parser)

library(tidyverse)
library(patchwork)

source("scripts/functions.R")
table_columns <- c('Proteome_ID', 'Tax_ID', 'count1', 'count2', 'count3', 
                   'genome', 'source', 'species', 'mnemo')
theme_set(theme_bw())

# opts <- NULL
# opts$ref <- "data/ids/UP000005640_rev.txt"
# opts$taxidmap <- "results/Hsap_LUCA/db/taxidmap"
# opts$blast <- "test/blast/blast.tsv"
# opts$blast_brh <- "results/Hsap_LUCA/homology/UP000005640_blast_brh.tsv"
# opts$fs <- "test/blast/fs.tsv"
# opts$fs_brh <- "results/Hsap_LUCA/homology/UP000005640_fs_brh.tsv"
# opts$seed <- 9606
# opts$meta <- "data/meta/Hsap_taxon.csv"
# opts$input <- "data/input_tables/LUCA.csv"

# in human only consider reference proteins 20k
# ref <- readLines(opts$ref)
# ref <- paste0("AF-",ref,"-F1")

taxidmap <- read_delim(opts$taxidmap, col_names = c("target", "Tax_ID"), show_col_types=F)

# blast results
bs <- read_delim(opts$blast, col_names=blast_columns, show_col_types=F) %>% 
  filter(query!=target)
bs_brh <- read_delim(opts$blast_brh, show_col_types=F,
                     col_names = c("query", "target"))

# foldseek results
fs <- read_delim(opts$fs, col_names=fs_columns, show_col_types=F) %>% 
  left_join(taxidmap) %>% 
  filter(query!=target)
fs_brh <- read_delim(opts$fs_brh, show_col_types=F, 
                     col_names = c("query", "target"))

# read taxonomic data
tax <- read_delim(opts$meta, show_col_types=F)
# read input table
lvls_tax <- rev(c(unique(tax$Clade), "archaea", "bacteria"))

table <- read_delim(opts$input, show_col_types=F, 
                    delim = "\t", col_names = table_columns) %>% 
  left_join(tax, by = "Tax_ID") %>% 
  mutate(Clade = ifelse(is.na(Clade), mnemo, Clade),
         Clade = factor(Clade, levels = lvls_tax))

evals <- c(Inf, 1e-2, 1e-3, 1e-5)
covs <- c(0, 30, 80)

final_df <- tibble()

for (eval in evals) {
  for (cov in covs) {
    print(paste(eval, "-",cov))
    add_bs <- bs %>%
      filter(evalue < eval, qcov > cov) %>% 
      left_join(table, by = c("staxids"="Tax_ID")) %>% 
      group_by(query, Clade) %>% 
      count() %>% 
      group_by(Clade) %>% 
      count(n>0) %>% 
      mutate(n=n/length(unique(bs$query)), evalue=eval, cov=cov, method="blast", brh=FALSE)
    add_bs_brh <- bs %>%
      filter(evalue < eval, qcov > cov) %>% 
      inner_join(bs_brh, by = c("query", "target")) %>% 
      left_join(table, by = c("staxids"="Tax_ID")) %>% 
      group_by(query, Clade) %>% 
      count() %>% 
      group_by(Clade) %>% 
      count(n>0) %>% 
      mutate(n=n/length(unique(bs$query)), evalue=eval, cov=cov, method="blast", brh=TRUE)
    add_fs <- fs %>%
      filter(evalue < eval, qcov > cov/100) %>% 
      left_join(table, by="Tax_ID") %>% 
      group_by(query, Clade) %>% 
      count() %>% 
      group_by(Clade) %>% 
      count(n>0) %>% 
      mutate(n=n/length(unique(fs$query)), evalue=eval, cov=cov, method="fs", brh=FALSE)
    add_fs_brh <- fs %>%
      filter(evalue < eval, qcov > cov/100) %>% 
      inner_join(fs_brh, by = c("query", "target")) %>% 
      left_join(table, by="Tax_ID") %>% 
      group_by(query, Clade) %>% 
      count() %>% 
      group_by(Clade) %>% 
      count(n>0) %>% 
      mutate(n=n/length(unique(fs$query)), evalue=eval, cov=cov, method="fs", brh=TRUE)
    final_df <- bind_rows(final_df, add_bs, add_bs_brh, add_fs, add_fs_brh)
  }
}

p1 <- ggplot(filter(final_df, !brh)) + 
  geom_tile(aes(y=Clade, method, fill=n)) + 
  # facet_grid(evalue~cov) +
  ggh4x::facet_nested(~evalue+cov) +
  scale_fill_gradientn(colours = wespal) + 
  ggtitle("All hits") +
  coord_cartesian(expand=0)

p2 <- ggplot(filter(final_df, brh)) + 
  geom_tile(aes(y=Clade, method, fill=n)) + 
  # facet_grid(evalue~cov) +
  ggh4x::facet_nested(~evalue+cov) +
  scale_fill_gradientn(colours = wespal) + 
  ggtitle("Only BRH") +
  coord_cartesian(expand=0)

final_plot <- p1/p2

ggsave(opts$outfile, final_plot, width=10, height=7)
