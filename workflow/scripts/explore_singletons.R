queries <- sample(unique(blast$query), 20)

fs_filt <- filter(fs, evalue<.01, query %in% queries, query!=target)
blast_filt <- filter(blast, evalue<.01, query %in% queries, query!=target)

blast_self <- read_delim("results/Ncra_opistho/homology/allvall/UP000001805_UP000001805_blast.tsv", col_names = blast_columns) %>% 
  filter(query==target) %>% 
  mutate(blast_selfbit=bitscore) %>% 
  select(query, blast_selfbit)

fs_self <- read_delim("results/Ncra_opistho/homology/allvall/UP000001805_UP000001805_fs.tsv", col_names = fs_columns) %>% 
  filter(query==target) %>% 
  mutate(fs_selfbit=bitscore) %>% 
  select(query, fs_selfbit)

self <- full_join(blast_self, fs_self)

df <- full_join(select(blast_filt, query, target, evalue, bitscore), 
          select(fs_filt, query, target, evalue, bitscore), 
          by = c("query", "target"), suffix=c("_blast", "_fs")) %>% 
  left_join(self) %>% 
  mutate(singleton = case_when(is.na(evalue_fs) ~ "only_blast", 
                               is.na(evalue_blast) ~ "only_fs", .default = "common"),
         blast_norm = bitscore_blast/blast_selfbit, fs_norm = bitscore_fs/fs_selfbit,
         logevalue_blast = -log10(evalue_blast), logevalue_fs = -log10(evalue_fs)) %>% 
  mutate(logevalue_blast = ifelse(is.infinite(logevalue_blast), 180, logevalue_blast),
         logevalue_fs = ifelse(is.infinite(logevalue_fs), 180, logevalue_fs))

bit_plot <- df %>% 
  select(query, target, blast_norm, fs_norm, singleton) %>% 
  pivot_longer(!c(query, target, singleton), values_to = "normalized_bitscore") %>% 
  ggplot(aes(x=normalized_bitscore, y=query, fill=singleton)) + 
  ggridges::geom_density_ridges(scale = 1, alpha=.4, quantile_lines = TRUE, quantiles = 2) + 
  facet_grid(~name)

bit_plot


# test for each gene if distribution is different
# look at % identity of foldseek singletons