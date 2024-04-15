source("scripts/functions.R")
theme_set(theme_bw())
opt$blast <- "results/Hsap_opistho//homology/UP000005640_blast.tsv"
opt$fs <- "results/Hsap_opistho/homology/UP000005640_fs.tsv"

fs <- read_delim(opt$fs, col_names = fs_columns, show_col_types = FALSE) %>% 
  mutate(evalue = ifelse(evalue==0, 1e-180, evalue)) %>%
  filter(query!=target) %>% 
  filter(evalue<1e-3)
  # left_join(taxidmap)
blast <- read_delim(opt$blast, col_names = blast_columns, show_col_types = FALSE) %>% 
  mutate(evalue = ifelse(evalue==0, 1e-180, evalue)) %>%
  filter(query!=target) %>% 
  filter(evalue<1e-3)

df_overlap <- tibble()

min_hits <- c(seq(50, 900, 100), 1000)

for (hits in min_hits){
      a <- group_by(blast, query) %>% 
        arrange(desc(bitscore)) %>% 
        slice_head(n = hits)
      b <- group_by(fs, query) %>% 
        arrange(desc(bitscore)) %>% 
        slice_head(n = hits)

      to_add <- full_join(select(a, query, target, pident, evalue, bitscore, qcov),
                          select(b, query, target, pident, evalue, bitscore, qcov),
                          by = c("query", "target"), suffix=c("_blast", "_fs")) %>% 
        mutate(singleton = case_when(is.na(evalue_fs) ~ "only_blast",
                                     is.na(evalue_blast) ~ "only_fs", .default = "common")) %>%
        group_by(query, singleton) %>% 
        count() %>%
        ungroup(singleton) %>% 
        mutate(n_tot = n/sum(n), n_hits = hits)

      df_overlap <- bind_rows(df_overlap, to_add)
}


df_overlap %>%
  group_by(n_hits, singleton) %>%
  summarize(mn=median(n_tot)) %>%
  ggplot(aes(n_hits, mn, color=singleton)) + 
  geom_point() + 
  # geom_smooth(method = lm, formula = y ~ splines::bs(x, 3)) +
  geom_vline(xintercept = 150) +
  geom_line(aes(group=singleton))
  # scale_y_continuous(limits = c(0,1))
