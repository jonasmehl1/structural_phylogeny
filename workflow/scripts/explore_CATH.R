library(tidyverse)
source("workflow/scripts/functions.R")
theme_set(theme_classic())

CATH_df <- read_delim(list.files("results/data/cath/", full.names = TRUE), delim = "\t", id = "proteome") %>% 
  filter(!is.na(Gene3D)) %>% 
  mutate(proteome=gsub("_cath.tsv", "", basename(proteome)), 
         n_domains = str_count(Gene3D, ";")) %>% 
  # filter(n_domains==1) %>% 
  mutate(#Gene3D = gsub(";", "", Gene3D),
         target=paste0("AF-", Entry, "-F1"))

nms <- read_delim("test/cath/cath-names.txt", comment = "#", delim = "    ",
                  col_names = c("Gene3D", "abbr", "name")) %>% 
  select(-abbr) %>% 
  mutate(name=gsub("^:", "", name))

ggplot(CATH_df, aes(n_domains)) + 
  geom_histogram(binwidth = 1, position = "dodge") +
  scale_x_continuous(breaks = seq(1,100,1))

# you can join filtered blast and foldseek information (the ones that reach the trees)
# Then see the distinct number of CATH, CAT, CA and C

blast <- read_delim(readLines("test/cath/testblast.files"), col_names = blast_columns)
fs <- read_delim(readLines("test/cath/testfs.files"), col_names = fs_columns)

df <- full_join(mutate(select(fs, query, target), method="fs"), mutate(select(blast, query, target), method="blast"),
          by = c("query", "target"), suffix=c("_fs", "_blast")) %>% 
  mutate(singleton = case_when(is.na(method_fs) ~ "only_blast",
                        is.na(method_blast) ~ "only_fs", 
                        .default = "common"))

inner_join(df, select(CATH_df, -proteome, -Entry)) %>%
  group_by(query, singleton, n_domains) %>% 
  summarise(n_distinct_CATH = n_distinct(Gene3D)) %>% 
  filter(n_domains==1) %>% 
  ggplot(aes(n_distinct_CATH, fill=singleton)) + 
  geom_histogram(binwidth = 1, position = "dodge") +
  # facet_grid(n_domains~., scales = "free_y") +
  scale_x_continuous(breaks = seq(1,100,1)) + 
  scale_y_continuous(expand = expansion(0,0)) + 
  scale_fill_manual(values=palette_singleton)

# stratify by C A T H

# Try to deal with multidomains

long_CATH_df <- CATH_df %>% 
  # filter(target %in% unique(blast$query)) %>% 
  mutate(Gene3D=strsplit(Gene3D, ";")) %>% 
  unnest(Gene3D) %>% 
  select(-proteome, -Entry) %>% 
  separate(col = Gene3D, into = c("C", "A", "T", "H"), sep = "\\.", remove = FALSE) %>% 
  mutate(CA = paste(C,A, sep = "."), CAT=paste(C,A,T,sep = ".")) %>% 
  left_join(nms, by=c("C"="Gene3D")) %>%  
  left_join(nms, by=c("CA"="Gene3D"), suffix = c("_C", "_CA")) %>% 
  left_join(nms, by=c("CAT"="Gene3D"), suffix = c("_CA", "_CAT"))

df %>% 
  select(-contains("method")) %>%
  left_join(long_CATH_df, by=c("query"="target")) %>% 
  left_join(long_CATH_df, by=c("target")) %>% 
  filter(!is.na(Gene3D.x), !is.na(Gene3D.y), query!=target)

# how to treat NA values?