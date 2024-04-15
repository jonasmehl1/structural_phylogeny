# library("optparse")
# 
# option_list <- list(
#   make_option(c("-t", "--time"),
#               type = "character", default = NULL,
#               help = "benchmark results"
#   ),
#   make_option(c("-a", "--alnstats"),
#               type = "character", default = NULL,
#               help = "alignment stats"
#   ),
#   make_option(c("-o", "--output"),
#               type = "character", default = NULL,
#               help = "output plot", metavar = "character", dest = "outfile"
#   )
# )
# 
# desc <- "Produces informative plot on trees comparison"
# 
# opt_parser <- OptionParser(option_list = option_list, description = desc)
# opt <- parse_args(opt_parser)

# opt <- NULL
# opt$alnstats <- "results/Hsap_opistho/stats/UP000005640_aln.stats"
# opt$time <- "results/Hsap_opistho/stats/UP000005640_runtime.stats"
# opt$outfile <- "test/stats/test.png"

library(tidyverse)
library(patchwork)
source("workflow/scripts/functions.R")

theme_set(theme_bw(base_size = 8))

cnames <- c("bn","step","seed","gene","method","alphabet",
            "s","h:m:s","max_rss","max_vms","max_uss","max_pss",
            "io_in","io_out","mean_load","cpu_time")

df <- read_delim(snakemake@input[["time"]], col_names = cnames) %>% 
  mutate(max_rss_gb = measurements::conv_unit(max_rss, from = "MB", to = "GB"))

aln_stats <- read_delim(snakemake@input[["aln"]]) %>%
  mutate(trim = ifelse(grepl("clean", file), "clean", "aln"),
         file = gsub("\\..*", "", file),
         gap_prop = sum_gap/(num_seqs*min_len)) %>% 
  separate(file, into = c("gene", "method", "alphabet"), "_")

# is blast slower than fs
plot_homo <- filter(df, step=="homology", is.na(method)) %>% 
  ggplot(aes(gene, y=`h:m:s`, fill=gene)) + 
  geom_point(size=4, pch=21) + 
  labs(subtitle = "Homology search", y = "Time", x="Method") + 
  scale_fill_manual(values = palettes_method) + 
  theme(legend.position = "none")

# is the alignment different?
plot_aln <- filter(df, step=="mafft") %>%
  left_join(filter(aln_stats, trim=="aln")) %>% 
  ggplot(aes(log(num_seqs*min_len), `h:m:s`, shape=method, fill=alphabet)) +
  geom_point(pch=21) + 
  # geom_smooth(se=F, aes(group=method)) +
  labs(y="Time", subtitle = "Mafft --auto") +
  scale_fill_manual(values = palettes_alphabet) + 
  theme(legend.position = "bottom")

# How does trimming change based on alphabets
# aln_stats %>%
#   ggplot(aes(trim, gap_prop, fill=alphabet)) +
#   # geom_boxplot() + 
#   geom_violin() +
#   scale_fill_manual(values=palettes_alphabet)

# is iqtree slower than ft?
plot_memtree <- filter(df, step %in% c("ft", "iqtree")) %>%
  ggplot(aes(max_rss_gb, `h:m:s`, fill=step)) + 
  geom_point(pch=21) + 
  labs(y="Time", x="Mem (GB)", subtitle="Foldtree vs IQTree") + 
  theme(legend.position = "bottom")

# if so how much does it change between alphabets?
plot_iqtree <- filter(df, step=="iqtree") %>%
  left_join(filter(aln_stats, trim=="clean")) %>% 
  ggplot(aes(log(min_len*num_seqs), `h:m:s`, fill=alphabet)) + 
  geom_point(pch=21) +
  geom_smooth(se=F) +
  labs(y="Time", subtitle = "IQtree runtimes") + 
  scale_fill_manual(values = palettes_alphabet) + 
  facet_grid(~method) +
  theme(legend.position = "bottom")


final_plot <- plot_homo + plot_aln + plot_memtree + plot_iqtree + 
  plot_layout(widths = c(1, 3)) + 
  plot_annotation(tag_levels = c("A"))


ggsave(snakemake@output[[1]], final_plot, width = 9, height = 6, dpi = 300)

