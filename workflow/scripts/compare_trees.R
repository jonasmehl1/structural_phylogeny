suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggpubr))
source("workflow/scripts/functions.R")

theme_set(theme_classic())

# order of models
models <- c("QT", "FT", "LG", "GTR", "3Di")

# all trees
df_trees <- read_delim(snakemake@input[["trees"]], show_col_types = FALSE,
                       col_names = c("gene", "target", "alphabet", "model", "tree")) %>% 
  mutate(model=factor(model, levels=models))

# ML trees with selected model and LL
df_ml <- read_delim(snakemake@input[["mltrees"]], delim = "\t", 
                    show_col_types=FALSE,
                    col_names = c("gene", "rate", "ll", "tree", "targets")) %>% 
  # mutate(fn=gsub("_mltrees.txt", "", basename(fn))) %>% 
  # separate(fn, c("gene", "target", "alphabet", "model"), sep = "_") %>% 
  mutate(model=gsub("GTR20", "GTR", gsub("\\+.*", "", rate)), 
         rate=sub("\\+", "", gsub("GTR20|3DI|LG", "", rate)))

# Read all trees
ts <- read.tree(text = df_trees$tree)
names(ts) <- paste(df_trees$gene, df_trees$target, df_trees$alphabet, df_trees$model, sep = "_")

# # get BS distribution
df_bs <- fortify(ts) %>%
  filter(!isTip) %>%
  mutate(support=as.numeric(label)) %>%
  separate(.id, c("gene", "target", "alphabet", "model")) %>%
  filter(!is.na(support))

# plot bootstrap
plot_bs <- ggplot(df_bs, aes(model, support, fill=target)) +
  # geom_violin(position = position_dodge(width = 1)) +
  geom_boxplot() +
  # ggdist::stat_pointinterval(position = position_dodge(width = 1), color="black") +
  scale_fill_manual(values = palettes_method) +
  theme_classic() +
  labs(y="Boostrap", x="") +
  theme(legend.position = "bottom")

# First of all let's check if they have the same heterogeneity model
plot_rate <- df_ml %>%
  select(-tree) %>%
  mutate(Freq=grepl("F", rate), Inv=grepl("I", rate),
         Gamma=grepl("G", rate), FreeRate=grepl("R", rate)) %>%
  select(-rate, -ll) %>%
  pivot_longer(cols = c(Freq,Inv,Gamma,FreeRate)) %>%
  inner_join(.,., by=c("gene", "targets", "name")) %>%
  filter(model.x>model.y) %>%
  mutate(comp=paste0(model.x, " vs ", model.y)) %>%
  group_by(name, model.x, model.y) %>%
  count(value.x, value.y) %>%
  ggplot(aes(value.x, value.y, fill=n)) +
  geom_tile() +
  geom_text(aes(label=n), color="white") +
  scale_fill_viridis_c(option = "C") +
  ggh4x::facet_nested(name~model.x+model.y) +
  coord_fixed(expand = F) +
  labs(x="Is same?") +
  theme(legend.position = "none", axis.title.y = element_blank())


plot_3di <- filter(df_ml, model!="LG") %>% 
  pivot_wider(id_cols = c(gene, targets), values_from = c(rate,ll), names_from = model) %>% 
  # print() %>% 
  mutate(diff_ll=ll_3DI-ll_GTR, better_3Di=diff_ll>0) %>% 
  ggplot(aes(diff_ll)) +
  # geom_density() + 
  ggdist::stat_slab(aes(fill = after_stat(level)), .width = c(.5, .66, .95, 1)) +
  scale_fill_brewer() +
  geom_text(data = . %>%
              mutate(x=quantile(diff_ll, .95)) %>% 
              group_by(targets, better_3Di, x) %>% 
              summarise(n=n()) %>% 
              pivot_wider(names_from = better_3Di, values_from = n) %>% 
              mutate(prop_better=`TRUE`/(`TRUE`+`FALSE`)), 
            mapping = aes(x=x, y=0.5, 
                          label=paste0(round(prop_better*100, 2),"%"))) +
  geom_vline(xintercept = 0) + 
  coord_cartesian(expand = 0) + 
  facet_wrap(~targets, nrow = 1) +
  labs(x="LL 3Di - LL GTR", y="Density", subtitle = "GTR vs 3Di LogLik") + 
  theme(legend.position = "bottom")

# Compare trees with RF distance
rf_df <- filter(df_trees, target=="common") %>% 
  select(-alphabet) %>% 
  inner_join(x=., y = ., by = c("gene", "target")) %>%
  filter(!is.na(tree.x), !is.na(tree.y)) %>%
  rowwise() %>%
  mutate(RF = phangorn::RF.dist(normalize = TRUE, read.tree(text = tree.x), read.tree(text=tree.y)),
         n_tips = length(read.tree(text = tree.x)$tip.label)) %>%
  select(-tree.x, -tree.y)

rf_text <- group_by(rf_df, model.x, model.y) %>% 
  filter(as.integer(model.x)<=as.integer(model.y)) %>% 
  summarise(value=median(RF, na.rm=T))

# plto the heatmap of dist distribution and median value
rf_plot <- rf_df %>% 
  ggplot(aes(RF)) +
  geom_tile(aes(x=0.5, y=0.5, fill=value), data = rf_text) +
  geom_text(aes(x=0.5, y=.5, label=round(value, 2)), color="white", data = rf_text) +
  geom_rug(linewidth=.1, data= . %>% filter(as.integer(model.x)>as.integer(model.y))) +
  geom_density(aes(y = after_stat(scaled)), data= . %>% filter(as.integer(model.x)>as.integer(model.y))) + 
  scale_fill_distiller(palette="YlGnBu", direction = 1, limits = c(0,1)) + 
  scale_x_continuous(breaks = c(0,.5,1)) +
  coord_fixed(xlim = c(0,1), expand=0) +
  # facet_grid(method~method2, switch = "both") +
  facet_grid(model.x~model.y) +
  # ggh4x::facet_nested(targets+model~targets2+model2, nest_line = element_line(linetype = 1)) +
  labs(x="RF Distance") +
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
scores <- read_delim(snakemake@input[["scores"]], show_col_types = FALSE) %>% 
  mutate(model=factor(model, levels=models))

reco <- read_delim(snakemake@input[["reco"]], show_col_types = FALSE)

scores <- scores %>% left_join(reco, by=c("id"="gene"))
print(scores)

# Ranger results
ranger_plot <- mutate(scores, n_events = (dups+losses)/n_tips) %>% 
  ggplot(aes(model, n_events, fill=targets)) + 
  geom_boxplot() +
  scale_fill_manual(values = palettes_method) +
  labs(x="", y="(D+L)/# Tips") +
  guides(fill = guide_legend(nrow = 3)) + 
  theme(legend.position = "none")

TCS_df <- scores %>% 
  select(id, targets, model, score)
  # pivot_wider(names_from = model, values_from = score)

plot_tcs <- inner_join(TCS_df, TCS_df, by = c("id", "targets"), 
                       multiple = "all", relationship = "many-to-many") %>% 
  filter(as.integer(model.x)>as.integer(model.y)) %>% 
  ggplot(aes(score.x, score.y, color=targets)) +
  geom_point(alpha=.4, size=1) + 
  geom_abline() +
  ggh4x::facet_grid2(model.x~model.y, render_empty = F, 
                     strip = ggh4x::strip_vanilla()) + 
  scale_color_manual(values = palettes_method) +
  labs(x="TCS Score", y="") +
  coord_fixed() +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))


plot_varr2t <- scores %>% 
  filter(model %in% c("3Di", "LG", "GTR")) %>% 
  ggplot(aes(variance_r2t, color=model)) +
  geom_density() + 
  facet_wrap(.~targets, nrow = 1) +
  scale_x_continuous(limits = c(0,.25), expand = c(0,0.01)) +
  scale_color_manual(values = palettes_model) +
  theme(legend.position = 'bottom')

reco <- ((plot_bs/ranger_plot)|plot_tcs) + plot_layout(widths = c(1,2.5))

model <- plot_3di/(plot_rate|rf_plot)/plot_varr2t + 
  plot_layout(heights = c(1,4,1))

ggsave(snakemake@output[["reco"]], reco, width = 10, height = 8)
ggsave(snakemake@output[["model"]], model, width = 15, height = 12)

