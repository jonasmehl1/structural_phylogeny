df <- read_delim("trees_sampling.txt", col_names = F)


ggplot(df, aes(X3, X4, color=X2/pmin(X3,X4))) + 
  geom_point() + 
  coord_fixed() +
  theme_bw() +
  geom_abline(slope = 1) +
  scale_color_viridis_c()
