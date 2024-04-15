blast_columns <- c("query","target","pident","length","mismatch","gapopen",
                   "qstart","qend","sstart","send","evalue","bitscore",
                   "qcov","qcovhsp","qlen","slen","staxids")

fs_columns <- c("query","target","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore",
                "lddt","alntmscore","rmsd","prob","qcov","tcov")

palettes_method <- c("#267365", "#F29F05", "#F23030")
names(palettes_method) <- c("common", "blast", "fs")

palettes_alphabet <- c("#327AD9", "#BF3475", "#F29F05")
names(palettes_alphabet) <- c("3Di", "aa", "fident")
wespal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")


add_marginal <- function(df, plot, x, y, fill, palette){
  xdens <- axis_canvas(plot, axis = "x")+
    geom_density(data = df, aes(x = .data[[x]], fill = .data[[fill]]),
                 alpha = 0.7, size = 0.2) +
    coord_cartesian(expand=0) + 
    scale_fill_manual(values = palette)
  ydens <- axis_canvas(plot, axis = "y", coord_flip = TRUE) +
    geom_density(data = df, aes(x = .data[[y]], fill = .data[[fill]]),
                 alpha = 0.7, size = 0.2)+
    coord_flip(expand = 0) +
    scale_fill_manual(values = palette)
  
  p1 <- insert_xaxis_grob(plot, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  ggdraw(p2)
}

