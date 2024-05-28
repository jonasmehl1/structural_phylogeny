interesting_types <- c("Topological domain", "Domain", 
                       "Region", "Helix", "Zinc finger",
                       "Beta strand", "Turn", "Coiled coil", "Transmembrane")

blast_columns <- c("query","target","pident","length","mismatch","gapopen",
                   "qstart","qend","sstart","send","evalue","bitscore",
                   "qcov","qcovhsp","qlen","slen","staxids")

fs_columns <- c("query","target","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore",
                "lddt","alntmscore","rmsd","prob","qcov","tcov")

models <- c("LG", "3Di", "GTR", "QT", "FT", "FTPY")

palettes_model <- c("#BF360C", "#75C323", "#40A79D", "#FFA000", "#D4E157","#646B00")
names(palettes_model) <- models
palettes_model_d <- colorspace::darken(palettes_model, amount = 0.6)
names(palettes_model_d) <- models

palettes_method <- c("#F23030", "#F29F05", "#267365", "#D6FAB7")
names(palettes_method) <- c("common", "blast", "fs", "union")

palette_singleton <- c("#F23030", "#F29F05", "#267365", "#D6FAB7")
names(palette_singleton) <- c("common", "only_blast", "only_fs", "union")

palette_brh <- c("#F23030","#F29F05", "#267365")
names(palette_brh) <- c("common","blast_brh", "fs_brh")

palettes_alphabet <- c("#327AD9", "#BF3475", "#F29F05")
names(palettes_alphabet) <- c("3Di", "aa", "fident")
wespal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

palettes_step <- c( "#F0E442", "#009E73", "#CC79A7")
names(palettes_step) <- c("foldtree", "iqtree", "quicktree")

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


# Function to find identical columns and return their indexes
find_identical_columns <- function(matrix1, matrix2) {
  identical_columns <- apply(matrix1, 2, function(col1) any(apply(matrix2, 2, function(col2) all(col1 == col2))))
  return(which(identical_columns))
}


get_idxs_df <- function(seqs, aln, trimmed) {
  kept_columns <- find_identical_columns(as.matrix(aln@unmasked), as.matrix(trimmed@unmasked)) - 1
  df_idxs <- NULL
  for (seq in names(seqs)) {
    str_untrimmed <- as.character(aln@unmasked[[seq]])
    gaps <- stringr::str_locate_all(str_untrimmed, '-')[[1]]
    idxs <- c()
    aln_idxs <- c()
    
    for (idx in 1:nchar(str_untrimmed)){
      if (!idx %in% gaps[,1]) {
        to_remove <- sum(gaps[,1]<idx)
        aln_idxs <- c(aln_idxs, idx)
        idxs <- c(idxs, idx-to_remove)
      }
    }
    
    df_seed <- tibble(label=seq, pos=idxs, pos_aln=aln_idxs) %>% 
      left_join(tibble(pos_aln=kept_columns, pos_trim=1:length(kept_columns)), by="pos_aln") %>% 
      mutate(is_in_trimmed=pos_aln %in% kept_columns)
    
    df_idxs <- bind_rows(df_idxs, df_seed)
  }
  return(df_idxs)
}
