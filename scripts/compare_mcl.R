blines <- readLines("results/mcl/blast.I12.test")
blines <- blines[grepl("\\t", blines)]
flines <- readLines("results/mcl/foldseek.I12.test")
flines <- flines[grepl("\\t", flines)]

list_b <- lapply(blines[1:1000], function(x) strsplit(x, "\\t")[[1]])
list_f <- lapply(flines[1:1000], function(x) strsplit(x, "\\t")[[1]])

# combns <- gtools::combinations(length(list_b), length(list_f), )
nrows <- length(list_b)
ncols <- length(list_f)

result_matrix <- matrix(0, nrow = nrows, ncol = ncols)

for (i in 1:nrows){
    set_blast <- list_b[i][[1]]
  for (j in 1:ncols){
    set_fs <- list_f[j][[1]]
    norm_value <- pmax(length(set_blast), length(set_fs))
    result_matrix[i,j] <- length(intersect(set_blast, set_fs))/norm_value
  }
}


heatmap(result_matrix)

