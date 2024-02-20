library("optparse")

option_list <- list(
  make_option(c("-s","--sptree"),
              type = "character", default = NULL,
              help = "species tree"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = "output filename", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

library(phytools)

sptree <- read.tree(opt$sptree)
sptree$edge.length <- rep(0.1, length(sptree$edge))

if (!dir.exists(dirname(opt$output))) {
  dir.create(dirname(opt$output), recursive = TRUE)
}

write.tree(sptree, opt$output)
