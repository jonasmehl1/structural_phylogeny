suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(tidyverse))
source("workflow/scripts/functions.R")

theme_set(theme_classic())

# order of models
models <- c("QT", "FT", "LG", "GTR", "3Di")


