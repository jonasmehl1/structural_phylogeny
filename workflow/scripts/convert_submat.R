blosum <- read.delim("test/msa/BLOSUM62", comment.char = "#", sep = "")
blosum <- blosum[1:20, 1:20]

# Pij <- exp(lambda*Sij)*Fi*Fj
# Sij <- (1/lambda)*ln(Pij/(Fi*Fj))
# Sij <- round(Sij)
# as its rounded its impossible to go from S to P!

freqs <- c(0.074,0.052,0.045,0.054,0.025,0.034,0.054,0.074,0.026,0.068,0.099,0.058,0.025,0.047,0.039,0.057,0.051,0.013,0.032,0.073)
names(freqs) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

lambda <- 0.347

outmat <- matrix(nrow = 20, ncol = 20)
rownames(outmat) <- rownames(blosum)
colnames(outmat) <- colnames(blosum)
row <- 18
col <- 18
for(row in 1:nrow(blosum)) {
  fi <- freqs[rownames(blosum)[row]]
    for(col in 1:ncol(blosum)) {
      fj <- freqs[colnames(blosum)[col]]
      Sij <- blosum[row,col]
      Pij <- exp(lambda*Sij)*fi*fj
      outmat[row, col] <- Pij
  }
}

scale(outmat)


x <- scan("test/msa/blosum62.parse")
dims <- floor(sqrt(length(x) * 2))

m <- matrix(NA, 20, 20)
m[upper.tri(m, diag = FALSE)] <- x
m <- t(m)
rownames(m) <- rownames(blosum)
colnames(m) <- colnames(blosum)

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

m <- makeSymm(m)
diag(m) <- 1

