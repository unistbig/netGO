BuildCutoffPPI <- function(PPI, Cutoff = 0.75) {
  for (i in 1:nrow(PPI)) {
    PPI[i, ] <- sapply(PPI[i, ], function(j) {
      if (j < Cutoff) {
        return(0)
      }
      return(j)
    })
  }
  rem <- which(rowSums(PPI) == 0)
  PPI <- PPI[-rem, -rem]
  PPI
}
