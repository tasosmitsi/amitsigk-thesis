perPosOfSpeChrom <- function(data) {
  res <- data.frame()
  for(chrom in levels(data$CHROM)) {
    res <- rbind(res, data.frame(PerPosOfSpeChrom = length(which(final_res$CHROM == chrom & final_res$Positive)) / length(which(final_res$CHROM == chrom))))
  }
  row.names(res) <- levels(data$CHROM)
  return(res)
}