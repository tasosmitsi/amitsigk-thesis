findPosIndexes <- function(vcf_pos, startPos, endPos) {
  return(which(startPos <= vcf_pos & vcf_pos <= endPos))
}