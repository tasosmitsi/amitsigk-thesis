findRelations <- function(templateRef, templateAlt, vfcRef, vfcAlt, posIndexes) {
  res <- data.frame()
  if (length(posIndexes) == length(vfcRef)) {
    for (i in 1 : length(vfcRef)) {
      for (posIndex in posIndexes[[i]]) {
        if (templateRef[posIndex] == vfcRef[i] & templateAlt[posIndex] == vfcAlt[i]) {
          #do something to store matching
          res <- rbind(res, data.frame(templateId = posIndex, vfcId = i))
          break
        }
      }
    }
  }
  return(res)
}