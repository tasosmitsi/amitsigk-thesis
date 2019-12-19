changeIndelsFormat <- function(ref, alt, type) {
  if (type == "insertion") {
    ref <- "_"
    alt <- str_extract(alt, "\\B[AGTC]*")
  }else if (type == "deletion") {
    ref <- str_extract(ref, "\\B[AGTC]*")
    alt <- "_"
  }
  return(data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE))
  #return(t(c(ref,alt)))
}