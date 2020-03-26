findIndelsType <- function(ref, alt) {
  if (nchar(toString(ref)) == nchar(toString(alt))) {
    return ("change")
  }else if (nchar(toString(ref)) > nchar(toString(alt))) {
    return ("deletion")
  }else if (nchar(toString(ref)) < nchar(toString(alt))) {
    return ("insertion")
  }
}