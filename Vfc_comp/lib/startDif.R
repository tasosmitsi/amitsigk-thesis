startDif <- function(string1, string2) {
  n = pmin(nchar(string1), nchar(string2))
  i = 1
  while (i <= n) {
    if (substr(string1, i, i) != substr(string2, i, i)) 
      return(i)
    i = i + 1
  }
}