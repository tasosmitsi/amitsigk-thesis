leftToRight_read <- function(L, m, mystring) {
  n = nchar(mystring)
  if (L >= (-m+2) & L <= (n)) {
    if (L < 1) {
      ns = abs(L) + 1
      read = paste(strrep("N", ns), substr(mystring, 1, m - ns), strrep("N", if ((m - ns - n) > 0) (m - ns - n)  else 0), sep="")
    }else if (L >= 1) {
      ns = m - (n - L + 1)
      read = paste(substr(mystring, L, n + ns), strrep("N", if (ns < 0) 0 else ns), sep="")
    }
    return(read)
  }else {
    return(FALSE)
  }
}