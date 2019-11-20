rightToLeft_read <- function(H, m, mystring) {
  n = nchar(mystring)
  if (H >= 1 & H <= (n + m - 1)) {
    temp = H - n
    ns = if (temp < 0) (0) else (temp)
    temp = H + 1 - m
    end = H - ns
    start = if (temp >= 1) (temp) else (1)
    read = paste(strrep("N", if (temp <=0) (abs(temp) + 1) else (0)), substr(mystring, start, end), strrep("N", ns), sep="")
    return(read)
  }else {
    return(FALSE)
  }
}