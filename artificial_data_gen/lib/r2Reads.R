r2Reads <- function(reps, m, mystring, x, y) {
  res = data.frame()
  n = nchar(mystring)
  lowerLim = 1
  higherLim = n + m - 1
  if (reps > 0) {
    randStart = sample(lowerLim:higherLim, reps, replace = TRUE)
    for (j in 1 : reps) {
      temp = randStart[j] - n
      ns = if (temp < 0) (0) else (temp)
      temp = randStart[j] + 1 - m
      end = randStart[j] - ns
      start = if (temp >= 1) (temp) else (1)
      read = paste(strrep("N", ns), stri_reverse(substr(mystring, start, end)), strrep("N", if (temp <=0) (abs(temp) + 1) else (0)), sep="")
      quality <- intToUtf8(round(uniToOthers(x, y, runif(m))))
      res <- rbind(res, data.frame(Sequence = read, Length = nchar(read), Q = quality))
    }
    return(res)
  }else {
    stop("reps argument must be a possitive integer.") 
  }
}