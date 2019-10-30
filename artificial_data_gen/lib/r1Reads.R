r1Reads <- function(reps, m, mystring, x, y) {
  res = data.frame()
  n = nchar(mystring)
  lowerLim = -m+2
  higherLim = n
  if (reps > 0) {
    randStart = sample(lowerLim:higherLim, reps, replace = TRUE)
    for (j in 1 : reps) {
      if (randStart[j] < 1) {
        ns = abs(randStart[j]) + 1
        read = paste(strrep("N", ns), substr(mystring, 1, m - ns), strrep("N", if ((m - ns - n) > 0) (m - ns - n)  else 0), sep="")
      }else if (randStart[j] >= 1) {
        ns = m - (n - randStart[j] + 1)
        read = paste(substr(mystring, randStart[j], n + ns), strrep("N", if (ns < 0) 0 else ns), sep="")
      }
      quality <- intToUtf8(round(uniToOthers(x, y, runif(m))))
      res <- rbind(res, data.frame(Sequence = read, Length = nchar(read), Q = quality))
    }
    return(res)
  }else {
    stop("reps argument must be a possitive integer.") 
  }
}