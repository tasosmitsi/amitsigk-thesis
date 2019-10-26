uniToOthers <- function(x, y, rndy) {
  res = c()
  for(p in rndy){
      idx = which.min(abs(y -  p))
      xrand = x[idx]
      #yrand = y[idx]
      res = c(res, xrand)
      #points(xrand, yrand, pch = 25, col = "blue")
  }
  return(res)
}