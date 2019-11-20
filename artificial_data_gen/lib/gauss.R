gauss <- function(x, mean, sd) {
  y <- (1/(sqrt(2*pi*sd^2))) * exp(-((x - mean)^2)/(2*sd^2))
  return(y)
}