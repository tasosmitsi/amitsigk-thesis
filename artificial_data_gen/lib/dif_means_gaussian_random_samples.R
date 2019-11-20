dif_means_gaussian_random_samples <- function(L_q, H_q, means, sds) {
  res <- c()
  xq = seq(L_q,H_q,by=0.1)
  for (i in 1 : length(means)) {
    yq <- pnorm(xq, mean = means[i], sd = sds[i])
    res <- c(res, uniToOthers(xq, yq, runif(1)))
  }
  #plot(xq,yq)
  return(res)
}