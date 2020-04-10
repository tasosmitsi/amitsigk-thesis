Reads <- function(AmpliconID, reps, m, i, mystring, means, sds, L_q, H_q) {
  resr1 = data.frame()
  resr2 = data.frame()
  if (reps > 0) {
    n = nchar(mystring)
    rand_start = sample(1:n, reps, replace = TRUE)
    rand_offset = sample(0:(m-1), reps, replace = TRUE)
    
    for (j in 1 : reps) {
      Hrtl = rand_start[j] - rand_offset[j] + m - 1
      Lltr = rand_start[j] + rand_offset[j] - m + 1
      
      read <- leftToRight_read(L = Lltr, m = m, mystring = mystring)
      quality <- intToUtf8(round(dif_means_gaussian_random_samples(L_q = L_q, H_q = H_q, means = means, sds = sds)))
      resr1 <- rbind(resr1, data.frame(
        AmpliconID = paste(AmpliconID, i, j, " 1:N:0:1", sep=""), 
        Sequence = read, Length = nchar(read), Q = quality))
      
      
      read <-rightToLeft_read(H = Hrtl, m = m, mystring = mystring)
      quality <- intToUtf8(round(dif_means_gaussian_random_samples(L_q = L_q, H_q = H_q, means = means, sds = sds)))
      resr2 <- rbind(resr2, data.frame(
        AmpliconID = paste(AmpliconID, i, j, " 2:N:0:1", sep=""),
        Sequence = read, Length = nchar(read), Q = quality))
    }
    return(list(R1 = resr1, R2 = resr2))
  }else {
    stop("reps argument must be a possitive integer.") 
  }
}