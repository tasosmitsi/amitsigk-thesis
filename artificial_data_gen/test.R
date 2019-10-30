data = read.delim("artificialDatasetTemplate.txt", header = TRUE, sep = "\t")
data <- data[1:5, ]
############################## load libraries ##############################
library(foreach)
library(doParallel)
#library(stringi)

############################## load sources ##############################
source("lib/r1Reads.R")
source("lib/r2Reads.R")
source("lib/uniToOthers.R")
source("lib/comb.R")

############################## parameters ##############################
cores = 4L
m = 150L


registerDoParallel(cores)
#length(i <- grep("[AGTC]+$", data[['Sequence']]))

############################## srtuctures and data init ##############################
#reads1 = data.frame(Sequence = character(0), Type = character(0), Length = integer(0), Idx = integer(0), Q =character(0), stringsAsFactors = FALSE)
#reads2 = data.frame(Sequence = character(0), Type = character(0), Length = integer(0), Idx = integer(0), Q =character(0), stringsAsFactors = FALSE)

x = seq(33,126,by=0.01)
#round(min(x))
y <- pnorm(x, mean = 79.5, sd = 15)
#plot(x,y)

############################## main start ##############################
start_time <- Sys.time()
message("Starting data pre-calculations...")

for (i in 1:length(data[,1])) {
  if (toString(data[i, 'Type']) == "Ref") {
    coverage = data[i, 'MeanCoverage']
    start_idx = i
    end_idx = i + 1
    repeat {
      if (end_idx > length(data[,1])) {
        end_idx = end_idx - 1
        break
      }else if (toString(data[end_idx, 'Type']) == "Ref") {
        end_idx = end_idx - 1
        break
      }else {
        end_idx = end_idx + 1
      }
    }
    reps_sum = 0L
    for (k in end_idx : start_idx) {
      if (k != start_idx){
        reps = floor(data[k, 'VAF'] * coverage)
        reps_sum = reps_sum + reps
        data[k, 'MeanCoverage'] <- reps
      }else {
        reps = coverage - reps_sum
        data[k, 'MeanCoverage'] <- reps
      }
      
    }
  }
}
message("Data pre-calculations ends in ", Sys.time() - start_time, " mins", sep="")

start_time <- Sys.time()
message("Starting data random readings...")
writeLines(c(""), "log.txt")
res <- foreach(i = 1:length(data[,1]), .combine='comb', .packages="stringi") %dopar% {
  r1 <- data.frame()
  r2 <- data.frame()
  if (data[i, 'MeanCoverage'] > 0) {
    tic <- Sys.time()
    r1 <- cbind(r1Reads(reps = data[i, 'MeanCoverage'], mystring = toString(data[i, 'Sequence']), m = m, x = x, y = y), Idx = i, Type = toString(data[i, 'Type']))
    tocr1 <- difftime(Sys.time(), tic, units = c("mins")) 
    tic <- Sys.time()
    r2 <- cbind(r2Reads(reps = data[i, 'MeanCoverage'], mystring = toString(data[i, 'Sequence']), m = m, x = x, y = y), Idx = i, Type = toString(data[i, 'Type']))
    tocr2 <- difftime(Sys.time(), tic, units = c("mins"))
    sink("log.txt", append=TRUE)
    cat(paste("Iteration ", i," r1 takes: ", tocr1, " mins, r2 takes: ", tocr2, "mins.", "\n"))
  }
  list(R1 = r1, R2 = r2)
}
stopImplicitCluster()
message("Data random readings ends in ", Sys.time() - start_time, " mins", sep="")
