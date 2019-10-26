data = read.delim("artificialDatasetTemplate.txt", header = TRUE, sep = "\t")

############################## load libraries ##############################
library(foreach)
library(doParallel)

############################## load sources ##############################
source("lib/r1Reads.R")
source("lib/r2Reads.R")

############################## parameters ##############################
cores = 4L
m = 150L


registerDoParallel(cores)
#length(i <- grep("[AGTC]+$", data[['Sequence']]))

############################## srtuctures and data init ##############################
reads1 = data.frame(Sequence = character(0), Type = character(0), Length = integer(0), Idx = integer(0), Q = character(0), stringsAsFactors = FALSE)
reads2 = data.frame(Sequence = character(0), Type = character(0), Length = integer(0), Idx = integer(0), Q = character(0), stringsAsFactors = FALSE)

x = seq(33,126,by=0.01)
#round(min(x))
y <- pnorm(x, mean = 79.5, sd = 15)
#plot(x,y)

############################## main start ##############################
for (i in 1:length(data[,1])){
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
      }else {
        reps = coverage - reps_sum
      }
      reads1 <- rbind(reads1, cbind(r1Reads(reps = reps, mystring = toString(data[k, 'Sequence']), m = m, x = x, y = y), Idx = i, Type = toString(data[k, 'Type'])))
      reads2 <- rbind(reads2, cbind(r2Reads(reps = reps, mystring = toString(data[k, 'Sequence']), m = m, x = x, y = y), Idx = i, Type = toString(data[k, 'Type'])))
    }
  }
}