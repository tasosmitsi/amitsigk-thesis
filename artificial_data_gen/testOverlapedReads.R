############################## load libraries ##############################
library(foreach)
library(doParallel)
library(ChIPsim)
library(stringr)

############################## load sources ##############################
source("lib/overlapedReads/leftToRight_read.R")
source("lib/overlapedReads/rightToLeft_read.R")
source("lib/overlapedReads/Reads.R")

source("lib/dif_means_gaussian_random_samples.R")
source("lib/uniToOthers.R")
source("lib/comb.R")

############################## parameters ##############################
cores = 4L

registerDoParallel(cores)
#length(i <- grep("[AGTC]+$", data[['Sequence']]))

############################## srtuctures and data init ##############################
m = 150L
x = seq(1,1000,by=1)

means <- dnorm(x, 500, 500)
means <- means[450:599]
x = seq(1,m,by=1)

sds = seq(0.8 , 5, by = 0.0279)
sds = sds[1:150]
kati = 700 / max(means) 
means <- (means * kati)-630
max(means)
plot(x,means)

L_q = 33
H_q = 75


data = read.delim("artificialDatasetTemplate.txt", header = TRUE, sep = "\t")
data <- data[1:3, ]

############################## main start ##############################
start_time <- Sys.time()
message("Starting data pre-calculations...")

data[ , 'AmpliconID'] <- str_remove(data[ , 'AmpliconID'], "tile_\\d")

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

res <- foreach(i = 1:length(data[,1]), .combine='comb') %dopar% {
  if (data[i, 'MeanCoverage'] > 0) {
    indexes = list(data.frame(Idx = i, Type = toString(data[i, 'Type'])), data.frame(Idx = i, Type = toString(data[i, 'Type'])))
    lst <- Reads(AmpliconID = data[i, 'AmpliconID'], reps = data[i, 'MeanCoverage'], mystring = toString(data[i, 'Sequence']), m = m, i = i, means = means, sds = sds, L_q, H_q)
    mapply(cbind, lst, indexes, SIMPLIFY=FALSE)
  }else{
    list(data.frame(), data.frame())
  }
}
z <- sapply(res[["R1"]], is.factor)
res[["R1"]][z] <- lapply(res[["R1"]][z], as.character)

z <- sapply(res[["R2"]], is.factor)
res[["R2"]][z] <- lapply(res[["R2"]][z], as.character)

rm(z)
stopImplicitCluster()

writeFASTQ(res[["R1"]][["Sequence"]], res[["R1"]][["Q"]], res[["R1"]][["AmpliconID"]], file="E:/WSL-shared/test/r1.fastq")
writeFASTQ(res[["R2"]][["Sequence"]], res[["R2"]][["Q"]], res[["R2"]][["AmpliconID"]], file="E:/WSL-shared/test/r2.fastq")


message("Data random readings ends in ", Sys.time() - start_time, " mins", sep="")
