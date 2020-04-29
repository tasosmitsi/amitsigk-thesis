############################## load libraries ##############################
library(foreach, quietly = T, warn.conflicts = F)
library(doParallel, quietly = T, warn.conflicts = F)
library(ChIPsim, quietly = T, warn.conflicts = F)
library(stringr, quietly = T, warn.conflicts = F)
library(Biostrings, quietly = T, warn.conflicts = F)
library(optparse, quietly = T, warn.conflicts = F)
############################## load sources & functions init ###############
# source("lib/overlapedReads/leftToRight_read.R")
# source("lib/overlapedReads/rightToLeft_read.R")
# source("lib/overlapedReads/Reads.R")
# source("lib/dif_means_gaussian_random_samples.R")
# source("lib/uniToOthers.R")
# source("lib/comb.R")

leftToRight_read <- function(L, m, mystring) {
  n = nchar(mystring)
  if (L >= (-m+2) & L <= (n)) {
    if (L < 1) {
      ns = abs(L) + 1
      read = paste(strrep("N", ns), substr(mystring, 1, m - ns), strrep("N", if ((m - ns - n) > 0) (m - ns - n)  else 0), sep="")
    }else if (L >= 1) {
      ns = m - (n - L + 1)
      read = paste(substr(mystring, L, n + ns), strrep("N", if (ns < 0) 0 else ns), sep="")
    }
    return(read)
  }else {
    return(FALSE)
  }
}

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

comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

############################## set arguments ##############################
option.list <- list(make_option("--template.file.path", type = "character", default = "artificialDatasetTemplate.txt", help = "Full path to template tab-delimited file.", metavar = "character"),
                    make_option("--cores", type = "integer", default = 1L, help = "Set number of cores to use when creating reads and read quals. Preferable. The default is 1 core.", metavar = "integer"),
                    make_option("--reads.length", type = "integer", default = 150L, help = "Set the length of reads and read quals. Preferable. The default is 150 bases.", metavar = "integer"),
                    make_option("--meanCov.mult", type = "integer", default = 1L, help = "Set the multiplier of template meanCoverage. Preferable. The default is 1 .", metavar = "integer"),
                    make_option("--results.folder.path", type = "character", default = "results", help = "Path of the results folder (ex. '/results').", metavar = "character"),
                    make_option("--compress", type = "character", default = "yes", help = "Compress output files? 'yes'/'no' (default: yes)", metavar = "character"))


opt.parser <- OptionParser(option_list = option.list)
opt <- parse_args(opt.parser)

# Check file.paths is given.
if (is.null(opt$template.file.path)){
  print_help(opt.parser)
  stop("The absolute path to template tab-delimited file needs to be provided.", call. = F)
}
if (is.null(opt$cores) | opt$cores <= 0){
  print_help(opt.parser)
  stop("The number of cores must be a positive integer.", call. = F)
}
if (is.null(opt$reads.length) | opt$reads.length <= 0){
  print_help(opt.parser)
  stop("The length of reads must be a positive integer.", call. = F)
}
if (is.null(opt$meanCov.mult) | opt$meanCov.mult <= 0){
  print_help(opt.parser)
  stop("The meanCov multiplier must be a positive integer.", call. = F)
}
if (is.null(opt$compress) | !(opt$compress == "no" | opt$compress == "yes")){
  print_help(opt.parser)
  stop("The compress parameter must be either yes or no.", call. = F)
}
if (is.null(opt$results.folder.path)){
  print_help(opt.parser)
  stop("The absolute path to results folder needs to be provided.", call. = F)
}

############################## Init variables and paths ##############################
template.file.path <- opt$template.file.path
cores <- opt$cores
m <- opt$reads.length
results.folder.path <- paste0("/", opt$results.folder.path)
mult <- opt$meanCov.mult
if (opt$compress == "yes") {
  compress_out <- TRUE
} else if (opt$compress == "no"){
  compress_out <- FALSE
}

# Check results folder exists, or create one.
if (! dir.exists(paste0(getwd(), results.folder.path))) {
  dir.create(path = paste0(getwd(), results.folder.path))
  message(paste0("*** Created [script_path]", results.folder.path, " directory."))
} else {
  message(paste0("*** Directory [script_path]", results.folder.path, " located."))
}

results.path <- paste0(getwd(), results.folder.path, "/")

name.R1 <- if (compress_out) paste0(results.path, "r1.compressed.fastq.gz") else paste0(results.path, "r1.fastq")
name.R2 <- if (compress_out) paste0(results.path, "r2.compressed.fastq.gz") else paste0(results.path, "r2.fastq")

############################## parameters ##############################

registerDoParallel(cores)
#length(i <- grep("[AGTC]+$", data[['Sequence']]))
############################## srtuctures and data init ##############################

# realistic Quality, quote if needs
# x = seq(1,1000,by=1)
# 
# means <- dnorm(x, 500, 500)
# means <- means[450:599]
# x = seq(1,m,by=1)
# 
# sds = seq(0.8 , 5, by = 0.0279)
# sds = sds[1:150]
# kati = 700 / max(means) 
# means <- (means * kati)-630
####################################

# best Quality, quote if needs
means <- rep(63,m)
sds <- rep(0.1,m)
max(means)
##############################
L_q = 33
H_q = 75


data = read.delim(template.file.path, header = TRUE, sep = "\t")
data$MeanCoverage <- data$MeanCoverage * mult
#data <- data[1:4, ]

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
message("Data pre-calculations ends in ", Sys.time() - start_time, " secs", sep=" ")

start_time <- Sys.time()
for (j in 1 : ceiling(length(data[,1]) / cores)) {
  res <- foreach(i = (((j - 1) * cores) + 1) : (j * cores), .combine='comb') %dopar% {
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
  
  seqs <- c(res[["R1"]][["Sequence"]])
  names(seqs) <- c(res[["R1"]][["AmpliconID"]])
  seqs <- BStringSet(seqs)
  quals <- c(res[["R1"]][["Q"]])
  quals <- BStringSet(quals)
  writeXStringSet(seqs, name.R1, format = "fastq", qualities = quals, compress = compress_out, append = TRUE)
  
  seqs <- c(res[["R2"]][["Sequence"]])
  names(seqs) <- c(res[["R2"]][["AmpliconID"]])
  seqs <- BStringSet(seqs)
  quals <- c(res[["R2"]][["Q"]])
  quals <- BStringSet(quals)
  writeXStringSet(seqs, name.R2, format = "fastq", qualities = quals, compress = compress_out, append = TRUE)
  
  message(paste0("Have been calculated ", (sum(data$MeanCoverage[1 : (j * cores)]) / sum(data$MeanCoverage[1 : length(data[,1])])) * 100), "% of the overall data.")
  
  rm(z)
}

stopImplicitCluster()
message("Random readings created in", Sys.time() - start_time, "h", sep = " ")