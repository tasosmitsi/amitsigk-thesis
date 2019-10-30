############################## load libraries ##############################
library(stringr)
library(PTXQC)

############################## load sources ##############################
source("lib/r1Reads.R")
source("lib/uniToOthers.R")
source("eval/r1_RandRead_eval.R")

############################## parameters ##############################
string <- "dsgvbksdbsbvhjsbdfvjhbdfjvhbdjfhvbav"
m <- 10
reps <- 1000
x = seq(33,126,by=0.01)
y <- pnorm(x, mean = 79.5, sd = 15)

result <- r1_RandRead_eval(r1Reads(reps = reps, m = m, mystring = string, x = x, y = y))