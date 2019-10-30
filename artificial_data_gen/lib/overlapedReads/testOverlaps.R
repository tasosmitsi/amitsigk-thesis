library(stringi)

source("lib/overlapedReads/rightToLeft_read.R")
source("lib/overlapedReads/leftToRight_read.R")

string = "tasosmitsigkolas"
m = 10
n = nchar(string)

rand_start = sample(1:n, 1)
rand_offset = sample(0:(m-1),1)

Hrtl = rand_start - rand_offset + m - 1
Lltr = rand_start + rand_offset - m + 1

leftToRight_read(L = Lltr, m = m, mystring = string)
rightToLeft_read(H = Hrtl, m = m, mystring = string)
