############################## load libraries ##############################
library(stringr)
library(PTXQC)

############################## load sources ##############################

############################## parameters ##############################
test <- res


test[['R1']][['Sequence']] <- str_replace_all(test[['R1']][['Sequence']], "[N]", "")
test[['R1']][['Length']] <- nchar(test[['R1']][['Sequence']])


test[['R2']][['Sequence']] <- str_replace_all(test[['R2']][['Sequence']], "[N]", "")
test[['R2']][['Length']] <- nchar(test[['R2']][['Sequence']])

r1 <- c()
r2 <- c()
for (i in 1 : length(test[['R1']][,1])) {
  r1 <- c(r1, if (LCS(toString(res[['R1']][['Sequence']][i]), test[['R1']][['Sequence']][i]) == test[['R1']][['Sequence']][i]) TRUE else FALSE)
  r2 <- c(r2, if (LCS(toString(res[['R2']][['Sequence']][i]), test[['R2']][['Sequence']][i]) == test[['R2']][['Sequence']][i]) TRUE else FALSE)
}

if (length(which(r1 == FALSE)) != 0) {
  message("Something went wrong...with r1")
}else {
  message("Everything is ok with r1!")
}

if (length(which(r2 == FALSE)) != 0) {
  message("Something went wrong...with r2")
}else {
  message("Everything is ok with r2!")
}