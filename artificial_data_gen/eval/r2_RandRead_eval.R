r2_RandRead_eval <- function (data) {
  testData <- data
  testData[['Sequence']] <- str_replace_all(testData[['Sequence']], "[N]", "")
  testData[['Length']] <- nchar(testData[['Sequence']])
  reads <- c()
  for (i in 1 : length(testData[,1])) {
    reads <- c(reads, if (LCS(toString(data[['Sequence']][i]), testData[['Sequence']][i]) == testData[['Sequence']][i]) TRUE else FALSE)
  }
  if (length(which(reads == FALSE)) != 0) {
    message("Something went wrong...with r2")
    FALSE
  }else {
    message("Everything is ok with r2!")
    TRUE
  }
}