############################## load libraries ##############################
library(optparse, quietly = T, warn.conflicts = F)
library(stringr)

############################## load sources ##############################
#source("lib/.R")

############################## set arguments ##############################
option.list <- list(make_option("--project.folder.path", type = "character", default = getwd(), help = "Full path to project folder.", metavar = "character"),
                    make_option("--result.folders.first.delim", type = "character", default = "result", help = "Type the first delimiter of result folders name.", metavar = "character"),
                    make_option("--columns", type = "character", default = "REF.data,ALT.Data,AF.Template,DiffAF,DP.Template,DiffDP,AD.Template,DiffAD,REF.Template,ALT.Template", help = "Type the column names separated with comma.", metavar = "character"))



opt.parser <- OptionParser(option_list = option.list)
opt <- parse_args(opt.parser)

# Check file.paths is given.
if (is.null(opt$project.folder.path)){
  print_help(opt.parser)
  stop("The absolute path to the project's folder needs to be provided.", call. = F)
}
if (is.null(opt$result.folders.first.delim)){
  print_help(opt.parser)
  stop("The first delimiter for the result folders needs to be provided.", call. = F)
}
if (is.null(opt$columns)){
  print_help(opt.parser)
  stop("The column names separated with comma needs to be provided.", call. = F)
}


############################## Init variables and paths ##############################
project.folder.path <- opt$project.folder.path
result.folders.first.delim <- opt$result.folders.first.delim

column_names <- unlist(strsplit(opt$columns, ",", fixed=TRUE))
template_column_names <- column_names[grep("\\.Template", column_names)]
data_column_names <- column_names[-grep("\\.Template", column_names)]

template_column_names <- unlist(lapply(template_column_names, function(x) paste0("\\.", x, "$")))
data_column_names <- unlist(lapply(data_column_names, function(x) paste0("\\.", x, "$")))

############################## copy files ##############################


############################## read files ##############################
paths <- list.files(pattern = "results.table_", recursive = TRUE, full.name = TRUE)

i <- 1 
dflist <- c()
for (path in paths) {
  tmp <- str_extract(basename(path), "MRPAS\\d+")
  dflist <- c(dflist, tmp)
  assign(tmp, read.delim(path, header = TRUE, sep = "\t"))
  i <- i + 1
}
rm(path)
dflist <- sort(dflist, decreasing = FALSE)
############################## main ##############################

overall_summary = data.frame(matrix(NA, nrow = length(get(dflist[1])[[1]]), ncol = 1))
for (i in dflist) {
  tmp <- get(i)
  colnames(tmp) <- paste(i, colnames(tmp), sep = ".")
  overall_summary <- cbind(overall_summary, tmp)
}
overall_summary <- overall_summary[-1]

template_result <- data.frame(tmp[unlist(lapply(template_column_names, grep, x = colnames(tmp)))])
data_result <- data.frame(overall_summary[unlist(lapply(data_column_names, grep, x = colnames(overall_summary)))])
data_result <- data_result[unlist(lapply(paste0(dflist, "\\."), grep, x = colnames(data_result)))] #sorting of the data_result data.frame
#data_result <- na.omit(data_result, na.action = "omit")
#rownames(data_result) <- NULL
result <- cbind(template_result, data_result)

result[] <- lapply(result, function(x) {  #apply abs and round to all non numeric cells
  if (is.numeric(x)) 
    abs(round(x, 3))
  else
    x
  })

write.table(result, "overall_summary.table", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))