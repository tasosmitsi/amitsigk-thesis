############################## load libraries ##############################
library(stringr)

############################## load sources ##############################
source("lib/.R")

############################## Init variables and paths ##############################
current_path <- getwd()
my_art_data_path <- "E:/WSL-shared/MyArtData"

############################## copy files ##############################
file.copy(paste(my_art_data_path, "/SampleName_GATK_4_1_0_0_variants.results.table", sep = ""), current_path, overwrite = TRUE)
file.copy(paste(my_art_data_path, "/SampleName_freebayes-v1_3_1.results.table", sep = ""), current_path, overwrite = TRUE)

############################## read files ##############################
template = read.delim("../artificial_data_gen/artificialDatasetTemplate.txt", header = TRUE, sep = "\t")
gatk_data = read.delim("SampleName_GATK_4_1_0_0_variants.results.table", header = TRUE, sep = "\t")
freebayes_data = read.delim("SampleName_freebayes-v1_3_1.results.table", header = TRUE, sep = "\t")

############# - Template post proccessing - #############
template <- template[-which(template$Type == 'Ref'), ]
template <- template[,!(names(template) %in% c("MeanCoverage"))]

#chrs <- str_match_all(toString(template[ , 'AmpliconID']), "chr(\\w+)\\.")[[1]][ ,2]
positions <- as.data.frame(str_match_all(toString(template[ , 'AmpliconID']), "chr(\\w+)\\.(\\d+)\\.(\\d+)")[[1]][ ,2:4],
                           stringsAsFactors = FALSE)
colnames(positions) = c("chr", "StartPos", "EndPos")