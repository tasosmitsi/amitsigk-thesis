############################## load libraries ##############################
library(stringr)
#library(vcfR)

############################## load sources ##############################
#source("lib/.R")
source("lib/findPosIndexes.R")
source("lib/isFalsePos.R")
source("lib/findIndelsType.R")
source("lib/changeIndelsFormat.R")
source("lib/findRelations.R")

############################## Init variables and paths ##############################
current_path <- getwd()
my_art_data_path <- "E:/WSL-shared/MyArtData"

colNames <- c("posIndexes", "isFalsePositive")
results <- data.frame(matrix(ncol = length(colNames), nrow = 0), stringsAsFactors = FALSE)
colnames(results) <- colNames
rm(colNames)

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


posIndexes <- lapply(gatk_data$POS, findPosIndexes, startPos = positions$StartPos, endPos = positions$EndPos)
falsePositives <- lapply(indexes, isFalsePos)
indels_type <- apply(gatk_data[,c('REF','ALT')], 1, function(y) findIndelsType(ref = y['REF'], alt = y['ALT']))

test <- apply(cbind(gatk_data[,c('REF','ALT')], indels_type), 1,
                function(y) changeIndelsFormat(ref = y['REF'], alt = y['ALT'], type = y['indels_type']))
indels <- data.frame(matrix(unlist(test), nrow=length(test), byrow=T), stringsAsFactors=FALSE)
colnames(indels) = c("REF", "ALT")
rm(test)


relations <- findRelations(templateRef = template[,'REF'], templateAlt = template[,'ALT'],
                           vfcRef = indels[,'REF'], vfcAlt = indels[,'ALT'], posIndexes = posIndexes)



