############################## load libraries ##############################
library(stringr)
library(stringdist)
library(optparse, quietly = T, warn.conflicts = F)
library(marray)

############################## load sources ##############################
#source("lib/.R")
source("lib/findPosIndexes.R")
source("lib/isInsideBound.R")
source("lib/findIndelsType.R")
source("lib/changeIndelsFormat.R")
source("lib/findRelations.R")
source("lib/startDif.R")
source("lib/perPosOfSpeChrom.R")
source("lib/try.setwd.R")


############################## set arguments ##############################
option.list <- list(make_option("--vcf.file.path", type = "character", default = "SampleName.results.table", help = "Full path to vcf results table file.", metavar = "character"),
                    make_option("--template.file.path", type = "character", default = "artificialDatasetTemplate.txt", help = "Full path to template tab-delimited file.", metavar = "character"),
                    make_option("--coordinates.file.path", type = "character", default = "Coordinates.txt", help = "Full path to coordinates.txt file.", metavar = "character"),
                    make_option("--results.folder.name", type = "character", default = "default", help = "Name of the results folder.", metavar = "character"))
                    # make_option("--min.read.len", type = "integer", default = 75, help = "Minimum length of produced reads (default: 75).", metavar = "integer"), 
                    # make_option("--max.read.len", type = "integer", default = 200, help = "Maximum length of produced reads (default: 200). max.read.len > amplicon lengths is supported by extending tile sequences.", metavar = "integer"), 
                    # make_option("--max.insert.len", type = "integer", default = 350, help = "Maximum length of produced inserts (default: 350). A max.insert.len  > amplicon lengths is supported by extending tile sequences.", metavar = "integer"), 
                    # make_option("--cores", type = "integer", default = 1L, help = "Set number of cores to use when creating reads and read quals. Preferable. Requires linux and R package 'parallel'.", metavar = "integer"),
                    # make_option("--mate1.suffix", type = "character", default = "1:N:0:1", help = "Set custom suffix (default: '1:N:0:1')", metavar = "character"), 
                    # make_option("--mate2.suffix", type = "character", default = "2:N:0:1", help = "Set custom suffix (default: '2:N:0:1')", metavar = "character"), 
                    # make_option("--compress", type = "character", default = "no", help = "Compress output files? 'yes'/'no' (default: no)", metavar = "character"))

opt.parser <- OptionParser(option_list = option.list)
opt <- parse_args(opt.parser)

# Check file.paths is given.
if (is.null(opt$vcf.file.path)){
  print_help(opt.parser)
  stop("The absolute path to the vcf results table needs to be provided.", call. = F)
}
if (is.null(opt$template.file.path)){
  print_help(opt.parser)
  stop("The absolute path to template tab-delimited file needs to be provided.", call. = F)
}
if (is.null(opt$coordinates.file.path)){
  print_help(opt.parser)
  stop("The absolute path to coordinates.txt file needs to be provided.", call. = F)
}

############################## Init variables and paths ##############################
vcf.file.path <- opt$vcf.file.path
template.file.path <- opt$template.file.path
coordinates.file.path <- opt$coordinates.file.path
results.folder.name <- opt$results.folder.name

# Check results folder exists, or create one.
if (! dir.exists(paste0(getwd(), "/results_", results.folder.name))) {
  dir.create(path = paste0(getwd(), "/results_", results.folder.name))
  message(paste0("*** Created [script_path]/results_", results.folder.name, " directory."))
} else {
  message(paste0("*** Directory [script_path]/results_", results.folder.name, " located."))
}

results.path <- paste0(getwd(), "/results_", results.folder.name, "/")

############################## copy files ##############################
#file.copy(Vcf.file.path, current_path, overwrite = TRUE)

############################## read files ##############################
template = read.delim(template.file.path, header = TRUE, sep = "\t")
gatk_data = read.delim(vcf.file.path, header = TRUE, sep = "\t")
coordinates = read.delim(coordinates.file.path, header = TRUE, sep = "\t")

templateRefs <- template[-which(template$Type != 'Ref'), ]

sanityCheck <- data.frame()
for (i in 1 : length(templateRefs[,1])){
  for (j in 1 : length(coordinates[,1])){
    if (is.null(startDif(toString(templateRefs[i,3]), toString(coordinates[j,10])))) {
      sanityCheck <- rbind(sanityCheck, data.frame(TemplateId = i, coordinatesId = j))
      break
    }
  }
}
#differences <- data.frame()
#for (i in 1 : length(coordinates[,1])) {
#  differences <- rbind(differences, data.frame(differenceOfPos = (coordinates$End.Position[i]-coordinates$Start.Position[i]-nchar(toString(coordinates[i,10])))))
#}
#coordinates <- cbind(coordinates, differences)

############# - Template post proccessing - #############
startPos <- vector()
startPos[as.numeric(rownames(templateRefs[sanityCheck$TemplateId,]))] <- coordinates$Start.Position[sanityCheck$coordinatesId]
maxLength <- max(c(length(template[,1]), length(startPos)))
length(startPos) <- maxLength
template <- cbind(template, startPos)

mutPos <- vector()
for (i in 1 : length(template[,1])) {
  if (template$Type[i] == 'Ref') {
    refIdx <- i
  }else{
    template$MeanCoverage[i] = template$MeanCoverage[refIdx]
    mutPos[i] <- template$startPos[refIdx] + startDif(toString(template$Sequence[refIdx]), toString(template$Sequence[i])) - 1
    template$startPos[i] <- template$startPos[refIdx]
  }
}
template <- cbind(template, mutPos)

rm(i)
rm(mutPos)

template <- template[-which(template$Type == 'Ref'), ]
#template <- template[,!(names(template) %in% c("MeanCoverage"))]

#chrs <- str_match_all(toString(template[ , 'AmpliconID']), "chr(\\w+)\\.")[[1]][ ,2]
positions <- as.data.frame(str_match_all(toString(template[ , 'AmpliconID']), "chr(\\w+)\\.(\\d+)\\.(\\d+)")[[1]][ ,2:4],
                           stringsAsFactors = FALSE)
colnames(positions) = c("chr", "StartPos", "EndPos")


posIndexes <- lapply(gatk_data$POS, findPosIndexes, startPos = positions$StartPos, endPos = positions$EndPos)
isInsideBoundaries <- lapply(posIndexes, isInsideBound)
indels_type <- apply(gatk_data[,c('REF','ALT')], 1, function(y) findIndelsType(ref = y['REF'], alt = y['ALT']))

temp <- apply(cbind(gatk_data[,c('REF','ALT')], indels_type), 1,
                function(y) changeIndelsFormat(ref = y['REF'], alt = y['ALT'], type = y['indels_type']))

indels <- data.frame(matrix(unlist(temp), nrow=length(temp), byrow=T), stringsAsFactors=FALSE)
colnames(indels) = c("REF", "ALT")
rm(temp)

relations <- findRelations(templateRef = template[,'REF'], templateAlt = template[,'ALT'],
                           vfcRef = indels[,'REF'], vfcAlt = indels[,'ALT'], posIndexes = posIndexes)

temp <- vector()
temp[relations$vfcId] <- relations$templateId
temp1 <- vector()
temp1[relations$vfcId] <- as.numeric(row.names(template[relations$templateId,]))

temp2 <- vector()
temp2[as.numeric(rownames(gatk_data[-relations$vfcId,]))] <- FALSE
temp2[as.numeric(rownames(gatk_data[relations$vfcId,]))] <- TRUE

AD_splited <- strsplit(as.character(gatk_data$SampleName.AD), ",", fixed=TRUE)
AD_splited <- lapply(data.frame(matrix(unlist(AD_splited), nrow=length(AD_splited), byrow=T)), function(x) { if(is.factor(x)) as.numeric(as.character(x)) else x })

results <- data.frame(CHROM = gatk_data$CHROM,
                      
                      IsInsideBoundaries.Data = unlist(isInsideBoundaries),
                      Positive = temp2,
                      POS.Data = gatk_data$POS,
                      POS.Template = template$mutPos[temp],
                      IndelType.Data = indels_type,
                      
                      REF.data = indels$REF,
                      REF.Template = template$REF[temp],
                      ALT.Data = indels$ALT,
                      ALT.Template = template$ALT[temp],
                      AF.Data = gatk_data$SampleName.AF,
                      AF.Template = template$VAF[temp],
                      AD.Data = gatk_data$SampleName.AD,
                      AD.Template = template$MeanCoverage[temp] * template$VAF[temp],
                      DP.Data = gatk_data$DP,
                      DP.Template = template$MeanCoverage[temp],
                      MBQ.Data = gatk_data$MBQ,
                      MMQ.Data = gatk_data$MMQ,
                      
                      DiffDP = template$MeanCoverage[temp] - gatk_data$DP,
                      DiffAD = (template$MeanCoverage[temp] * template$VAF[temp]) - AD_splited$X2,
                      DiffAF = template$VAF[temp] - gatk_data$SampleName.AF,
                      
                      PerDiffInfo.DP = 100 * ((template$MeanCoverage[temp] * template$VAF[temp]) - gatk_data$DP) / (template$MeanCoverage[temp] * template$VAF[temp]),
                      PerDiffFormat.DP = 100 * ((template$MeanCoverage[temp] * template$VAF[temp]) - gatk_data$SampleName.DP) / (template$MeanCoverage[temp] * template$VAF[temp]),
                      PerDiffAD = 100 * ((template$MeanCoverage[temp] * template$VAF[temp]) - (AD_splited$X1 + AD_splited$X2)) / (template$MeanCoverage[temp] * template$VAF[temp]), 
                      
                      AbsoluteMatchingTemplateId = temp,
                      AbsoluteMatchingOfficialTemplateId = temp1,
                      
                      stringsAsFactors = FALSE
                      )
rm(temp)
rm(temp1)
rm(temp2)

hasNoRelations <- template[-relations$templateId,]
hasNoRelations <- hasNoRelations[,!(names(hasNoRelations) %in% c("Type"))]
hasNoRelations <- cbind(hasNoRelations, data.frame(Negative = FALSE, AD = hasNoRelations$MeanCoverage * hasNoRelations$VAF))
colnames(hasNoRelations) = c("CHROM", "Sequence", "DP.Template", "AF.Template", "REF.Template", "ALT.Template", "StartPos.Template", "POS.Template", "Negative", "AD.Template")

hasNoRelations$CHROM <- positions[-relations$templateId, 'chr']

final_res <- merge(results, hasNoRelations, all = TRUE)

############################## summary ##############################


summary <- list(
  "Template number of mutexes" = length(template[,1]),
  "Template overal depth" = sum(template$MeanCoverage),
  
  "GATK number of mutexes" = length(gatk_data[,1]),
  "GATK overal INFO DP" = sum(final_res$DP.Data[which(!is.na(final_res$DP.Data))]),
  "GATK overal FORMAT DP" = sum(gatk_data$SampleName.DP),
  "GATK overal AD" = sum(AD_splited$X1 + AD_splited$X2),
  "Uninformative reads (INFO.DP - AD)" = sum(final_res$DP.Data[which(!is.na(final_res$DP.Data))]) - sum(AD_splited$X1 + AD_splited$X2),
  
  "TruePositives" = length(which(final_res$Positive)),
  "FalsePositives" = length(which(!final_res$Positive)),
  "FalseNegatives" = length(which(!final_res$Negative)),
  "Data that is inside boundaries" = length(which(final_res$IsInsideBoundaries.Data)),
  "Data that isn't inside boundaries" = length(which(!final_res$IsInsideBoundaries.Data)),
  
  "GATK % of insertions" = length(which(final_res$IndelType.Data == "insertion")) / length(gatk_data[,1]),
  "GATK % of deletions" = length(which(final_res$IndelType.Data == "deletion")) / length(gatk_data[,1]),
  "GATK % of changes" = length(which(final_res$IndelType.Data == "change")) / length(gatk_data[,1]),
  
  "% TruePositives" = length(which(final_res$Positive)) / length(gatk_data[,1]),
  "% FalsePositives" = length(which(!final_res$Positive)) / length(gatk_data[,1]),
  "% FalseNegatives" = length(which(!final_res$Negative)) / length(gatk_data[,1]),
  
  "% of positives of specific CHROM" = perPosOfSpeChrom(data.frame(CHROM = final_res$CHROM, Positive = final_res$Positive)),
  
  "% of Positives that is inside boundaries" = length(which(final_res$IsInsideBoundaries.Data)) / length(which(!is.na(final_res$Positive))),
  
  "Mean of absolutes values of difference of AF" = mean(abs(final_res$DiffAF[which(!is.na(final_res$DiffAF))])),
  "SD of absolutes values of difference of AF" = sd(abs(final_res$DiffAF[which(!is.na(final_res$DiffAF))])),
  
  "Mean of absolutes values of percentage DP - INFO.DP" = mean(abs(final_res$PerDiffInfo.DP[which(!is.na(final_res$PerDiffInfo.DP))])),
  "SD of absolutes values of percentage DP - INFO.DP" = sd(abs(final_res$PerDiffInfo.DP[which(!is.na(final_res$PerDiffInfo.DP))])),
  
  "Mean of absolutes values of percentage DP - FORMAT.DP" = mean(abs(final_res$PerDiffFormat.DP[which(!is.na(final_res$PerDiffFormat.DP))])),
  "SD of absolutes values of percentage DP - FORMAT.DP" = sd(abs(final_res$PerDiffFormat.DP[which(!is.na(final_res$PerDiffFormat.DP))])),
  
  "Mean of absolutes values of percentage DP - AD" = mean(abs(final_res$PerDiffAD[which(!is.na(final_res$PerDiffAD))])),
  "SD of absolutes values of percentage DP - AD" = sd(abs(final_res$PerDiffAD[which(!is.na(final_res$PerDiffAD))]))
)

keep.date <- date()
keep.date.paste <- paste(unlist(strsplit(date(), " ")), collapse = "_")

write.table(final_res, paste0(results.path, "results.table_", "keep.date.paste", "_", results.folder.name), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

write.list(summary, filename = paste0(results.path, "summary_", "keep.date.paste", "_", results.folder.name), append = FALSE, closefile = TRUE)