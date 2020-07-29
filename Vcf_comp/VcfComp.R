############################## load libraries ##############################
library(stringr, quietly = T, warn.conflicts = F)
library(stringdist, quietly = T, warn.conflicts = F)
library(optparse, quietly = T, warn.conflicts = F)
library(marray, quietly = T, warn.conflicts = F)

############################## load sources ##############################
#source("lib/.R")
#source("lib/findPosIndexes.R")
#source("lib/isInsideBound.R")
#source("lib/findIndelsType.R")
#source("lib/changeIndelsFormat.R")
#source("lib/findRelations.R")
# source("lib/startDif.R")
# source("lib/perPosOfSpeChrom.R")
# source("lib/try.setwd.R")

findPosIndexes <- function(vcf_pos, startPos, endPos) {
  return(which(startPos <= vcf_pos & vcf_pos <= endPos))
}

isInsideBound <- function(indexes) {
  if (length(indexes) == 0) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}

findIndelsType <- function(ref, alt) {
  if (nchar(toString(ref)) == nchar(toString(alt))) {
    return ("change")
  }else if (nchar(toString(ref)) > nchar(toString(alt))) {
    return ("deletion")
  }else if (nchar(toString(ref)) < nchar(toString(alt))) {
    return ("insertion")
  }
}

changeIndelsFormat <- function(ref, alt, type) {
  if (type == "insertion") {
    ref <- "_"
    alt <- str_extract(alt, "\\B[AGTC]*")
  }else if (type == "deletion") {
    ref <- str_extract(ref, "\\B[AGTC]*")
    alt <- "_"
  }
  return(data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE))
  #return(t(c(ref,alt)))
}

findRelations <- function(templateRef, templateAlt, vfcRef, vfcAlt, posIndexes) {
  res <- data.frame(templateId = integer(), vfcId = integer())
  if (length(posIndexes) == length(vfcRef)) {
    for (i in 1 : length(vfcRef)) {
      for (posIndex in posIndexes[[i]]) {
        if (templateRef[posIndex] == vfcRef[i] & templateAlt[posIndex] == vfcAlt[i]) {
          #do something to store matching
          res <- rbind(res, data.frame(templateId = posIndex, vfcId = i))
          break
        }
      }
    }
  }
  return(res)
}

startDif <- function(string1, string2) {
  n = pmin(nchar(string1), nchar(string2))
  i = 1
  while (i <= n) {
    if (substr(string1, i, i) != substr(string2, i, i)) 
      return(i)
    i = i + 1
  }
}

perPosOfSpeChrom <- function(data) {
  res <- data.frame()
  for(chrom in levels(data$CHROM)) {
    res <- rbind(res, data.frame(PerPosOfSpeChrom = length(which(final_res$CHROM == chrom & final_res$Positive)) / length(which(final_res$CHROM == chrom))))
  }
  row.names(res) <- levels(data$CHROM)
  return(res)
}

try.setwd <- function(x) {
  out <- tryCatch( {
    message("*** Setting working directory in [script_path]/results.")
    setwd(x)
  },
  error = function(cond){
    message("Cannot set working directory correctly. Please check directory rights.")
    message(cond)
    return(NA)
  },
  warning = function(cond){
    message("Warning issued:")
    message(cond)
    return(NULL)
  }, 
  finally = {message("*** Working directory is:")}
  )
  
  return(getwd())
}

############################## set arguments ##############################
option.list <- list(make_option("--vcf.file.path", type = "character", default = "SampleName.results.table", help = "Full path to vcf results table file.", metavar = "character"),
                    make_option("--template.file.path", type = "character", default = "artificialDatasetTemplate.txt", help = "Full path to template tab-delimited file.", metavar = "character"),
                    make_option("--coordinates.file.path", type = "character", default = "Coordinates.txt", help = "Full path to coordinates.txt file.", metavar = "character"),
                    make_option("--results.folder.name", type = "character", default = "default", help = "Name of the results folder.", metavar = "character"))

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

keep.date <- date()
keep.date.paste <- paste(unlist(strsplit(str_replace_all(keep.date, ":", "-"), " ")), collapse = "#")

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

if (nrow(gatk_data) == 0) {
  stop("Now mutation was found... The input table is empty!", call. = F)
}

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

endPos <- vector()
endPos[as.numeric(rownames(templateRefs[sanityCheck$TemplateId,]))] <- coordinates$End.Position[sanityCheck$coordinatesId]
maxLength <- max(c(length(template[,1]), length(endPos)))
length(endPos) <- maxLength
template <- cbind(template, endPos)



mutPos <- vector()
for (i in 1 : length(template[,1])) {
  if (template$Type[i] == 'Ref') {
    refIdx <- i
  }else{
    template$MeanCoverage[i] = template$MeanCoverage[refIdx]
    mutPos[i] <- template$startPos[refIdx] + startDif(toString(template$Sequence[refIdx]), toString(template$Sequence[i])) - 1
    template$startPos[i] <- template$startPos[refIdx]
    template$endPos[i] <- template$endPos[refIdx]
  }
}
maxLength <- max(c(length(template[,1]), length(mutPos)))
length(mutPos) <- maxLength
template <- cbind(template, mutPos)

rm(i, startPos, mutPos, endPos)

template <- template[-which(template$Type == 'Ref'), ]
#template <- template[,!(names(template) %in% c("MeanCoverage"))]

#chrs <- str_match_all(toString(template[ , 'AmpliconID']), "chr(\\w+)\\.")[[1]][ ,2]


############################################
# positions from AmpliconId
# positions <- as.data.frame(str_match_all(toString(template[ , 'AmpliconID']), "chr(\\w+)\\.(\\d+)\\.(\\d+)")[[1]][ ,2:4],
#                            stringsAsFactors = FALSE)
#colnames(positions) = c("chr", "StartPos", "EndPos")


###########################################
# positions from coordinates fille
positions <- data.frame(StartPos = template$startPos, EndPos = template$endPos)

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

AD_splited <- strsplit(as.character(gatk_data$SampleName.AD), ",", fixed=TRUE)
AD_splited <- lapply(data.frame(matrix(unlist(AD_splited), nrow=length(AD_splited), byrow=T)), function(x) { if(is.factor(x)) as.numeric(as.character(x)) else x })

if (nrow(relations) > 0) {
  temp <- vector()
  temp[relations$vfcId] <- relations$templateId
  temp1 <- vector()
  temp1[relations$vfcId] <- as.numeric(row.names(template[relations$templateId,]))
  
  temp2 <- vector()
  temp2[as.numeric(rownames(gatk_data[-relations$vfcId,]))] <- FALSE
  temp2[as.numeric(rownames(gatk_data[relations$vfcId,]))] <- TRUE
  
  results <- data.frame(CHROM = gatk_data$CHROM,
                        
                        IsInsideBoundaries.Data = unlist(isInsideBoundaries),
                        Positive = temp2,
                        POS.Data = gatk_data$POS,
                        POS.Template = template$mutPos[temp],
                        Position.diff = template$mutPos[temp] - gatk_data$POS,
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
} else {
  results <- data.frame(CHROM = gatk_data$CHROM,
                        
                        IsInsideBoundaries.Data = unlist(isInsideBoundaries),
                        Positive = FALSE,
                        POS.Data = gatk_data$POS,
                        POS.Template = NA,
                        Position.diff = NA,
                        IndelType.Data = indels_type,
                        
                        REF.data = indels$REF,
                        REF.Template = NA,
                        ALT.Data = indels$ALT,
                        ALT.Template = NA,
                        AF.Data = gatk_data$SampleName.AF,
                        AF.Template = NA,
                        AD.Data = gatk_data$SampleName.AD,
                        AD.Template = NA,
                        DP.Data = gatk_data$DP,
                        DP.Template = NA,
                        MBQ.Data = gatk_data$MBQ,
                        MMQ.Data = gatk_data$MMQ,
                        
                        DiffDP = NA,
                        DiffAD = NA,
                        DiffAF = NA,
                        
                        PerDiffInfo.DP = NA,
                        PerDiffFormat.DP = NA,
                        PerDiffAD = NA,
                        
                        AbsoluteMatchingTemplateId = NA,
                        AbsoluteMatchingOfficialTemplateId = NA,
                        
                        stringsAsFactors = FALSE)
  hasNoRelations <- template
}

hasNoRelations <- hasNoRelations[,!(names(hasNoRelations) %in% c("Type"))]
hasNoRelations <- cbind(hasNoRelations, data.frame(Negative = logical(length = length(hasNoRelations[,1])), AD = hasNoRelations$MeanCoverage * hasNoRelations$VAF))
colnames(hasNoRelations) = c("CHROM", "Sequence", "DP.Template", "AF.Template", "REF.Template", "ALT.Template", "StartPos.Template", "EndPos.Template", "POS.Template", "Negative", "AD.Template")

hasNoRelations$CHROM <- str_match_all(toString(hasNoRelations[ , 'CHROM']), "chr(\\w+)")[[1]][,2]
final_res <- merge(results, hasNoRelations, all = TRUE)

############################## summary ##############################

if (nrow(relations) > 0) {
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
  write.list(summary, filename = paste0(results.path, "summary_", keep.date.paste, "_", results.folder.name), append = FALSE, closefile = TRUE)
}

write.table(final_res, paste0(results.path, "results.table_", keep.date.paste, "_", results.folder.name), append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

##################################test#########################################################################
if (!is.null(final_res$Positive)) {
  if (length(which(final_res$Positive == TRUE)) > 0) {
    of_template = read.delim(template.file.path, header = TRUE, sep = "\t")
    test <- final_res[which(final_res$Positive == TRUE), c('AbsoluteMatchingOfficialTemplateId', 'Position.diff', "IndelType.Data")]
    test <- cbind(test, data.frame(AmpliconID = of_template$AmpliconID[test$AbsoluteMatchingOfficialTemplateId]))
    
    for (i in 1:length(test[,1])) {
      for (j in 1:length(test[,1])) {
        if (all(toString(test$AmpliconID[i]) == toString(test$AmpliconID[j])) & i != j) {
          test$Position.diff[i] <- paste(test$Position.diff[i], test$Position.diff[j], sep = ";")
          test$IndelType.Data[i] <- paste(test$IndelType.Data[i], test$IndelType.Data[j], sep = ";")
          test <- test[-j,] 
        }
      }
    }
    write.table(test[,c('AmpliconID', 'Position.diff', 'IndelType.Data')], paste0(results.path, "pos_difs_per_amplicon_", keep.date.paste, "_", results.folder.name), append = FALSE, sep = "\t", dec = ".",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}
##############################################################################################################