# This code will install required packages if they are not already installed

if (!require("foreach")) {
  install.packages("foreach")
  #library(foreach)
}
if (!require("doParallel")) {
  install.packages("doParallel")
  #library(doParallel)
}
if (!require("optparse")) {
  install.packages("optparse")
  #library(optparse)
}
if (!require("stringr")) {
  install.packages("stringr")
  #library(stringr)
}
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("ChIPsim")) {
  BiocManager::install("ChIPsim")
  #library(ChIPsim)
}
if (!require("Biostrings")) {
  BiocManager::install("Biostrings")
  #library(Biostrings)
}
