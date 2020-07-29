# This code will install required packages if they are not already installed

if (!require("optparse")) {
  install.packages("optparse")
  #library(optparse)
}
if (!require("stringr")) {
  install.packages("stringr")
  #library(stringr)
}
if (!require("stringdist")) {
  install.packages("stringdist")
  #library(stringdist)
}
if (!require("marray")) {
  install.packages("marray")
  #library(marray)
}