#!/bin/bash
usage() {                                      	# Function: Print a help message.
  echo "Usage: cmd [-v Absolute path to the vcf file] [-t Absolute path to the template file] [-c Absolute path to the coordinates file] [-o Output file name] [-l Results folder name]" 1>&2 
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}
while getopts ":v:t:c:o:l:" opt; do
  case ${opt} in
    v ) vcf=${OPTARG}
		;;
	t ) template=${OPTARG}
		;;
	c ) coordinates=${OPTARG}
		;;
	o ) outName=${OPTARG}
		;;
	l ) resName=${OPTARG}
		;;
    : ) exit_abnormal
		;;
	* ) exit_abnormal
		;;
  esac
done

gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F DP -F ECNT -F MBQ -F MFRL -F MMQ -F MPOS -F POPAF -F TLOD -GF GT -GF AD -GF AF -GF DP -GF F1R2 -GF F2R1 -GF SB -O $outName 

Rscript  VcfComp.R --vcf.file.path=$outName --template.file.path=$template --coordinates.file.path=$coordinates --results.folder.name=$resName