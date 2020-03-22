#!/bin/bash
usage() {                                      	# Function: Print a help message.
  echo "Usage: cmd [-p Absolute path to the project folder] [-d first delimiter of name of result folders] [-c Columns names]" 1>&2 
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}
while getopts ":p:d:c:" opt; do
  case ${opt} in
    p ) path=${OPTARG}
		;;
	d ) delimiter=${OPTARG}
		;;
	c ) columnNames=${OPTARG}
		;;
    : ) exit_abnormal
		;;
	* ) exit_abnormal
		;;
  esac
done


Rscript  summary_of_iterations.R --project.folder.path=$path --result.folders.first.delim=$delimiter --columns=$columnNames