# VCF Files Comparator And Summarization Tool

Compare the Vcf files with foundings mutations and the amplicons' template from which vcf was created.

Furthermore, the second script creates a summary table of all mutations (TP, FP, TN, FN).


### Prerequisites

R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"

### Installing

Just run the VcfComp_install.R and all required packages will be installed automaticaly.

## Executing the Main Script


#### VCF Files Comparator

Just execute the wraper bash script as it is mentioned bellow:

./VcfComp.sh -v **$pref/SampleName_GATK_4_1_0_0_variants_MRPAS$i.vcf** -t **artificialDatasetTemplate.txt** -c **Coordinates.txt** -o **$pref/SampleName.results.table.MRPAS$i** -l **MRPAS$i**

#### Summarization Tool

Just execute the wraper bash script as it is mentioned bellow:

./VcfComp-summary.sh -p **$pref/** -d **result** -c **CHROM,POS.Template,REF.Template,ALT.Template,DP.Template,AD.Template,AF.Template,POS.Data,REF.data,ALT.Data,DiffDP,DiffAD,DiffAF,Position.diff**

**-c** flag is a string with comma separated values each of them indicates the column of the vcf file that we want to include in the summary table.

## Versioning

First version of Vcf Comparator Tools

## Authors

* **Anastasios Mitsigkolas** 
* **Fotis Psomopoulos** 


## License

This project is licensed under the MIT License.