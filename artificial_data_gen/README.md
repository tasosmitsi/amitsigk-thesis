# Artificial Data Generator

It is a set of algorithms written in R to generate synthetic next-generation sequencing reads by using an amplicons' template. 

### Prerequisites

R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"

### Installing

Just run the adg_overlap_install.R and all required packages will be installed automaticaly.

## Executing the Main Script

Execute the following script:

Rscript  adg_overlap.R --cores=**$cores** --meanCov.mult=**1** --reads.length=**150** --template.file.path=**artificialDatasetTemplate.txt** --compress=**yes** --results.folder.path=**$pref**

## Versioning

First version of Artificial Data Generator

## Authors

* **Anastasios Mitsigkolas** 
* **Fotis Psomopoulos** 


## License

This project is licensed under the MIT License.