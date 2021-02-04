# Code for generating figures in "Protein structure-based gene expression signatures" 
## Rahman, R et al. 2021, PNAS

### Directories

`figures/` Contains R code and scripts used to generate figures for the paper. These scripts can run using an interactive R session 

`figures/figure-1/` code for figures 1C (supplemental figure S2) 

`figures/figure-3/` code for figures 3A-D (supplemental figure S4)  

`figures/figure-4/` code for figure 4 (supplemental figure S6) 

`figures/figure-5/` code for figures 5A-D (supplemental figures S8-11) 

`figures/figure-6/` code for figures 6A-E (supplemental figures S12-17) 

`figures/figure-7/` code for figures 7A-C 

`figures/sfigures/` gcode for supplemental figures 1,5,7,18

**All of the code references a `data/` directory** This is a directory that houses data used to generate figures. Due to the size of the generated data (>4GB), it is not included in this repo, but is available upon request. 

`scripts/autoencoder` scripts to train autoencoder models and output reconstruction errors for test data from GES or sGES data

`utilities/` simple shell scripts to run strucutural signatures on GTEX and ARCHS4 data  


### R libraries: 

data.table

magrittr

tidyr

vegan

ggplot2

ranger

pROC

plyr

dplyr

tidyverse

rhdf5

tensorflor

keras

tools

You can use the following command to install dependancies in R:

`install.packages(c("tidyverse", "data.table" ,"magrittr", "vegan", "ggplot2", "ranger", "pROC", "plyr", "dplyr", "rhdf5", "tools" ))`


#### Python Libraries 

Tensorflow 
Keras 