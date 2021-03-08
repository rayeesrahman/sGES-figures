## Code for generating figures in "Protein structure-based gene expression signatures" 
Rahman, R et al. 2021, PNAS

### Instructions 
To generate figures from the article you can run the scripts found in each directory during an interactive R session, where each line of the script is sent to the R terminal. 

Each script for each figure is standalone and can be run individually. However all scripts require the data directory in the figure directory to be populated. 

There are a few utility scripts included in this repository to parse GTEx data or to create and train autoencoders in R. 

This code has been validated to run on the following platforms: Ubuntu 20.04 LTS, Ubuntu 18.04 LTS and Windows Subsystem for linux. Your mileage may vary for other platforms. 

### Data directory 

**All of the code references a `figures/data/` directory** This is a directory that houses data used to generate figures. 

Due to the size of the generated data (>4GB), it is not included in this repo, but is available to download from here: http://iyengarlab.org/dtoxs/PNAS-sGES.tar.bz2 

Please extract the files and place folder into `figures/data/`. The directory structure should be as follows: 
``` 
data/
├── TISSUES-harmonizome
│   ├── inte-gl
│   ├── results
│   └── rna-gl
├── archs
│   ├── ARCHS4_r
│   ├── embeddings
│   └── structural-signatures
│       └── formatted
├── dtox
├── gtex
│   ├── stable_signature_data
│   │   ├── selected_gene_lists
│   │   └── stable_signature_data_old
│   └── structural-signatures
├── humanprotatlas
│   ├── genelists
│   ├── raw
│   └── results
├── l1000
├── msigdb
│   ├── genelists
│   ├── raw-gmt
│   └── results
│       └── ss_output
├── reconstruction_errors
└── rocs
```

### Directories

`figures/` Contains R code and scripts used to generate figures for the paper. 

`figures/figure-1/` code for figures 1C (supplemental figure S2) 

`figures/figure-3/` code for figures 3A-D (supplemental figure S4)  

`figures/figure-4/` code for figure 4 (supplemental figure S6) 

`figures/figure-5/` code for figures 5A-D (supplemental figures S8-11) 

`figures/figure-6/` code for figures 6A-E (supplemental figures S12-17) 

`figures/figure-7/` code for figures 7A-C 

`scripts/autoencoder` scripts to train autoencoder models and output reconstruction errors for test data from GES or sGES data

`utilities/` simple shell scripts to run strucutural signatures on GTEX and ARCHS4 data  

`figures/data` generated data that is analyzed by the scripts 


### Dependancies 

#### R libraries: 

tools

BBmisc

RColorBrewer

Rtsne

data.table

dplyr

ggdendro

ggplot2

gplots

grid

gridExtra

iheatmapr

keras

magrittr

pROC

plotly

plyr

ranger

rhdf5

stringr

tidyr

tidyverse

vegan

#### Python Libraries (will be installed during R keras installation )

Tensorflow 

Keras 

#### Dependancy installation 

Use the following code to install the R dependancies. Python tensorflow will also be installed during the installation of the repsective R versions 

`library("tools","BBmisc","RColorBrewer","Rtsne","data.table","dplyr","ggdendro","ggplot2","gplots","grid","gridExtra","keras","magrittr","pROC","plotly","plyr","ranger",,"stringr","tidyr","tidyverse","vegan")`

`if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")` 

`BiocManager::install("rhdf5", "iheatmapr")`  

#### Contact 

For any issues relating to this repository contact the corresponding authors: rayees(dot)rahman(at)icahn(dot)mssm(dot)edu or avner(dot)schlessinger(at)mssm(dot)edu 
