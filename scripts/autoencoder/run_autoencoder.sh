#!/bin/bash

##arg 1  jobname 
##arg 2  feature ##logfc counts_observed presense pvalue 
##arg 3 autoencoder size 

#./generate_autoencoder_signature_model.R $1.genes  data_dir/allcombined.250.genes.gtex  gene $3 presense 50 

#./generate_autoencoder_signature_model.R $1.fold  data_dir/allcombined.250.fold.csv.3.gtex  fold $3 $2 50 

#./generate_autoencoder_signature_model.R $1.superfamily  data_dir/allcombined.250.superfamily.csv.3.gtex  superfamily $3 $2 50 

#./generate_autoencoder_signature_model.R $1.family  data_dir/allcombined.250.family.csv.3.gtex  family $3 $2 50 

#./generate_autoencoder_signature_model.R $1.domain  data_dir/allcombined.250.domain.csv.3.gtex  domain $3 $2 50 

#./embed_and_combine_signatures.R $1.gtex training $1.genes $1.domain $1.family $1.superfamily $1.fold

./embed_and_combine_signatures.R archs-no-outliers.archs testing archs-with-outliers.genes archs-with-outliers.domain archs-with-outliers.family archs-with-outliers.superfamily archs-with-outliers.fold  data_dir/archs.tissue.genelist.formated.no.outlier  data_dir/allcombined.archs.250.domain.csv.formatted.no.outlier data_dir/allcombined.archs.250.family.csv.formatted.no.outlier data_dir/allcombined.archs.250.superfam.csv.formatted.no.outlier data_dir/allcombined.archs.250.fold.csv.formatted.no.outlier

#sed -E 's/.GTEX$//' $1.gtex.csv >  $1.gtex.csv.2 

#mv $1.gtex.csv.2 $1.gtex.csv 

./train_and_evaluate_model.R archs-no-outliers.model archs-with-outliers.gtex.csv archs-no-outliers.archs.csv no yes yes yes