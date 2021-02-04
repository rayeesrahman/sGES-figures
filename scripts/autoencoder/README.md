# Embedding Signatures using an Stacked Denoising Autoencoder

## Scripts 

Running the `generate_autoencoder_signature_model.R` script will generate 

1) an hdf5 file with the embedding model 
2) an Rdata file containining training data, headers, performance, etc  

The `embed_and_combined_signatures.R` script uses the outputs of  `generate_autoencoder_signature_model.R` to generate embeddings of either the training data or new validation data. It will create a single csv output with the embeddings. 

`train_and_evaluate_model.R` will take training data and testing data embeded by `embed_and_combined_signatures.r` and train a neural network to predict classes depending on the training 'class' column. 

`generate_roc_curve.R` will take the roc output from train_and_evaluate_model and plot ROC curves based on the input data

## Getting help

Running the scripts without any arguments will output a help with all the required arguments needed to run the code 

## Directories 
The `*_rda/` directories contain autoencoder models and rdata files for the databases used. 

The `embed_models/` directory contains embbedings 

