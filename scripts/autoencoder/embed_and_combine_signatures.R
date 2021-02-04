#!/usr/bin/Rscript  

args <- commandArgs(TRUE)

if (length(args) < 7 ) {
stop("\n\t\033[31mThe following inputs are required\033[0m:\n[1] Job Name
[2] Embed [a] 'training' or [b] 'testing' 
[3] gene-model job name 
[4] domain-model job name
[5] family-model job name
[6] superfamily-model job name
[7] fold-model job name
\tIf embeding testing data the following options are requred 
\t(formatted the same as needed by generate_autoencoder_models.R): 
[8] testing genes data
[9] testing domain data
[10] testing family data
[11] testing superfamily data 
[12] testing fold data
\033[34mGenerated output: \033[35m'Job Name'.csv\033[0m", call.=FALSE)
} 
  
# args[1] = "test"
# args[2] = "training"
# args[3] = "gtex_rda/gtex.genes.model"
# args[4] = "gtex_rda/gtex.domain.model"
# args[5] = "gtex_rda/gtex.family.model"
# args[6] = "gtex_rda/gtex.superfamily.model"
# args[7] = "gtex_rda/gtex.fold.model"

# args[8]  = "go-analysis/allcombined.GO.genelist.csv.format"
# args[9] = "go-analysis/allcombined.GO.domain.csv.format"
# args[10] = "go-analysis/allcombined.GO.fam.csv.format" 
# args[11] = "go-analysis/allcombined.GO.sfam.csv.format"
# args[12] = "go-analysis/allcombined.GO.fold.csv.format"

jobid = as.character(args[1] ) 
type  = as.character(args[2] )  #Either a gene list, or the outputs from structural signatures 
if (type == "testing" && length(args) < 12 ) 
{
stop("
\t\033[31mIf embeding testing data the following arguments are requred 
\t(formatted the same as needed by generate_autoencoder_models.R):\033[0m 
[8] testing genes data
[9] testing domain data
[10] testing family data
[11] testing superfamily data 
[12] testing fold data" )
}
library(data.table)
library(tidyverse)
library(keras)
library(BBmisc)  
print("loading models")



load(paste0(args[3], ".rda")) ## gene
geneinfo = out_data
genemod <- load_model_hdf5(paste0(args[3],".h5"))

load(paste0(args[4], ".rda")) ## domain
dominfo = out_data
dommod <- load_model_hdf5(paste0(args[4],".h5"))

load(paste0(args[5], ".rda")) ## family 
faminfo = out_data
fammod <- load_model_hdf5(paste0(args[5],".h5"))

load(paste0(args[6], ".rda")) ## superfamily 
sfaminfo = out_data
sfammod <- load_model_hdf5(paste0(args[6],".h5"))

load(paste0(args[7], ".rda")) ## fold 
foldinfo = out_data
foldmod <- load_model_hdf5(paste0(args[7],".h5"))


generate_embedings = function(model, data , meta  )
{
    embed <- predict(model, data.matrix(data)) %>% as.data.frame
    embed = cbind(embed, meta) %>% as.data.frame
    embed
}

combined_embedings = function(geneembed , domainembed, famembed, sfamembed , foldembed )
{
    gene.fold = merge(geneembed, foldembed , by = "ID" )
    gene.fold$class.x = NULL
    gene.fold$class.y = NULL
    gene.fold.domain = merge(gene.fold, domainembed, by = "ID" )
    gene.fold.domain$class = NULL
    gene.fold.domain.sfam =   merge(gene.fold.domain, sfamembed , by = "ID" )
    gene.fold.domain.sfam$class = NULL
    gene.fold.domain.sfam.fam = merge(gene.fold.domain.sfam, famembed, by = "ID" )
    gene.fold.domain.sfam.fam 
}

format_data_to_training = function(testing.data , training.info , type )
{
    reference_table = 
        cbind(
            rep("setupvariable", 
                length(training.info$header[3:length(training.info$header)])) , 
            rep("setupvariable", 
                length(training.info$header[3:length(training.info$header)])) , 
            training.info$header[3:length(training.info$header)], 
            rep(1, 
                length(training.info$header[3:length(training.info$header)]))
            )  %>% as.data.table()
    if (type == "gene")
    {
        header_genes = c(
           "ID",
            "class" , 
            "gene", 
            "rank"
        )
        names(testing.data) = header_genes
    } else if ( type == "ss")
    {
        header_ss = c(
            "structure", 
            "counts_observed",
            "background_counts", 
            "number_of_genes_in_set",
            "total_number_proteins_proteome",
            "pvalue",
            "fdr",
            "bonforroni_cutoff",
            "log_fold_change",
            "ID",
            "class"
        )
        names(testing.data) = header_ss
    }
    testing.data$presense = rep(1, nrow(testing.data))  
    dattype = training.info$params[[1]]
    testing.data$datacol = testing.data[, ..dattype]
    
    if ( type == "ss")
    {   
        testing.data = testing.data[,c("ID","class","structure", "datacol")]
        testing.data = testing.data[which(testing.data$structure %in% training.info$header[3:length(training.info$header)]) , ] 
    }
    else if ( type == "gene")
    {
        testing.data = testing.data[which(testing.data$gene %in% training.info$header[3:length(training.info$header)]) , ] 
    }
    if ( dattype == "pvalue" )
    {
        testing.data$datacol = -log10(testing.data$datacol ) 
    }
    if ( type == "ss" )
    {
        testing.data = testing.data[,c("ID","class","structure", "datacol")]
        names(reference_table) = c("ID", "class", "structure", "datacol")
        testing.data = rbind(reference_table, testing.data)
        testing.data.wide  = spread(testing.data, key = structure, fill = 0 , value = datacol )
    }   else if (type == "gene")
    {
        testing.data = testing.data[,c("ID","class","gene", "datacol")]
        names(reference_table) = c("ID", "class", "gene", "datacol")
        testing.data = rbind(reference_table, testing.data)
        testing.data.wide  = spread(testing.data, key = gene, fill = 0 , value = datacol )
    }
    testing.data.wide = testing.data.wide[! which(testing.data.wide$ID == "setupvariable"), ]
    row.names(testing.data.wide) = testing.data.wide$ID
    testing.data.y = testing.data.wide[,c("ID", "class")]
    testing.data.x = testing.data.wide[,3:ncol(testing.data.wide)] %>% data.matrix()  %>% normalize( method = "range", range = c(0, 1), margin = 1) 
    output = list(testing.y = testing.data.y , testing.x = testing.data.x )
    output 
}


if (type == "training") 
{
    print("generating training embedding")
    
    geneembed = generate_embedings(genemod, geneinfo$training_data, geneinfo$training_meta)
    names(geneembed) =c(paste0("gene", 1:as.numeric(geneinfo$params[2])), "ID", "class")

    domainembed = generate_embedings(dommod, dominfo$training_data, dominfo$training_meta)
    names(domainembed) =c(paste0("domain", 1:as.numeric(dominfo$params[2])), "ID", "class")
    
    famembed = generate_embedings(fammod, faminfo$training_data, faminfo$training_meta)
    names(famembed) =c(paste0("family", 1:as.numeric(faminfo$params[2])), "ID", "class")
    
    sfamembed  = generate_embedings(sfammod, sfaminfo$training_data, sfaminfo$training_meta)
    names(sfamembed ) =c(paste0("superfamily", 1:as.numeric(sfaminfo$params[2])), "ID", "class")
    
    foldembed = generate_embedings(foldmod, foldinfo$training_data, foldinfo$training_meta)
    names(foldembed) =c(paste0("fold", 1:as.numeric(foldinfo$params[2])), "ID", "class")
    
    combined = combined_embedings(geneembed, domainembed, famembed, sfamembed  ,  foldembed)
    
    write.table(combined, paste0(jobid,".csv"), col.names = T, eol = "\n", quote = F, sep = "," ,
     row.names = F )
} else if (type == "testing") 
{
    
    print("formating testing data")
    train.gene = fread(args[8], header = F, sep = "," , stringsAsFactors= F )
    train.domain = fread(args[9], header = F, sep = "," , stringsAsFactors= F )
    train.fam = fread(args[10], header = F, sep = "," , stringsAsFactors= F )
    train.sfam = fread(args[11], header = F, sep = "," , stringsAsFactors= F )
    train.fold =  fread(args[12], header = F, sep = "," , stringsAsFactors= F )
    print(paste0("working on ", args[8]))
    train.gene.format = format_data_to_training(  train.gene, geneinfo, "gene"  )
    print(paste0("working on ", args[9]))
    train.domain.format = format_data_to_training(  train.domain, dominfo, "ss"  )
    print(paste0("working on ", args[10]))
    train.fam.format = format_data_to_training(  train.fam  , faminfo, "ss"  )
    print(paste0("working on ", args[11]))
    train.sfam.format = format_data_to_training(  train.sfam, sfaminfo, "ss"  )
    print(paste0("working on ", args[12]))
    train.fold.format = format_data_to_training(  train.fold, foldinfo, "ss"  )
    print("generating testing embedding")
    
    geneembed = generate_embedings(genemod, train.gene.format$testing.x, train.gene.format$testing.y)
    names(geneembed) =c(paste0("gene", 1:as.numeric(geneinfo$params[2])), "ID", "class")
    
    domainembed = generate_embedings(dommod, train.domain.format$testing.x, train.domain.format$testing.y)
    names(domainembed) =c(paste0("domain", 1:as.numeric(dominfo$params[2])), "ID", "class")
    
    famembed = generate_embedings(fammod, train.fam.format$testing.x, train.fam.format$testing.y)
    names(famembed) =c(paste0("family", 1:as.numeric(faminfo$params[2])), "ID", "class")

    sfamembed= generate_embedings(sfammod, train.sfam.format$testing.x, train.sfam.format$testing.y)
    names(sfamembed) =c(paste0("superfamily", 1:as.numeric(sfaminfo$params[2])), "ID", "class")
    
    foldembed = generate_embedings(foldmod, train.fold.format$testing.x, train.fold.format$testing.y)
    names(foldembed) =c(paste0("fold", 1:as.numeric(foldinfo$params[2])), "ID", "class")
    
    combined = combined_embedings( geneembed, domainembed, famembed, sfamembed , foldembed)
    write.table(combined, paste0(jobid,".csv"), col.names = T, eol = "\n", quote = F, sep = "," ,
     row.names = F )
}