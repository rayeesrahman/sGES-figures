library(data.table)
library(tidyverse)
library(keras)
library(BBmisc)  
library(data.table)
rm(list = ls())

##Figure 6A/S12, GES autoencoder reconstruction errors
gen_encoder_arch = function(inshape, bottle,  lossfx = "mean_squared_error")
{
    input_layer <- 
    layer_input(shape = inshape)  
    encoder <- 
        input_layer %>% 
        layer_dense(units = 100, activation = "relu") %>%
        layer_dense(units = 50, activation = "relu") %>%
        layer_dense(units = 25, activation = "relu") %>%
        layer_dense(units = bottle,  name = "bottle") # 3 dimensions for the output layer
    decoder <- 
        encoder %>% 
        layer_dense(units = 25, activation = "relu") %>% 
        layer_dense(units = 50, activation = "relu") %>% 
        layer_dense(units = 100, activation = "sigmoid") %>% 
        layer_dense(units = inshape)
    autoencoder_model <- keras_model(inputs = input_layer, outputs = decoder)
    autoencoder_model %>% compile(
        loss = lossfx, 
        optimizer = "adam"
    )
    return(autoencoder_model)
}  
 
embed_ss_signatures = function(data , type = "pvalue" , datatype  , bottle = 3,  epochs = 25 )
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
    names(data) = header_ss
    ref_header = data.table
    if ( datatype == "domain")
    {
        ref_header = fread("../files/entry.list", sep = "\t" ,  header = F, stringsAsFactors= F )
        ref_header = ref_header[V2 == "Domain", V1 ] 
    }  else if ( datatype == "fold")
    {
        ref_header = fread("../files/scope_total_2.06-stable.txt", sep = "|", header = F, stringsAsFactors = F, quote="")
        ref_header = ref_header[V1 == "folds", V2 ]
    }  else if ( datatype == "family")
    {
        ref_header = fread("../files/scope_total_2.06-stable.txt", sep = "|", header = F, stringsAsFactors = F, quote="")
        ref_header = ref_header[V1 == "families", V2 ]
    }   else if ( datatype == "superfamily")
    {
        ref_header = fread("../files/scope_total_2.06-stable.txt", sep = "|", header = F, stringsAsFactors = F, quote="")
        ref_header = ref_header[V1 == "superfamilies", V2 ]
    }
    if ( type == "presense")
    {
        data$presense = rep(1, nrow(data))  
    }
    datacol = data[, ..type] 
    data = data[,c("ID","class","structure")] 
    ####pvalues need to rescaled by their -log since the fill is 0 
    if ( type == "pvalue" )
    {
        data$datacol = -log10(datacol ) 
    }
    reference_table = cbind(rep("setupvariable", length(ref_header)) , 
        rep("setupvariable", length(ref_header)) , 
        ref_header , 
        rep(1, length(ref_header)))  %>% as.data.table()
    names(reference_table) = c("ID", "class", "structure", "datacol" )      
    data = data[which(data$structure %in% ref_header) , ] 
    data = rbind(reference_table, data) %>% as.data.table 
    print("Converting data from long to wide")
    data.wide  = spread(data, key = structure, fill = 0 , value = datacol )
    data.wide = data.wide[! which(data.wide$ID == "setupvariable"), ]
    print("Training autoencoder model")
    row.names(data.wide) = data.wide$id 
    data.y = data.wide[,c("ID", "class")]
    data.x = data.wide[,3:ncol(data.wide)] %>% data.matrix() # %>% normalize( method = "range", range = c(0, 1), margin = 1) 
    ##to make a stacked denoising autoencoder 
    data.x.noise = apply(data.x, 2, function(x) { noise = rnorm(length(x) , .01 , .005) ;  x + noise} ) 
    autoencoder_model = gen_encoder_arch( inshape = ncol(data.x), bottle = bottle ) 
    history = autoencoder_model %>% 
        fit(
            x = data.matrix(data.x.noise) , 
            y = data.matrix(data.x) , 
            epochs = epochs, 
            batch_size =  32,
            validation_split = 0.2,
            view_metrics = TRUE,
        )
    header = names(data.wide)
    intermediate_layer_model <- keras_model(inputs = autoencoder_model$input, outputs = get_layer(autoencoder_model, "bottle")$output)
    return_data = list(embed_model =intermediate_layer_model , header = header , training_data = data.x , training_meta = data.y , params = c(type = type, bottleneck = bottle, epochs = epochs) , performance = history )
    return(return_data)
}


gtex.gene250 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.250.csv", header = F , stringsAsFactors = F)
archs.gene250 = fread("figures/data/archs/archs.tissue.genelist.formated", header = F , stringsAsFactors = F)
gtex.gene250$V2 =NULL 
gtex.n = gtex.gene250$V3 %>% unique 
archs.n = archs.gene250$V2 %>% unique 
shared = gtex.n[which(gtex.n %in% archs.n)]
gtex.gene250 = gtex.gene250[gtex.gene250$V3 %in% shared, ]
archs.gene250 = archs.gene250[archs.gene250$V2 %in% shared, ]
ids  = gtex.gene250$V1 %>% unique 
valsize = length(ids) *.2
valsize = ceiling(valsize)
gtex.val.set = sample(ids, valsize )
gtex.gene250$V5 = rep(1, nrow(gtex.gene250)) 
gtex.gene250$type = rep("gtex", nrow(gtex.gene250)) 
archs.gene250$type = rep("archs", nrow(archs.gene250))
header_genes = c(
        "id",
        "class" , 
        "gene", 
        "datacol", 
        "type"
    )
names(gtex.gene250) <- header_genes
names(archs.gene250) <- header_genes
combined = rbind(gtex.gene250, archs.gene250)

combined.wide  = spread(combined , key = gene, fill = 0 , value = datacol )
row.names(combined.wide ) = combined.wide$id 
gtex.gene250.wide = combined.wide[combined.wide$type == "gtex", ]
archs.gene250.wide = combined.wide[combined.wide$type == "archs", ]

gtex.validation =  gtex.gene250.wide[ gtex.gene250.wide$id %in%  gtex.val.set , ]
gtex.training =  gtex.gene250.wide [ ! gtex.gene250.wide $id %in%  gtex.val.set , ]

gtex.validation.y = gtex.validation[,c("id", "class","type")]
gtex.validation.x = gtex.validation[,4:ncol(gtex.validation )] %>% data.matrix()

gtex.training.y = gtex.training[,c("id", "class","type")]
gtex.training.x = gtex.training[,4:ncol(gtex.training)] %>% data.matrix()

gtex.training.x.noise = apply(gtex.training.x, 2, function(x) { noise = rnorm(length(x) , .1 , .05) ;  x + noise} ) 

archs.y = archs.gene250.wide[,c("id", "class","type")]
archs.x = archs.gene250.wide[,4:ncol(archs.gene250.wide)] %>% data.matrix()

autoencoder_model = gen_encoder_arch( inshape = ncol(gtex.training.x), bottle = 20) 
history = autoencoder_model %>% 
        fit(
            x = data.matrix(gtex.training.x.noise) , 
            y = data.matrix(gtex.training.x) , 
            epochs = 100, 
            batch_size =  32,
            validation_split = 0.2,
            view_metrics = FALSE,
        )
intermediate_layer_model <- keras_model(inputs = autoencoder_model$input, outputs = get_layer(autoencoder_model, "bottle")$output)

pred_val =  predict(autoencoder_model , gtex.validation.x)
se_val = apply((data.matrix(gtex.validation.x) - pred_val)^2, 1, sum  )
mse_val = se_val / ncol(pred_val)
gtex.val.mse = cbind(gtex.validation.y, mse_val) 


pred_archs =  predict(autoencoder_model , archs.x)
archs.se_val = apply((data.matrix(archs.x) - pred_archs)^2, 1, sum  )
archs.val.mse = archs.y
archs.val.mse$mse_val = archs.se_val / ncol(pred_archs)
archs.val.mse$id =  str_split_fixed(archs.val.mse$id, "\\-", 2)[,1] 

archs.gtex = rbind(archs.val.mse, gtex.val.mse)
cl = c(
  "#7e1e9c","#15b01a","#0343df","#ff81c0","#653700","#e50000","#95d0fc","#f97306","#029386","#c20078",
      "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
      "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00"
      )
library(ggplot2)
ggplot(archs.gtex, aes( factor(type), mse_val, color = factor(type))) + 
    geom_jitter( size = 2.5, width = .3, aes(alpha = .5)) + 
    theme_bw() + 
    geom_hline( yintercept = mean(gtex.val.mse$mse_val) + (2*sd(gtex.val.mse$mse_val)), color = "red", size = 1) +
    facet_wrap(~class,  scales = "free_x" , nrow = 3  ) +
    scale_color_manual(values = cl) + 
    theme(axis.text = element_text(size = 20), legend.position="bottom" ,
        strip.text.x = element_text(size =20))
sel = c("Spleen", "Ovary", "Muscle", "Pancreas", "Heart")
archs.gtex.sel = archs.gtex[archs.gtex$class %in% sel]
ggplot(archs.gtex.sel, aes( factor(type), mse_val, color = factor(type))) + 
    geom_jitter( size = 2.5, width = .3, aes(alpha = .5)) + 
    theme_bw() + 
    geom_hline( yintercept = mean(gtex.val.mse$mse_val) + (2*sd(gtex.val.mse$mse_val)), color = "red", size = 1) +
    facet_wrap(~class,  scales = "free_x" , nrow = 1  ) +
    scale_color_manual(values = cl) + 
    theme(axis.text = element_text(size = 20), legend.position="bottom" ,
        strip.text.x = element_text(size =20))
write.table(archs.gtex, file = paste0("reconstruction.error.", "gene", "archs.gtex.csv"), sep = ",", 
            col.names = F, quote = F, eol = "\n" , row.names = F )



##Figure 6B/S13 reconstructions errors with GES

gen_encoder_arch = function(inshape, bottle,  lossfx = "mean_squared_error")
{
    input_layer <- 
    layer_input(shape = inshape)  
    encoder <- 
        input_layer %>% 
        layer_dense(units = 100, activation = "relu") %>%
        layer_dense(units = 50, activation = "relu") %>%
        layer_dense(units = 25, activation = "relu") %>%
        layer_dense(units = bottle,  name = "bottle") # 3 dimensions for the output layer
    decoder <- 
        encoder %>% 
        layer_dense(units = 25, activation = "relu") %>% 
        layer_dense(units = 50, activation = "relu") %>% 
        layer_dense(units = 100, activation = "sigmoid") %>% 
        layer_dense(units = inshape)
    autoencoder_model <- keras_model(inputs = input_layer, outputs = decoder)
    autoencoder_model %>% compile(
        loss = lossfx, 
        optimizer = "adam"
    )
    return(autoencoder_model)
}  

ss_to_wide = function(data , type = "pvalue" , datatype  , bottle = 3,  epochs = 25 )
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
    names(data) = header_ss
    ref_header = data.table
    if ( datatype == "domain")
    {
        ref_header = fread("files/entry.list", sep = "\t" ,  header = F, stringsAsFactors= F )
        ref_header = ref_header[V2 == "Domain", V1 ] 
    }  else if ( datatype == "fold")
    {
        ref_header = fread("files/scope_total_2.06-stable.txt", sep = "|", header = F, stringsAsFactors = F, quote="")
        ref_header = ref_header[V1 == "folds", V2 ]
    }  else if ( datatype == "family")
    {
        ref_header = fread("files/scope_total_2.06-stable.txt", sep = "|", header = F, stringsAsFactors = F, quote="")
        ref_header = ref_header[V1 == "families", V2 ]
    }   else if ( datatype == "superfamily")
    {
        ref_header = fread("files/scope_total_2.06-stable.txt", sep = "|", header = F, stringsAsFactors = F, quote="")
        ref_header = ref_header[V1 == "superfamilies", V2 ]
    }
    if ( type == "presense")
    {
        data$presense = rep(1, nrow(data))  
    }
    datacol = data[, ..type] 
    data = data[,c("ID","class","structure")] 
    ####pvalues need to rescaled by their -log since the fill is 0 
    if ( type == "pvalue" )
    {
        data$datacol = -log10(datacol ) 
    }
    else 
    {
        data$datacol  = datacol
    }
    reference_table = cbind(rep("setupvariable", length(ref_header)) , 
        rep("setupvariable", length(ref_header)) , 
        ref_header , 
        rep(1, length(ref_header)))  %>% as.data.table()
    names(reference_table) = c("ID", "class", "structure", "datacol" )      
    data = data[which(data$structure %in% ref_header) , ] 
    data = rbind(reference_table, data) %>% as.data.table 
    print("Converting data from long to wide")
    data.wide  = spread(data, key = structure, fill = 0 , value = datacol )
    data.wide = data.wide[! which(data.wide$ID == "setupvariable"), ]
    row.names(data.wide) = data.wide$id 
    data.y = data.wide[,c("ID", "class")]
    data.x = data.wide[,3:ncol(data.wide)] %>% data.matrix()  %>% normalize( method = "range", range = c(0, 1), margin = 1) 
    return_data = list(data_x = data.x, data_y = data.y )
}

embed_ss_signatures = function(data.x , data.y , type = "pvalue" , datatype  , bottle = 3,  epochs = 25, noise_m = .01, noise_sd = .005 )
{
    print("Training autoencoder model")
    ##to make a stacked denoising autoencoder 
    data.x.noise = apply(data.x, 2, function(x) { noise = rnorm(length(x) , noise_m , noise_sd ) ;  x + noise} ) 
    autoencoder_model = gen_encoder_arch( inshape = ncol(data.x), bottle = bottle ) 
    history = autoencoder_model %>% 
        fit(
            x = data.matrix(data.x.noise) , 
            y = data.matrix(data.x) , 
            epochs = epochs, 
            batch_size =  32,
            validation_split = 0.2,
            view_metrics = FALSE,
        )
    header = names(data.x)
    intermediate_layer_model <- keras_model(inputs = autoencoder_model$input, outputs = get_layer(autoencoder_model, "bottle")$output)
    return_data = list(complete_model = autoencoder_model , embed_model =intermediate_layer_model , header = header , training_data = data.x , training_meta = data.y , params = c(type = type, bottleneck = bottle, epochs = epochs) , performance = history )
    return(return_data)
}
##re run with different structural levels 

gtex.ss250 = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.domain.csv", header = F , stringsAsFactors = F)
archs.ss250 = fread("figures/data/archs/structural-signatures/allcombined.archs.250.domain.csv", header = F , stringsAsFactors = F)
re.gene = fread("figures/data/reconstruction_errors/reconstruction.error.gene.archs.gtex.csv", 
            header = F , stringsAsFactors = F)
gtex.meta = str_split_fixed(gtex.ss250$V10, "\\.", 2) %>% as.data.frame
gtex.meta$V1 = str_split_fixed(gtex.meta$V1, "\\-", 2)[,1]
gtex.ss250$V10 = gtex.meta$V2
gtex.ss250$V11 = gtex.meta$V1
gtex.n = gtex.ss250$V11 %>% unique %>% as.character
archs.n = archs.ss250$V11 %>% unique %>% as.character
shared = gtex.n[which(gtex.n %in% archs.n)] 
gtex.ss250 = gtex.ss250[gtex.ss250$V11 %in% shared, ]
archs.ss250 = archs.ss250[archs.ss250$V11 %in% shared, ]

gtex.val.set = re.gene[V3 == "gtex", V1]
gtex.validation = gtex.ss250[gtex.ss250$V10 %in% gtex.val.set ]
gtex.training = gtex.ss250[! gtex.ss250$V10 %in% gtex.val.set ]

gtex.training.wide = ss_to_wide( gtex.training, type = "counts_observed",  datatype = "domain")

gtex.training.ae = embed_ss_signatures( gtex.training.wide$data_x ,     noise_m = .1, noise_sd = .05,  
                        gtex.training.wide$data_y, datatype = "domain",  bottle = 20, epochs = 100)


gtex.validation.wide = ss_to_wide( gtex.validation, type = "counts_observed" , datatype = "domain")
archs.wide = ss_to_wide( archs.ss250, type = "counts_observed" ,  datatype = "domain")


pred_val =  predict(gtex.training.ae$complete_model , gtex.validation.wide$data_x)
se_val = apply((data.matrix(gtex.validation.wide$data_x) - pred_val)^2, 1, sum  )
mse_val = se_val / ncol(pred_val)
gtex.val.mse = cbind(gtex.validation.wide$data_y, mse_val) 
gtex.val.mse$type = rep("gtex", nrow(gtex.val.mse))

pred_archs =  predict(gtex.training.ae$complete_model, archs.wide$data_x)
archs.se_val = apply((data.matrix(archs.wide$data_x) - pred_archs)^2, 1, sum  )
archs.val.mse = archs.wide$data_y
archs.val.mse$mse_val = archs.se_val / ncol(pred_archs)
archs.val.mse$type = rep("archs", nrow(archs.val.mse))
archs.gtex = rbind(archs.val.mse, gtex.val.mse)

cl = c(  "#7e1e9c","#15b01a","#0343df","#ff81c0","#653700","#e50000","#95d0fc","#f97306","#029386","#c20078",
      "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
      "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00" )
library(ggplot2)
ggplot(archs.gtex, aes( factor(type), mse_val, color = factor(type))) + 
    geom_jitter( size = 2.5, width = .3, aes(alpha = .5)) + 
    theme_bw() + 
    #geom_hline( yintercept = sd(archs.gtex$mse_val) + mean(archs.gtex$mse_val) , color = "red", size = 1) +
    facet_wrap(~class,  scales = "free_x" , nrow = 3  ) +
    scale_color_manual(values = cl) + 
    theme(axis.text = element_text(size = 20), legend.position="bottom" ,
        strip.text.x = element_text(size =20))

###Figure 6C/S14
###integrate reconstructions errors
library(data.table)
library(tidyverse)
library(keras)
library(BBmisc)  

re.gene = fread("figures/data/reconstruction_errors/reconstruction.error.gene.archs.gtex.csv", 
            header = F , stringsAsFactors = F)
re.domain = fread("figures/data/reconstruction_errors/reconstruction.error.domain.archs.gtex.csv", 
            header = F , stringsAsFactors = F)
re.family = fread("figures/data/reconstruction_errors/reconstruction.error.family.archs.gtex.csv", 
            header = F , stringsAsFactors = F)
re.superfamily = fread("figures/data/reconstruction_errors/reconstruction.error.superfamily.archs.gtex.csv", 
            header = F , stringsAsFactors = F)
re.fold = fread("figures/data/reconstruction_errors/reconstruction.error.fold.archs.gtex.csv", 
            header = F , stringsAsFactors = F)                                    
gene.err = re.gene$V4
re.gene$V4 = re.gene$V3
re.gene$V3 = gene.err
header = c("id", "class", "error" , "type")
names(re.gene) = header 
names(re.domain) = header 
names(re.family) = header 
names(re.superfamily) = header 
names(re.fold) = header 
re.gene$gene_norm_error = re.gene$error %>%  normalize( method = "range", range = c(0, 1), margin = 1)
re.domain$domain_norm_error = re.domain$error %>%  normalize( method = "range", range = c(0, 1), margin = 1)
re.family$family_norm_error = re.family$error %>%  normalize( method = "range", range = c(0, 1), margin = 1)
re.superfamily$superfamily_norm_error = re.superfamily$error %>%  normalize( method = "range", range = c(0, 1), margin = 1)
re.fold$fold_norm_error = re.fold$error %>%  normalize( method = "range", range = c(0, 1), margin = 1)

combined_error = Reduce(function (x,y) merge(x = x, y = y , by = "id"), list(
        re.gene, re.domain[,c(1,5)], re.family[,c(1,5)], re.superfamily[,c(1,5)], re.fold[,c(1,5)]))

combined_error$mean_normal_error =  apply(combined_error[,c(5:9)], 1 , mean) 
cl = c(
  "#7e1e9c","#15b01a","#0343df","#ff81c0","#653700","#e50000","#95d0fc","#f97306","#029386","#c20078",
      "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
      "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00"
      )
ggplot(combined_error, aes( factor(type), mean_normal_error, color = factor(type))) + 
    #geom_violin(aes(alpha = .5)) + 
    geom_jitter( size = 2.5, width = .3, aes(alpha = .5)) + 
    theme_bw() + 
    #geom_hline( yintercept = .0125, color = "red", size = 1) +
    facet_wrap(~class,  scales = "free_x" , nrow = 3  ) +
    scale_color_manual(values = cl) + 
    theme(axis.text = element_text(size = 20), legend.position="bottom" ,
        strip.text.x = element_text(size =20))


####select outlier archs samples

combined_error.sp = split(combined_error, f = combined_error$type)
combined_error.sp.gtex = split(combined_error.sp$gtex , combined_error.sp$gtex$class )
combined_error.sp.archs = split(combined_error.sp$archs , combined_error.sp$archs$class )
gtex.error.sd = lapply(combined_error.sp.gtex, function(x) { 
    1*sd(x$mean_normal_error) + mean(x$mean_normal_error)})

outliers.id.all = data.frame()
remain.id.all = data.frame()
percent_cnt = data.frame()
for (t in 1:length(gtex.error.sd))
{
    tiss = names(gtex.error.sd[t])
    outliers = combined_error.sp.archs[[t]][which(combined_error.sp.archs[[t]]$mean_normal_error > gtex.error.sd[[t]] ),]
    outliers.id = outliers[,1:2] %>% unique  
    remain = combined_error.sp.archs[[t]][which(combined_error.sp.archs[[t]]$mean_normal_error <= gtex.error.sd[[t]] ),]
    remain.id = remain[,1:2] %>% unique  
    outlier_cnt = nrow(outliers[,1])
    remain_cnt = nrow(remain[,1])
    total_cnt = nrow(combined_error.sp.archs[[t]])
    pnt_out = (outlier_cnt/ total_cnt ) * 100
    pnt_remain = (remain_cnt/ total_cnt) * 100
    df = c(tiss,total_cnt,remain_cnt,outlier_cnt,pnt_out,pnt_remain) %>% t() %>% as.data.frame
    names(df) = c("tissue", "total.counts","remaining", "outliers", "percent_outliers","percent_remain")
    percent_cnt =  rbind(percent_cnt, df)
    outliers.id.all = rbind(outliers.id.all, outliers.id)
    remain.id.all = rbind(remain.id.all, remain.id)
}

archs.outliers = outliers.id.all$id 
archs.remain = remain.id.all$id 

write.table(archs.outliers, file = "archs.outliers.id.2sd.csv", sep = ",",
             col.names = F, quote = F, row.names = F)
write.table(archs.remain, file = "archs.remain.id.2sd.csv", sep = ",",
             col.names = F, quote = F, row.names = F)


###Figure 6D/S15 Performance of a classifier before and after Outlier removal 

###Figure S17 
###ARchs consistancy with /wo outliers 
library(data.table)
library(magrittr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ranger)
library(pROC)
library(plyr)
library(dplyr)
library(stringr)
gtex.ss250 = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.fold.csv", header = F , stringsAsFactors = F)
archs.ss250 = fread("figures/data/archs/structural-signatures/allcombined.archs.250.fold.csv", header = F , stringsAsFactors = F)
re.gene = fread("figures/data/reconstruction_errors/reconstruction.error.gene.archs.gtex.csv", 
            header = F , stringsAsFactors = F)
gtex.meta = str_split_fixed(gtex.ss250$V10, "\\.", 2) %>% as.data.frame
gtex.meta$V1 = str_split_fixed(gtex.meta$V1, "\\-", 2)[,1]
gtex.ss250$V10 = gtex.meta$V2
gtex.ss250$V11 = gtex.meta$V1
gtex.n = gtex.ss250$V11 %>% unique %>% as.character
archs.n = archs.ss250$V11 %>% unique %>% as.character
shared = gtex.n[which(gtex.n %in% archs.n)] 
gtex.ss250 = gtex.ss250[gtex.ss250$V11 %in% shared, ]
archs.ss250 = archs.ss250[archs.ss250$V11 %in% shared, ]

gtex.gene250 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.250.csv", header = F , stringsAsFactors = F)
archs.gene250 = fread("figures/data/archs/archs.tissue.genelist.formated", header = F , stringsAsFactors = F)

colnames(gtex.gene250) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(archs.gene250) =c("SID", "Tissue", "Gene", "Rank")

gtex.gene250$size = rep(250, nrow(gtex.gene250))
archs.gene250$size = rep(250, nrow(archs.gene250))
archs.outliers = read.table("archs.outliers.id.2sd.csv", sep =",")
archs.gene250$gsm = str_split_fixed(archs.gene250$SID, "-", 2)[,1]
archs.gene250.no.outliers = archs.gene250[! archs.gene250$gsm %in% archs.outliers$V1,]
archs.gene250.no.outliers$gsm = NULL 
archs.gene250$gsm = NULL
archs.gene250.sp = split(archs.gene250, f= archs.gene250$`Tissue` )
archs.gene250.no.outliers.sp = split(archs.gene250.no.outliers, f= archs.gene250.no.outliers$`Tissue` )
compute_jaccard  = function(df.sp ){
    df.final = data.frame()
    for ( i in df.sp)
    {
        i$val = rep(1, nrow(i))
        size = i$size %>%  unique 
        tiss = i$`Tissue` %>%  unique  
        overalltissue = i$Tissue %>% unique
        print(tiss)
        i = i[,c("SID","Gene", "size", "val")]
        i.wide = spread(data = i, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
        row.names(i.wide) = i.wide$SID 
        i.wide$SID = NULL
        i.dist = vegdist(data.matrix(i.wide), method = "jaccard", binary = T ,
                         upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
        distances = i.dist[upper.tri(i.dist)] %>%  as.numeric() 
        distances = 1 - distances 
        df = data.frame()
        df = cbind(rep(tiss, length(distances)), distances) %>%  as.data.frame()
        df$size = rep(size , nrow(df))
        df$tissue = rep(overalltissue , nrow(df))
        df.final = rbind(df.final, df)
    }
    df.final$V1 = df.final$V1 %>%  as.factor()
    df.final$distances = df.final$distances %>%  as.character() %>% as.numeric()
    df.final$size = df.final$size %>%  as.factor()
    return(df.final)
}

archs.df.250 = compute_jaccard(archs.gene250.sp)
archs.df.250.no.outliers = compute_jaccard(archs.gene250.no.outliers.sp )
archs.df.250$type = rep("ARCHS4 with outliers", nrow(archs.df.250))
gtex.gene250.sp = split(gtex.gene250, f= gtex.gene250$`Tissue` )

archs.df.250 = compute_jaccard(archs.gene250.sp)
archs.df.no.outlier.250 = compute_jaccard(archs.gene250.no.outliers.sp)
archs.df.no.outlier.250$type = rep("ARCHS4 without outliers", nrow(archs.df.no.outlier.250))
archs.df.250$type = rep("ARCHS4 with outliers", nrow(archs.df.250))
df.final = rbind( archs.df.250, archs.df.no.outlier.250)
#write.table(df.final, file = "Jc-gtex-archs-pairwise-comparisons.csv", sep = ",", row.names = F, col.names = T, quote = F)

cl = c(  "#7e1e9c","#15b01a","#0343df","#ff81c0","#653700","#e50000","#95d0fc","#f97306","#029386","#c20078",
      "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
      "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00"    )

archs.n = names(archs.gene250.sp)

archs.no.outliers.n = names(archs.gene250.no.outliers.sp)

shared = gtex.n[which(gtex.n %in% archs.no.outliers.n)]
df.final = df.final[df.final$V1 %in% shared, ]
ggplot(df.final, 
       aes(x = factor(type), 
           y = as.numeric(distances), fill = factor(type))) +
    geom_hline(yintercept = .25 , color ="red", size = .75) + 
    facet_wrap( ~V1 , scales = "free_x", ncol = 8 ) + 
    geom_violin(aes(alpha = .5)) + 
    #scale_fill_manual(values = cl) + 
    scale_x_discrete("Gene list size") + 
    scale_y_continuous("Jaccard Coeficient") + 
    theme_bw() + 
    theme(axis.title = element_text(size = 15), 
          panel.grid = element_line(size = 1),
          axis.text.x = element_blank() 
          )

# gtex.x.archs = jc_gtex_X_archs(gtex.gene250.sp , archs.gene250.sp)
# gtex.x.archs.no.outliers = jc_gtex_X_archs(gtex.gene250.sp , archs.gene250.no.outliers.sp)
# gtex.x.archs$type = rep("GTEX X ACRCHS4 with outliers", nrow(gtex.x.archs))
# gtex.x.archs.no.outliers$type = rep("GTEX X ACRCHS4 without outliers", nrow(gtex.x.archs.no.outliers))
# gtex.x.archs$iter = NULL
# gtex.x.archs.no.outliers$iter = NULL


###Figure 6E 
##consistancy of overlap between archs and gtex before and after outlier removal 

archs.outliers = read.table("archs.outliers.id.2sd.csv", sep =",")

gtex.ss250.f = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.fold.csv", header = F , stringsAsFactors = F)
archs.ss250.f = fread("figures/data/archs/structural-signatures/allcombined.archs.250.fold.csv", header = F , stringsAsFactors = F)

gtex.ss250.sf = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.superfam.csv", header = F , stringsAsFactors = F)
archs.ss250.sf = fread("figures/data/archs/structural-signatures/allcombined.archs.250.superfam.csv", header = F , stringsAsFactors = F)

gtex.ss250.fm = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.family.csv", header = F , stringsAsFactors = F)
archs.ss250.fm = fread("figures/data/archs/structural-signatures/allcombined.archs.250.family.csv", header = F , stringsAsFactors = F)

gtex.ss250.d = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.domain.csv", header = F , stringsAsFactors = F)
archs.ss250.d = fread("figures/data/archs/structural-signatures/allcombined.archs.250.domain.csv", header = F , stringsAsFactors = F)


gtex.meta.f = str_split_fixed(gtex.ss250.f$V10, "\\.", 2) %>% as.data.frame
gtex.meta.f$V1 = str_split_fixed(gtex.meta.f$V1, "\\-", 2)[,1]
gtex.meta.sf = str_split_fixed(gtex.ss250.sf$V10, "\\.", 2) %>% as.data.frame
gtex.meta.sf$V1 = str_split_fixed(gtex.meta.sf$V1, "\\-", 2)[,1]
gtex.meta.fm = str_split_fixed(gtex.ss250.fm$V10, "\\.", 2) %>% as.data.frame
gtex.meta.fm$V1 = str_split_fixed(gtex.meta.fm$V1, "\\-", 2)[,1]
gtex.meta.d = str_split_fixed(gtex.ss250.d$V10, "\\.", 2) %>% as.data.frame
gtex.meta.d$V1 = str_split_fixed(gtex.meta.d$V1, "\\-", 2)[,1]

gtex.ss250.f$V10 = gtex.meta$V2
gtex.ss250.f$V11 = gtex.meta$V1
gtex.ss250.sf$V10 = gtex.meta.sf$V2
gtex.ss250.sf$V11 = gtex.meta.sf$V1
gtex.ss250.fm$V10 = gtex.meta.fm$V2
gtex.ss250.fm$V11 = gtex.meta.fm$V1
gtex.ss250.d$V10 = gtex.meta.d$V2
gtex.ss250.d$V11 = gtex.meta.d$V1

gtex.ss250.f = gtex.ss250.f[,c(10,11,1)]
archs.ss250.f = archs.ss250.f[,c(10,11,1)]
gtex.ss250.sf = gtex.ss250.sf[,c(10,11,1)]
archs.ss250.sf = archs.ss250.sf[,c(10,11,1)]
gtex.ss250.fm = gtex.ss250.fm[,c(10,11,1)]
archs.ss250.fm = archs.ss250.fm[,c(10,11,1)]
gtex.ss250.d = gtex.ss250.d[,c(10,11,1)]
archs.ss250.d = archs.ss250.d[,c(10,11,1)]


gtex.ss250.f$rank = rep(1, nrow(gtex.ss250.f))
archs.ss250.f$rank = rep(1, nrow(archs.ss250.f))
gtex.ss250.sf$rank = rep(1, nrow(gtex.ss250.sf))
archs.ss250.sf$rank = rep(1, nrow(archs.ss250.sf))
gtex.ss250.fm$rank = rep(1, nrow(gtex.ss250.fm))
archs.ss250.fm$rank = rep(1, nrow(archs.ss250.fm))
gtex.ss250.d$rank = rep(1, nrow(gtex.ss250.d))
archs.ss250.d$rank = rep(1, nrow(archs.ss250.d))

gtex.ss250.f$size = rep(250, nrow(gtex.ss250.f))
archs.ss250.f$size = rep(250, nrow(archs.ss250.f))
gtex.ss250.sf$size = rep(250, nrow(gtex.ss250.sf))
archs.ss250.sf$size = rep(250, nrow(archs.ss250.sf))
gtex.ss250.fm$size = rep(250, nrow(gtex.ss250.fm))
archs.ss250.fm$size = rep(250, nrow(archs.ss250.fm))
gtex.ss250.d$size = rep(250, nrow(gtex.ss250.d))
archs.ss250.d$size = rep(250, nrow(archs.ss250.d))

names(gtex.ss250.f) = c("SID", "Tissue", "Gene", "Rank", "size")
names(archs.ss250.f) = c("SID", "Tissue", "Gene", "Rank", "size")
names(gtex.ss250.sf) = c("SID", "Tissue", "Gene", "Rank", "size")
names(archs.ss250.sf) = c("SID", "Tissue", "Gene", "Rank", "size")
names(gtex.ss250.fm) = c("SID", "Tissue", "Gene", "Rank", "size")
names(archs.ss250.fm) = c("SID", "Tissue", "Gene", "Rank", "size")
names(gtex.ss250.d) = c("SID", "Tissue", "Gene", "Rank", "size")
names(archs.ss250.d) = c("SID", "Tissue", "Gene", "Rank", "size")

##parse outliers 
archs.gene250$SID = str_split_fixed(archs.gene250$SID, "-", 2)[,1]
archs.gene250.n.outliers = archs.gene250[! archs.gene250$SID %in% archs.outliers$V1]
archs.ss250.f.n.outliers = archs.ss250.f[! archs.ss250.f$SID %in% archs.outliers$V1]
archs.ss250.sf.n.outliers = archs.ss250.sf[! archs.ss250.sf$SID %in% archs.outliers$V1 ]
archs.ss250.fm.n.outliers = archs.ss250.fm[! archs.ss250.fm$SID %in% archs.outliers$V1]
archs.ss250.d.n.outliers = archs.ss250.d[! archs.ss250.d$SID %in% archs.outliers$V1]

archs.ss250.f.sp = split(archs.ss250.f, f= archs.ss250.f$`Tissue` )
gtex.ss250.f.sp = split(gtex.ss250.f, f= gtex.ss250.f$`Tissue` )
archs.ss250.sf.sp = split(archs.ss250.sf, f= archs.ss250.sf$`Tissue` )
gtex.ss250.sf.sp = split(gtex.ss250.sf, f= gtex.ss250.sf$`Tissue` )
archs.ss250.fm.sp = split(archs.ss250.fm, f= archs.ss250.fm$`Tissue` )
gtex.ss250.fm.sp = split(gtex.ss250.fm, f= gtex.ss250.fm$`Tissue` )
archs.ss250.d.sp = split(archs.ss250.d, f= archs.ss250.d$`Tissue` )
gtex.ss250.d.sp = split(gtex.ss250.d, f= gtex.ss250.d$`Tissue` )
archs.gene250.sp = split(archs.gene250, f= archs.gene250$`Tissue` )
gtex.gene250.sp = split(gtex.gene250, f= gtex.gene250$`Tissue` )

archs.gene250.n.outliers.sp = split(archs.gene250.n.outliers , f= archs.gene250.n.outliers$`Tissue` )
archs.ss250.f.n.outliers.sp = split(archs.ss250.f.n.outliers , f= archs.ss250.f.n.outliers$`Tissue` )
archs.ss250.sf.n.outliers.sp = split(archs.ss250.sf.n.outliers , f= archs.ss250.sf.n.outliers$`Tissue` )
archs.ss250.fm.n.outliers.sp = split(archs.ss250.fm.n.outliers , f= archs.ss250.fm.n.outliers$`Tissue` )
archs.ss250.d.n.outliers.sp = split(archs.ss250.d.n.outliers , f= archs.ss250.d.n.outliers$`Tissue` )
### overlap calc 


jc_gtex_X_archs = function(gtex.sp , archs.sp )
{
  gtex.n = names(gtex.sp)
  archs.n = names(archs.sp)
  shared = gtex.n[which(gtex.n %in% archs.n)]
  df.final = data.frame()
  for (tiss in shared ) 
  {
    print(paste0("working on ", tiss))
    gtex.t = gtex.sp[[tiss]]
    archs.t = archs.sp[[tiss]]
    gtex.t$val = rep(1, nrow(gtex.t))
    archs.t$val = rep(1, nrow(archs.t))
    gtex.t = gtex.t[,c("SID", "Gene", "size", "val")]
    archs.t = archs.t[,c("SID", "Gene", "size", "val")]
    total = rbind(gtex.t, archs.t)
    gtex.sid = gtex.t$SID %>% unique 
    archs.sid = archs.t$SID %>% unique 
    total.wide = spread(data = total, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
    row.names(total.wide) = total.wide$SID 
    total.wide$SID = NULL
    size = total.wide$size %>% unique 
    total.wide$size = NULL
    total.dist = vegdist(data.matrix(total.wide), method = "jaccard", binary = T ,
                         upper = T, diag = T )  %>%  data.matrix() %>%  as.data.frame() 
    distances = total.dist[row.names(total.dist) %in% gtex.sid , names(total.dist) %in% archs.sid ] %>% as.matrix %>% as.numeric
    distances = 1 - distances  
    df.return = data.frame()
    df.return = cbind(rep(tiss, length(distances)), distances) %>%  as.data.frame()
    df.return$size = rep(size , nrow(df.return))
    df.return$tissue = rep(tiss, nrow(df.return))
    df.return$iter = rep(1 , nrow(df.return))
    df.final = rbind(df.final, df.return)    
  }
  return(df.final)
}

gtex.x.atchs.g = jc_gtex_X_archs(archs.gene250.sp  , gtex.gene250.sp)
gtex.x.atchs.f = jc_gtex_X_archs(archs.ss250.f.sp  , gtex.ss250.f.sp)
gtex.x.atchs.sf = jc_gtex_X_archs(archs.ss250.sf.sp  , gtex.ss250.sf.sp)
gtex.x.atchs.fm = jc_gtex_X_archs(archs.ss250.fm.sp  , gtex.ss250.fm.sp)
gtex.x.atchs.d = jc_gtex_X_archs(archs.ss250.d.sp  , gtex.ss250.d.sp)

gtex.x.atchs.g.n.out = jc_gtex_X_archs(archs.gene250.n.outliers.sp  , gtex.gene250.sp)
gtex.x.atchs.f.n.out = jc_gtex_X_archs(archs.ss250.f.n.outliers.sp  , gtex.ss250.f.sp)
gtex.x.atchs.sf.n.out = jc_gtex_X_archs(archs.ss250.sf.n.outliers.sp , gtex.ss250.sf.sp)
gtex.x.atchs.fm.n.out = jc_gtex_X_archs(archs.ss250.fm.n.outliers.sp , gtex.ss250.fm.sp)
gtex.x.atchs.d.n.out = jc_gtex_X_archs(archs.ss250.d.n.outliers.sp , gtex.ss250.d.sp)


gtex.x.atchs.g$type = rep("gene", nrow(gtex.x.atchs.g))
gtex.x.atchs.f$type = rep("fold", nrow(gtex.x.atchs.f))
gtex.x.atchs.fm$type = rep("family", nrow(gtex.x.atchs.fm))
gtex.x.atchs.sf$type = rep("superfamily", nrow(gtex.x.atchs.sf))
gtex.x.atchs.d$type = rep("domain", nrow(gtex.x.atchs.d))

gtex.x.atchs.g$type2 = rep("WithOutliers", nrow(gtex.x.atchs.g))
gtex.x.atchs.f$type2 = rep("WithOutliers", nrow(gtex.x.atchs.f))
gtex.x.atchs.fm$type2 = rep("WithOutliers", nrow(gtex.x.atchs.fm))
gtex.x.atchs.sf$type2 = rep("WithOutliers", nrow(gtex.x.atchs.sf))
gtex.x.atchs.d$type2 = rep("WithOutliers", nrow(gtex.x.atchs.d))


gtex.x.atchs.g.n.out$type = rep("gene", nrow(gtex.x.atchs.g.n.out))
gtex.x.atchs.f.n.out$type = rep("fold", nrow(gtex.x.atchs.f.n.out))
gtex.x.atchs.sf.n.out$type = rep("superfamily", nrow(gtex.x.atchs.sf.n.out))
gtex.x.atchs.fm.n.out$type = rep("family", nrow(gtex.x.atchs.fm.n.out))
gtex.x.atchs.d.n.out$type = rep("domain", nrow(gtex.x.atchs.d.n.out))

gtex.x.atchs.g.n.out$type2 = rep("WithoutOutliers", nrow(gtex.x.atchs.g.n.out))
gtex.x.atchs.f.n.out$type2 = rep("WithoutOutliers", nrow(gtex.x.atchs.f.n.out))
gtex.x.atchs.sf.n.out$type2 = rep("WithoutOutliers", nrow(gtex.x.atchs.sf.n.out))
gtex.x.atchs.fm.n.out$type2 = rep("WithoutOutliers", nrow(gtex.x.atchs.fm.n.out))
gtex.x.atchs.d.n.out$type2 = rep("WithoutOutliers", nrow(gtex.x.atchs.d.n.out))


df.final = rbind( gtex.x.atchs.g,
                  gtex.x.atchs.d, 
                  gtex.x.atchs.fm, 
                  gtex.x.atchs.sf, 
                  gtex.x.atchs.f,
                  gtex.x.atchs.g.n.out,
                  gtex.x.atchs.f.n.out,
                  gtex.x.atchs.sf.n.out,
                  gtex.x.atchs.fm.n.out,
                  gtex.x.atchs.d.n.out )

df.final = df.final[df.final$V1 %in% shared, ]
df.final$distances = df.final$distances %>% as.character %>% as.numeric
sel.tiss = c("Muscle", "WholeBlood")
shared = unique(archs.ss250.f.n.outliers$Tissue)
df.final.sel = df.final[df.final$V1 %in% sel.tiss,] 
cl = list(WithoutOutliers = "#54B8A3", WithOutliers = "#2E2E2E")
#cl = list( gene = "#51B7A1", domain = "#2274A5" , family = "#F75C03" , superfamily = "#00CB64" , fold = "#F5D34E"    )     

ggplot(df.final, 
       aes(x = factor(type), 
           y = as.numeric(distances), 
           fill = factor(type2))) +
    geom_hline(yintercept = .25 , color ="red", size = .75) + 
    facet_wrap( ~V1 , scales = "free_x", ncol = 2 ) + 
    geom_violin(aes(alpha = .5)) + 
    scale_fill_manual(values = cl) + 
    scale_x_discrete("Signature Type", limits = c("gene", "domain", "family", "superfamily", "fold")) + 
    scale_y_continuous("Jaccard Coeficient") + 
    theme_bw() + 
    geom_vline(xintercept=c(1.5,2.5,3.5,4.5),color="grey20")+
    theme(axis.title = element_blank(), 
        axis.text = element_text(size = 15), legend.position = "none")

ggsave(paste0("figures/S17-overlap-wo-outliers.png"), device =  "png", width = 8.5 , height = 11)


##determine concensus GES and sGES for muscle and whole blood tissue across ARCHS and GTEx
gtex.gene250.sp$WholeBlood$Gene 
archs.gene250.n.outliers.sp$WholeBlood$Gene

get_gtex_X_archs_signature = function(gtex.sp , archs.sp , pct_pres = .99 , type  )
{
  gtex.n = names(gtex.sp)
  archs.n = names(archs.sp)
  shared = gtex.n[which(gtex.n %in% archs.n)]
  df.final = data.frame()
  for (tiss in shared ) 
  {
    print(paste0("working on ", tiss))
    gtex.t = gtex.sp[[tiss]]
    archs.t = archs.sp[[tiss]]
    gtex.t$val = rep(1, nrow(gtex.t))
    archs.t$val = rep(1, nrow(archs.t))
    gtex.t = gtex.t[,c("SID", "Gene", "size", "val")]
    archs.t = archs.t[,c("SID", "Gene", "size", "val")]
    total = rbind(gtex.t, archs.t)
    gtex.sid = gtex.t$SID %>% unique 
    archs.sid = archs.t$SID %>% unique 
    total.wide = spread(data = total, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
    row.names(total.wide) = total.wide$SID 
    total.wide$SID = NULL
    size = total.wide$size %>% unique 
    total.wide$size = NULL
    total = nrow(total.wide)
    observed.freq = apply(total.wide, 2 , FUN = function(x) {sum(x)/ total })     
    signature = names(which(observed.freq >= pct_pres )) %>% as.data.frame 
    signature$tissue = rep(tiss, nrow(signature))
    signature$sigtype = rep(type, nrow(signature))
    names(signature) = c("signature","tissue")
    df.final = rbind(df.final, signature)
  }
  return(df.final)
}

signature.gene =  get_gtex_X_archs_signature(gtex.gene250.sp, archs.gene250.n.outliers.sp , 1, "gene")
signature.fold = get_gtex_X_archs_signature(gtex.ss250.f.sp, archs.ss250.f.n.outliers.sp, 1 , "fold")
signature.spr = get_gtex_X_archs_signature(gtex.ss250.sf.sp , archs.ss250.sf.n.outliers.sp, 1, "superfamily")
signature.fm = get_gtex_X_archs_signature(gtex.ss250.fm.sp , archs.ss250.fm.n.outliers.sp, 1, "family")
signature.domain = get_gtex_X_archs_signature(gtex.ss250.d.sp , archs.ss250.d.n.outliers.sp, 1 , "domain" )
