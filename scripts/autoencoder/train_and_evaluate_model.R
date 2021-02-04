#!/usr/bin/Rscript  


args <- commandArgs(TRUE)

if (length(args) < 6 ) {
stop("\n\033[31mThe following inputs are required\033[0m:
[1] job name
[2] training data (assumes a csv with ID and class column + n features )
[3] testing data (assumes a csv with ID and class column + n features )
[4] Output model hdf5?
    [a] yes or [b] no
[5] Output predictions?
    [a] yes or [b] no
[6] Generate confusion matrix? 
    [a] yes or [b] no
[7] Set common class? (That is training set will only have classes observed in the test set)
    [a] yes or [b] no
\033[34mGenerated output: 
    \033[33m'Job Name'.rocdata.csv 
\033[34mAdditional outputs:
    \033[35m'Job Name'.hdfs + (optional) 
    'Job name'.predictions.csv + (optional)
    'Job name'.confusion.csv (optional) \033[0m", call.=FALSE)
} 


library(tidyr)
library(dplyr)
library(ggridges)
library(readr)
library(purrr)
library(magrittr)
library(tidyverse)
library(keras)
library(ranger)
library(pROC)
library(data.table)
library(caret)

# args[1] = "go-tissue-predictions"
# args[2] = "embed_models/archs.combined.training.csv"
# args[3] = "go-analysis/go-analysis.combined.signatures.csv"
# args[4] = "yes"
# args[5] = "yes"
# args[6] = "yes"
# args[7] = "yes"
getroc_data = function(x, type  )
{
  x.roc.sp = split(x, x$`observed_label`)
  df = data.frame()
  for ( num in 1:length(x.roc.sp))
  {
    i = x.roc.sp[[num]]
    class = unique(i$`observed_label`) %>%  as.character()
    prob.positive = i[, c( which(names(i) == class ))] 
    negativecases = x[! x$observed_label == class, ]
    prob.negative = negativecases[, c( which(names(negativecases) == class ))]
    
    roc_dat = cbind(c(rep(1, length(prob.positive)),rep(0, length(prob.negative)) ) ,
                    c(prob.positive, prob.negative)) %>%  as.data.frame()
    names(roc_dat) = c("Observed", "Predicted")
    print(class)
    roc.dat = roc(roc_dat$Observed, roc_dat$Predicted)
    df.to.return = cbind(roc.dat$sensitivities , roc.dat$specificities) %>%  as.data.frame()
    names(df.to.return) = c("sensitivity", "specificity")
    df.to.return$comparison = rep( type, nrow(df.to.return))
    df.to.return$class = rep( class, nrow(df.to.return))
    df.to.return$auc = rep( roc.dat$auc %>%  as.numeric() , nrow(df.to.return))
    df = rbind(df, df.to.return)
  }
  return(df)
}

jn = as.character(args[1])
training = fread(args[2], header = T, stringsAsFactors =F , sep =",")
testing = fread(args[3], header = T, stringsAsFactors =F , sep =",")
h5 = as.character(args[4])
savpred = as.character(args[5])
conf = as.character(args[6])
common = as.character(args[7])
if (common == "yes")
{
    training = training[which(training$class %in% testing$class ), ] 
    testing = testing[which(testing$class %in% training$class ), ] 
}
train.y = training[,c("ID","class")]
train.x = training[,-c("ID","class")]
test.y = testing[,c("ID","class")]
test.x = testing[,-c("ID","class")]

train.classes = levels(as.factor(train.y$class)) 
train.y.n = as.factor(train.y$class) %>% as.numeric
train.y.oh = to_categorical(train.y.n) %>% as.data.frame
train.y.oh$V1 = NULL
names(train.y.oh) = train.classes

model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 100, input_shape = c(ncol(train.x))) %>% 
  layer_activation('relu') %>% 
  layer_dropout(0.) %>%
  layer_dense(units = ncol(train.y.oh)) %>% 
  layer_activation('softmax')

model %>% compile(
  optimizer = 'adam',
  loss = 'categorical_crossentropy',
  metrics = c('accuracy')
)
history = model %>% 
        fit(
            x = data.matrix(train.x) , 
            y = data.matrix(train.y.oh) , 
            epochs = 100, 
            batch_size =  100,
            validation_split = 0.1,
            view_metrics = FALSE,
        )
pred = predict(model, data.matrix(test.x )) %>% as.data.frame
comp.cases = ! rowSums(is.na(pred)) > 0
pred = pred[comp.cases, ]
names(pred) = train.classes
prob_class = apply(pred, 1 , max)
pred$predicted_label = names(pred)[ apply(pred, 1 , which.max) ]
pred$observed_label = test.y[comp.cases,]$class
pred$predicted_prob = prob_class 
pred.conf = pred[,c("predicted_label", "observed_label")] %>% table() %>% data.matrix

pred.roc = getroc_data(pred, jn)
write.table(pred.roc, paste0(jn,".rocdata.csv"), col.names = T, row.names = F, quote = F, sep = "," , eol ="\n")


if (h5 == "yes" )
{
    save_model_hdf5(model, paste0(jn,".h5"), overwrite = TRUE)
}

if (savpred == "yes" ) 
{
    write.table(pred, paste0(jn,".predictions.csv"), col.names = T, row.names = F, quote = F, sep = "," , eol ="\n")
}
if (conf == "yes")
{
    write.table(pred.conf, paste0(jn,".confusion.csv"), col.names = T, row.names = T, quote = F, sep = "," , eol ="\n")
}
