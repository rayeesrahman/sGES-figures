rm(list = ls())
setwd("sGES-figures/")
library(ggplot2)
library(magrittr)
library(vegan)
library(tidyverse)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(ranger)
library(pROC)
## Figure 4 - ------------------------------------------------------------
### 10X CV of  genes and sGES

genes = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.250.csv", header = F , stringsAsFactors = F, data.table = F )
domain = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.domain.csv", sep = "," , stringsAsFactors = F,  header = F, data.table = F )
fold = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.fold.csv", sep = "," , stringsAsFactors = F,  header = F, data.table = F )
family = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.family.csv", sep = "," , stringsAsFactors = F,  header = F, data.table = F )
sfam = fread("figures/data/gtex/structural-signatures/allcombined.gtex.250.superfam.csv", sep = "," , stringsAsFactors = F, header = F, data.table = F  )
domain[, c("V11","V12") ] = str_split_fixed(domain$V10, pattern = '\\.' , n = 2) 
fold[, c("V11","V12") ] = str_split_fixed(fold$V10, pattern = '\\.' , n = 2) 
family[, c("V11","V12") ] = str_split_fixed(family$V10, pattern = '\\.' , n = 2) 
sfam[, c("V11","V12") ] = str_split_fixed(sfam$V10, pattern = '\\.' , n = 2) 

domain[, c("V10") ] = str_split_fixed(domain$V10, pattern = '\\-' , n = 2)[,1] %>%  str_replace(".GTEX", "")
fold[, c("V10") ] = str_split_fixed(fold$V10, pattern = '\\-' , n = 2)[,1] %>%  str_replace(".GTEX", "")
family[, c("V10") ] = str_split_fixed(family$V10, pattern = '\\-' , n = 2)[,1] %>%  str_replace(".GTEX", "")
sfam[, c("V10") ] = str_split_fixed(sfam$V10, pattern = '\\-' , n = 2)[,1] %>%  str_replace(".GTEX", "")
header = c(
    "Gene",
    "counts",
    "bg-counts",
    "genesetsize",
    "bg-genesetsize",
    "pvalue",
    "qvalue",
    "bonforroni",
    "logfoldchange",
    "Tissue",
    "Sub-Tissue",
    "SID"
)
names(domain) = header 
names(family) = header 
names(sfam) = header 
names(fold) = header 
names(genes) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")


domain = domain[,c("SID", "Sub-Tissue" , "Tissue", "Gene" )]
fold = fold[,c("SID", "Sub-Tissue" , "Tissue", "Gene")]
family = family[,c("SID", "Sub-Tissue" , "Tissue", "Gene")]
sfam = sfam[,c("SID", "Sub-Tissue" , "Tissue", "Gene")]

genes.sp = split(x = genes, f = genes$`Sub-Tissue`)
domain.sp = split(x = domain, f = domain$`Sub-Tissue`)
fold.sp = split(x = fold, f = fold$`Sub-Tissue`)
family.sp = split(x = family, f = family$`Sub-Tissue`)
sfam.sp = split(x = sfam, f = sfam$`Sub-Tissue`)

getroc_data = function(x, type  )
{
    x.roc.sp = split(x, x$`Sub-Tissue`)
    df = data.frame()
    for ( num in 1:length(x.roc.sp))
    {
        i = x.roc.sp[[num]]
        tiss = unique(i$`Sub-Tissue`)
        prob.positive = i[, c( which(names(i) == make.names(tiss) ))] 
        negativecases = x[! x$`Sub-Tissue` == make.names(tiss), ]
        prob.negative = negativecases[, c( which(names(negativecases) == make.names(tiss) ))]
        
        roc_dat = cbind(c(rep(1, length(prob.positive)),rep(0, length(prob.negative)) ) ,
                        c(prob.positive, prob.negative)) %>%  as.data.frame()
        names(roc_dat) = c("Observed", "Predicted")
        print(tiss)
        roc.dat = roc(roc_dat$Observed, roc_dat$Predicted)
        df.to.return = cbind(roc.dat$sensitivities , roc.dat$specificities) %>%  as.data.frame()
        df.to.return$size = rep( type, nrow(df.to.return))
        df.to.return$tissue = rep( tiss, nrow(df.to.return))
        df.to.return$auc = rep( roc.dat$auc %>%  as.numeric() , nrow(df.to.return))
        co = coords(roc.dat, "all", ret = c("recall", "precision"), transpose = T) %>% as.data.frame
        names(co) = c(1:ncol(co))
        co = data.frame(t(co))
        df.to.return$recall = co$recall 
        df.to.return$precision = co$precision
        df.to.return = df.to.return[complete.cases(df.to.return),] %>% as.data.frame
        df = rbind(df, df.to.return)
    }
    return(df)
}
generate_roc <- function(x, split = .3  , fold = 10 , type )
{
    x$pres = rep(1, nrow(x))
    #x$pres = NULL
    #x.wide = spread(x, key = Gene, value = Rank , fill = 0 )
    x$Rank = NULL
    x.wide = spread(x, key = Gene, value = pres , fill = 0 )
    x.wide.sp = split(x.wide, x.wide$`Sub-Tissue`)
    rocdata = data.frame()
    for (iter in 1:fold )
    {
        print(paste0("iteration ", iter))
        testing = data.frame()
        training = data.frame()
        ###get x% of training and testing cases, per tissue
        for ( i in x.wide.sp )
        {
            #print(unique(i$`Sub-Tissue`))
            row.names(i) = 1:nrow(i)
            len = 1:nrow(i)
            testing.cases = sort(sample(1:nrow(i), size = round(nrow(i)* split), 
                                        replace = F )) ##random sample of values of size x
            training.cases = len[! len %in% testing.cases ]
            training.sub = i[training.cases, ]
            testing.sub = i[ testing.cases, ]
            testing = rbind(testing,testing.sub ) 
            training = rbind(training, training.sub)
        }
        ##so that model doesnt learn order, everything is shuffled
        training = training[sample(1:nrow(training) , nrow(training), replace = F  ),]
        testing = testing[sample(1:nrow(testing) , nrow(testing), replace = F  ),]
        training.labels = training[,1:4]
        test.labels = testing[,1:4]
        training[,1:4] = NULL
        testing[,1:4] =  NULL 
        names(training) = names(training) %>%  make.names()
        names(testing) = names(training) %>%  make.names()
        training$tissue = training.labels$`Sub-Tissue` %>% make.names() %>%  factor()
        model = ranger( tissue ~ .  , data = training , mtry = 10, 
                        num.threads = 4, verbose = T , num.trees = 100, oob.error = T, 
                        probability = T )
        predictions = predict(model,data = testing )
        predictions = predictions$predictions %>%  as.data.frame()
        #    predict.votes = data.frame()
        #    for (i in 1:nrow(predictions))
        #    {
        #      rw = predictions[i,]
        #      predictedlabel = which.max(rw) %>%  names() 
        #      prob = max(rw)
        #      rw.2 = cbind(predictedlabel, prob)
        #      predict.votes = rbind(predict.votes, rw.2)
        #    }
        predictions = cbind(predictions, test.labels)
        rocvalues = getroc_data(predictions , type  )
        rocvalues$iter = rep(iter, nrow(rocvalues))
        rocdata = rbind(rocdata, rocvalues )
    }
    names(rocdata) = c("Sensitivity", "Specificity", "size", "tissue" , "auc", "iter")
    return(rocdata)
}

genes.roc = generate_roc(genes, split = .5, fold = 10, 250  )
fold.roc = generate_roc(fold, split = .5, fold = 10, "fold" )
family.roc = generate_roc(family, split = .5, fold = 10, "family"  )
domain.roc = generate_roc(domain, split = .5, fold = 10, "domain"  )
sfam.roc = generate_roc(sfam, split = .5, fold = 10, "sfam" )
roc.total = rbind(genes.roc, fold.roc, sfam.roc ,family.roc ,  domain.roc )
write.table(roc.total, "roc.total.presence.asbsence.structural.sigs.csv", sep = ",", quote = F, eol = "\n" ,row.names = F, 
            col.names = T )
roc.total$size = as.factor(roc.total$size)
roc.total$tissue = as.factor(roc.total$tissue)
roc.total$Sensitivity = as.numeric(as.character(roc.total$Sensitivity))
roc.total$Specificity = as.numeric(as.character(roc.total$Specificity))
roc.total$auc = as.numeric(as.character(roc.total$auc))
names(roc.total) = c("Sensitivity", "Specificity", "size", "tissue", "auc", "iter", "Recall", "Precision")
roc.total.sp = split(roc.total, f =  roc.total$tissue )
average.aucs = roc.total[roc.total$iter == 1, ] %>%  group_by(tissue, size ) %>%  summarise(average_auc=(mean(auc))) %>%  as.data.frame()
cl = c("#20A387FF", "#2274A5", "#F75C03", "#F1C40F" , "#00CC66")
##plot rocs 
p = ggplot(roc.total[roc.total$iter == 2, ], aes( 1- Specificity ,Sensitivity , color = factor(size) )) + 
    geom_line(size = 1.5, alpha = .75)  +
    geom_text(data = average.aucs[average.aucs$size==250, ], aes(.7, .65, label = paste0("AUC250", ": ", 
                                                                                          sprintf('%.2f', round(average_auc, digits=2 ) ) ), 
                                                                group = tissue), color = "black", size = 5) + 
    geom_text(data = average.aucs[average.aucs$size=="domain", ], aes(.6, .50, label = paste0("AUC Domain", ": ", 
                                                                                               sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                 group = tissue), color = "black", size = 5) +
    geom_text(data = average.aucs[average.aucs$size=="family", ], aes(.64, .35, label = paste0("AUC Family", ": ", 
                                                                                               sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                  group = tissue), color = "black", size = 5) +
    geom_text(data = average.aucs[average.aucs$size=="sfam", ], aes(.63, .2, label = paste0("AUC Sfamily", ": ", 
                                                                                            sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                      group = tissue), color = "black", size = 5) +
    geom_text(data = average.aucs[average.aucs$size=="fold", ], aes(.65, .05, label = paste0("AUC Fold", ": ", 
                                                                                             sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                    group = tissue), color = "black", size = 5) +
    facet_wrap(~tissue, nrow = 4) + 
    geom_abline( slope = 1, intercept = 0 , color ="red") + 
    theme_bw() +
    scale_color_manual(values = cl) + 
    theme(axis.title = element_text(size  = 20), 
          strip.text.x = element_text(size = 12)) 
g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g$layout$name))
cl = c("#7e1e9c","#15b01a","#0343df","#ff81c0","#653700","#e50000","#95d0fc","#f97306","#029386","#c20078",
       "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
       "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00")
fills =  c("#0a888a", "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00", 
           "#c20078", "#c20078","#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368",
           "#ff81c0","#653700","#653700","#e50000","#e50000","#95d0fc","#95d0fc","#ff5b00","#ff5b00","#ff5b00","#029386",
           "#7e1e9c","#15b01a","#15b01a","#15b01a","#0343df","#0343df","#0343df","#0343df","#0343df","#0343df","#0343df")
k <- 1
for (i in strip_both) {
    if ( g$grobs[[i]]$grobs %>% is.null() ) 
    {
        next()
    }
    else 
    {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        print(paste(i,  k,  j , sep =" "))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
        k <- k+1
    }
}
grid.draw(g)
