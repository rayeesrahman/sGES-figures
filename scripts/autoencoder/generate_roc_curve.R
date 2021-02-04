#!/usr/bin/Rscript  
args <- commandArgs(TRUE)
if (length(args) < 1 ) {
  stop("\n\033[31mThe following inputs are required\033[0m:\n[1] 'Job Name'.rocdata.csv")
} 


library(ggplot2)
library(data.table)
library(tidyverse)
library(tcltk)
X11()

args[1] = "gtex.archs.complete.rocdata.csv.ss"
roc.1= fread("generate_autoencoder_models/archs-no-outliers.model.rocdata.csv", header =T, stringsAsFactors = F, sep =",")
roc.2= fread("generate_autoencoder_models/archs-with-outliers.model.rocdata.csv", header =T, stringsAsFactors = F, sep =",")
roc.1$comparison = rep("Archs outliers removed", nrow(roc.1))
roc.2$comparison = rep("Archs outliers retained", nrow(roc.2))
roc.2 = roc.2[which(roc.2$class %in% roc.1$class), ] 
roc.total = rbind(roc.1, roc.2 )
names(roc.total) = c("sensitivity","specificity","comparison", "class","auc")
average.aucs = roc.total %>%  group_by(class, comparison ) %>%  summarise(average_auc=(mean(auc))) %>%  as.data.frame()

# prompt  <- "hit spacebar to close plots"
# extra   <- "some extra comment"

# capture <- tk_messageBox(message = prompt, detail = extra)


cl = c("#20A387FF", "black", "#2274A5", "#F75C03", "#F1C40F" , "#00CC66")
hist(roc.total$auc )
ggplot(roc.total, aes( 1- specificity ,sensitivity , color = factor(comparison) )) + 
    geom_line(size = 1.5, alpha = .75)   + 
    facet_wrap(~class, nrow = 4) +
    geom_abline( slope = 1, intercept = 0 , color ="red") + 
    geom_text(data = average.aucs[average.aucs$comparison=="Archs outliers removed",], aes(.60, .25, label = paste0("Outliers removed ", "AUC: ", 
                        sprintf('%.2f', round(average_auc, digits=2 ) ) ), 
                            group = class), color = "black", size = 6)   +
    geom_text(data = average.aucs[average.aucs$comparison=="Archs outliers retained",], aes(.60, .05, label = paste0("Outliers retained ", "AUC: ", 
                        sprintf('%.2f', round(average_auc, digits=2 ) ) ), 
                            group = class), color = "black", size = 6)   +
    
    theme_bw() +
    scale_color_manual(values = cl) + 
    theme(axis.title = element_text(size  = 20), 
    strip.text.x = element_text(size = 12)) 

        
cl = c("#20A387FF", "black", "#2274A5", "#F75C03", "#F1C40F" , "#00CC66")
 ggplot(roc.total, aes( 1- Specificity ,Sensitivity , color = factor(size) )) + 
  geom_line(size = 1.5, alpha = .75) + 
  facet_wrap(~tissue, nrow = 4) +
    geom_text(data = average.aucs[average.aucs$size=="domain", ], aes(.73, .50, label = paste0("AUC Domain", ": ", 
                                                                                            sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                    group = tissue), color = "black", size = 5) +
  geom_text(data = average.aucs[average.aucs$size=="family", ], aes(.745, .35, label = paste0("AUC Family", ": ", 
                                                                                             sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                    group = tissue), color = "black", size = 5) +
  geom_text(data = average.aucs[average.aucs$size=="superfamily", ], aes(.735, .2, label = paste0("AUC Sfamily", ": ", 
                                                                                          sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                  group = tissue), color = "black", size = 5) +
  geom_text(data = average.aucs[average.aucs$size=="fold", ], aes(.77, .05, label = paste0("AUC Fold", ": ", 
                                                                                           sprintf('%.2f',round(average_auc, digits=2) )), 
                                                                  group = tissue), color = "black", size = 5) +
  
  geom_abline( slope = 1, intercept = 0 , color ="red") + 
  theme_bw() +
  scale_color_manual(values = cl) + 
  theme(axis.title = element_text(size  = 20), 
        strip.text.x = element_text(size = 12)) 

         