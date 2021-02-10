rm(list = ls())
setwd("sGES-figures/")
library(ggplot2)
library(magrittr)
library(vegan)
library(tidyverse)
library(ranger)
library(pROC)
library(plyr)
library(dplyr)
library(grid)
library(data.table)
library(Rtsne)
library(plotly)

#### read data 
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
    "Structure",
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
genes$Presense = rep(1, nrow(genes))
domain$Presense = rep(1, nrow(domain))
family$Presense = rep(1, nrow(family))
sfam$Presense = rep(1, nrow(sfam))
fold$Presense = rep(1, nrow(fold))
convert_to_wide = function(x, value_feat, key_feat)
{
    x = x[,c(key_feat, value_feat, "SID", "Tissue", "Sub-Tissue")]
    x.wide = spread(data = x, key = key_feat ,  value = value_feat,  fill = 0 ) %>%  as.data.frame()  
    row.names(x.wide) = x.wide$SID
    x.wide$SID = NULL
    return(x.wide)
}
plot_3d = function(tsne, cl , tissue, subtissue, out.table = F, table.name = "" )
{
    tsne.y = tsne$Y %>%  as.data.frame()
    tsne.y$subtissue = subtissue
    tsne.y$tissue = tissue
    p = plot_ly(tsne.y, 
        x = as.numeric(tsne.y$V1),
        y = as.numeric(tsne.y$V2), 
        z = as.numeric(tsne.y$V3), 
        color = factor(tsne.y$tissue) , 
        colors = cl ) %>% 
        layout(
            scene = list(
                xaxis = list(title = "TSNE1"),
                yaxis = list(title = "TSNE2"),
                zaxis = list(title = "TSNE3")
            )
        ) %>%  
        add_trace(
            text = paste(tsne.y$subtissue, "subtissue", sep = " "),
            hoverinfo = 'text',
            showlegend = T
        )
    if ( out.table == T ) 
    {
        if ( table.name == "")
        {
            print("No table name provided, using default name: tsne.out.csv")
            table.name = "tsne.out.csv"
        }
        write.table(table.name, table.name, 
                    sep = ",",  row.names = F, quote = F, 
                    eol = "\n" )
    }
    return(p)
}
cl = c(
    "#7e1e9c","#15b01a","#0343df","#ff81c0","#653700", "#e50000","#95d0fc","#f97306","#029386","#c20078",
       "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
       "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00"
       )

## Figure 21C -------------------------------------------------------------
### t-SNE clustering of GTeX tissues based on gene presense or absense 
genes.w = convert_to_wide(genes,  "Presense", "Gene")
g.tissues = genes.w$Tissue %>%  as.character()
g.subtissues = genes.w$`Sub-Tissue` %>% as.character()
genes.w$`Sub-Tissue` = NULL
genes.w$Tissue = NULL
genes.tsne = Rtsne( genes.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE, 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(genes.tsne, cl ,  g.tissues, g.subtissues, out.table = F)
### t-SNE clustering of GTeX tissues based on domain presense or absense
domain.w = convert_to_wide(domain,  "Presense", "Structure")
d.tissues = domain.w$Tissue %>%  as.character()
d.subtissues = domain.w$`Sub-Tissue` %>% as.character()
domain.w$`Sub-Tissue` = NULL
domain.w$Tissue = NULL
domain.tsne = Rtsne( domain.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE, 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(domain.tsne, cl ,  d.tissues, d.subtissues, out.table = F)
### t-SNE clustering of GTeX tissues based on fold presense or absense
fold.w = convert_to_wide(fold,  "Presense", "Structure")
f.tissues = fold.w$Tissue %>%  as.character()
f.subtissues = fold.w$`Sub-Tissue` %>% as.character()
fold.w$`Sub-Tissue` = NULL
fold.w$Tissue = NULL
fold.tsne = Rtsne( fold.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE,  check_duplicates = F , 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(fold.tsne, cl ,  f.tissues, f.subtissues, out.table = F)


## SFigure 2 -------------------------------------------------------------
### t-SNE clustering of GTeX tissues based on family and sfam presense or absense
family.w = convert_to_wide(family,  "Presense", "Structure")
family.tissues = family.w$Tissue %>%  as.character()
family.subtissues = family.w$`Sub-Tissue` %>% as.character()
family.w$`Sub-Tissue` = NULL
family.w$Tissue = NULL
family.tsne = Rtsne( family.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE,  check_duplicates = F , 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(family.tsne, cl ,  family.tissues, family.subtissues, out.table = F)

sfam.w = convert_to_wide(sfam,  "Presense", "Structure")
sfam.tissues = sfam.w$Tissue %>%  as.character()
sfam.subtissues = sfam.w$`Sub-Tissue` %>% as.character()
sfam.w$`Sub-Tissue` = NULL
sfam.w$Tissue = NULL
sfam.tsne = Rtsne( sfam.w, dims = 3, perplexity = 30 , 
                    partial_pca=TRUE,  check_duplicates = F , 
                    theta =.5 ,  max_iter = 1000, verbose = T )
plot_3d(sfam.tsne, cl ,  sfam.tissues, sfam.subtissues, out.table = F)