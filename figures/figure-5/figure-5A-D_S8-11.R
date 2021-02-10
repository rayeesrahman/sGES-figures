rm(list = ls())
setwd("sGES-figures/")
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
## read data -------------------------------------------------------------

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

## Figure 5A/S8 -------------------------------------------------------------
### Distributions of pairwise jaccard coefficients within tissues and across tissues 
#### compute pairwise jaccard coefficients
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

generate_random_Jc_distributions_per_tissue = function(data.sp, reps = 1, subset = 10){
  df.final = data.frame()
  for ( i in 1:reps )
  {
    print(paste0("bootstrap number ", i ))
    num.samp = 1:length(data.sp)
    samp.ids = lapply(data.sp, function(x) {
      unique(x$SID)
    }) 
    for (ii in 1:length(data.sp ))
    {
      ##get random sample of tissues equal to number of samples for a specific tissue 
      bootstrap.tissue = table(sample(num.samp[-ii], subset, replace = T)) %>% 
        as.data.frame
      bootstrap.data = apply(bootstrap.tissue, 1 , function(x) {
        b.tiss = x[1] %>% as.numeric 
        b.samps = x[2] %>% as.numeric
        ids = sample(samp.ids[[b.tiss]],b.samps, replace = T)  
        df = data.sp[[b.tiss]][SID %in% ids , ]
        return(df)
      })
      bootstrap.data = do.call(rbind, bootstrap.data)
      ##split by samples in the test data 
      z = data.sp[[ii]]$Tissue %>% unique 
      print(z)
      test = split(data.sp[[ii]], f = data.sp[[ii]]$SID) 
      test = test[unique(sample(length(test), subset, replace = T )) ] 
      for (iii in test)
      { 
        test.sid = iii$SID %>% unique 
        test.b.bind = rbind(iii,bootstrap.data)
        test.b.bind$val = rep(1, nrow(test.b.bind))
        size = iii$size %>%  unique 
        tiss = iii$`Tissue` %>%  unique  
        overalltissue = iii$Tissue %>% unique
        test.b.bind = test.b.bind[,c("SID", "Gene", "size", "val")]
        test.b.bind.wide = spread(data = test.b.bind, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
        row.names(test.b.bind.wide) = test.b.bind.wide$SID 
        test.b.bind.wide$SID = NULL
        test.b.bind.wide$size = NULL
        test.b.bind.dist = vegdist(data.matrix(test.b.bind.wide), method = "jaccard", binary = T ,
                         upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
        distances = test.b.bind.dist[row.names(test.b.bind.dist)== test.sid , names(test.b.bind.dist) != test.sid  ] %>% as.numeric
        distances = 1 - distances 
        df.return = data.frame()
        df.return = cbind(rep(tiss, length(distances)), distances) %>%  as.data.frame()
        df.return$size = rep(size , nrow(df.return))
        df.return$tissue = rep(z, nrow(df.return))
        df.return$iter = rep(i , nrow(df.return))
        df.final = rbind(df.final, df.return)
      }
    }
  }
  df.final$type = rep("random-bootstrap", nrow(df.final) )
  df.final$iter = NULL
  return(df.final )
}

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

gtex.ss250 = gtex.ss250[,c(10,11,1)]
archs.ss250 = archs.ss250[,c(10,11,1)]
gtex.ss250$rank = rep(1, nrow(gtex.ss250))
archs.ss250$rank = rep(1, nrow(archs.ss250))
gtex.ss250$size = rep(250, nrow(gtex.ss250))
archs.ss250$size = rep(250, nrow(archs.ss250))
names(gtex.ss250) = c("SID", "Tissue", "Gene", "Rank", "size")
names(archs.ss250) = c("SID", "Tissue", "Gene", "Rank", "size")


archs.gene250.sp = split(archs.ss250, f= archs.ss250$`Tissue` )
gtex.gene250.sp = split(gtex.ss250, f= gtex.ss250$`Tissue` )

gtex.df.250 = compute_jaccard(gtex.gene250.sp)
archs.df.250 = compute_jaccard(archs.gene250.sp)
gtex.df.250$type = rep("gtex", nrow(gtex.df.250))
archs.df.250$type = rep("archs", nrow(archs.df.250))

gtex.x.atchs = jc_gtex_X_archs(gtex.gene250.sp , archs.gene250.sp)
gtex.x.atchs$iter = NULL
gtex.x.atchs$type = rep("gtex.x.atchs", nrow(gtex.x.atchs))

#gtex.random.bootstrap = generate_random_Jc_distributions_per_tissue (gtex.gene250.sp, reps = 20, subset = 10)
df.final = rbind( gtex.df.250, archs.df.250, gtex.x.atchs )
#write.table(df.final, file = "Jc-gtex-archs-pairwise-comparisons.csv", sep = ",", row.names = F, col.names = T, quote = F)

cl = c(
  "#7e1e9c","#15b01a","#0343df","#ff81c0","#653700","#e50000","#95d0fc","#f97306","#029386","#c20078",
      "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
      "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00"
      )

gtex.n = names(gtex.gene250.sp)

archs.n = names(archs.gene250.sp)

shared = gtex.n[which(gtex.n %in% archs.n)]
df.final = df.final[df.final$V1 %in% shared, ]
ggplot(df.final, 
       aes(x = factor(type), 
           y = as.numeric(distances), fill = factor(type))) +
    geom_hline(yintercept = .3333 , color ="red", size = .75) + 
    facet_wrap( ~V1 , scales = "free_x", ncol = 8 ) + 
    geom_violin(aes(alpha = .5)) + 
    scale_fill_manual(values = cl) + 
    scale_x_discrete("Gene list size") + 
    scale_y_continuous("Jaccard Coeficient") + 
    theme_bw() + 
    theme(axis.title = element_text(size = 15), 
          panel.grid = element_line(size = 1))


## Figure 5B/S9 -------------------------------------------------------------
### Distributions of pairwise jaccard coefficients across GTEx and ARCHS for sGES 
#### compute pairwise jaccard coefficients

re.gene = fread("figures/data/reconstruction_errors/reconstruction.error.gene.archs.gtex.csv", 
            header = F , stringsAsFactors = F)
            
gtex.gene250 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.250.csv", header = F , stringsAsFactors = F)
archs.gene250 = fread("figures/data/archs/archs.tissue.genelist.formated", header = F , stringsAsFactors = F)
colnames(gtex.gene250) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(archs.gene250) =c("SID", "Tissue", "Gene", "Rank")


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

gtex.n.f = gtex.ss250.f$V11 %>% unique %>% as.character
archs.n.f = archs.ss250.f$V11 %>% unique %>% as.character
gtex.n.sf = gtex.ss250.sf$V11 %>% unique %>% as.character
archs.n.sf = archs.ss250.sf$V11 %>% unique %>% as.character
gtex.n.fm = gtex.ss250.fm$V11 %>% unique %>% as.character
archs.n.fm = archs.ss250.fm$V11 %>% unique %>% as.character
gtex.n.d = gtex.ss250.d$V11 %>% unique %>% as.character
archs.n.d = archs.ss250.d$V11 %>% unique %>% as.character



shared.f = gtex.n.f[which(gtex.n.f %in% archs.n.f)] 
shared.sf = gtex.n.sf[which(gtex.n.sf %in% archs.n.sf)] 
shared.fm = gtex.n.fm[which(gtex.n.fm %in% archs.n.fm)] 
shared.d = gtex.n.d[which(gtex.n.d %in% archs.n.d)] 

gtex.ss250.f = gtex.ss250.f[gtex.ss250.f$V11 %in% shared.f, ]
archs.ss250.f = archs.ss250.f[archs.ss250.f$V11 %in% shared.f, ]
gtex.ss250.sf = gtex.ss250.sf[gtex.ss250.sf$V11 %in% shared.sf, ]
archs.ss250.sf = archs.ss250.sf[archs.ss250.sf$V11 %in% shared.sf, ]
gtex.ss250.fm = gtex.ss250.fm[gtex.ss250.fm$V11 %in% shared.fm, ]
archs.ss250.fm = archs.ss250.fm[archs.ss250.fm$V11 %in% shared.fm, ]
gtex.ss250.d = gtex.ss250.d[gtex.ss250.d$V11 %in% shared.d, ]
archs.ss250.d = archs.ss250.d[archs.ss250.d$V11 %in% shared.d, ]

colnames(gtex.gene250) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(archs.gene250) =c("SID", "Tissue", "Gene", "Rank")

gtex.gene250$size = rep(250, nrow(gtex.gene250))
archs.gene250$size = rep(250, nrow(archs.gene250))

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

gtex.x.atchs.g = jc_gtex_X_archs(archs.gene250.sp  , gtex.gene250.sp)
gtex.x.atchs.f = jc_gtex_X_archs(archs.ss250.f.sp  , gtex.ss250.f.sp)
gtex.x.atchs.sf = jc_gtex_X_archs(archs.ss250.sf.sp  , gtex.ss250.sf.sp)
gtex.x.atchs.fm = jc_gtex_X_archs(archs.ss250.fm.sp  , gtex.ss250.fm.sp)
gtex.x.atchs.d = jc_gtex_X_archs(archs.ss250.d.sp  , gtex.ss250.d.sp)

gtex.x.atchs.g$type = rep("gene", nrow(gtex.x.atchs.g))
gtex.x.atchs.f$type = rep("fold", nrow(gtex.x.atchs.f))
gtex.x.atchs.fm$type = rep("family", nrow(gtex.x.atchs.fm))
gtex.x.atchs.sf$type = rep("superfamily", nrow(gtex.x.atchs.sf))
gtex.x.atchs.d$type = rep("domain", nrow(gtex.x.atchs.d))

df.final = rbind( gtex.x.atchs.g,
                  gtex.x.atchs.d, 
                  gtex.x.atchs.fm, 
                  gtex.x.atchs.sf, 
                  gtex.x.atchs.f )

gtex.n = names(gtex.gene250.sp)
archs.n = names(archs.gene250.sp)
shared = gtex.n[which(gtex.n %in% archs.n)]
df.final = df.final[df.final$V1 %in% shared, ]
df.final$distances = df.final$distances %>% as.character %>% as.numeric
sel.tiss = c("Spleen", "Ovary", "Muscle", "Pancreas", "Heart")
df.final.sel = df.final[df.final$V1 %in% sel.tiss,] 

cl = list( gene = "#51B7A1", domain = "#2274A5" , family = "#F75C03" , superfamily = "#00CB64" , fold = "#F5D34E"    )     

ggplot(df.final.sel, 
       aes(x = factor(type), 
           y = as.numeric(distances), fill = factor(type))) +
    geom_hline(yintercept = .25 , color ="red", size = .75) + 
    facet_wrap( ~V1 , scales = "free_x", ncol = 8 ) + 
    geom_violin(aes(alpha = .5)) + 
    scale_fill_manual(values = cl) + 
    scale_x_discrete("Signature Type", limits = c("gene", "domain", "family", "superfamily", "fold")) + 
    scale_y_continuous("Jaccard Coeficient") + 
    theme_bw() + 
    theme(axis.title = element_blank(), 
        axis.text = element_text(size = 15),
        legend.position = "none" , 
          panel.grid = element_line(size = 1))

##Figure S9
ggplot(df.final, 
       aes(x = factor(type), 
           y = as.numeric(distances), fill = factor(type))) +
    geom_hline(yintercept = .25 , color ="red", size = .75) + 
    facet_wrap( ~V1 , scales = "free_x", ncol = 8 ) + 
    geom_violin(aes(alpha = .5)) + 
    scale_fill_manual(values = cl) + 
    scale_x_discrete("Signature Type", limits = c("gene", "domain", "family", "superfamily", "fold")) + 
    scale_y_continuous("Jaccard Coeficient") + 
    theme_bw() + 
    theme(axis.title = element_blank(), 
        axis.text = element_text(size = 15),
        legend.position = "none" , 
          panel.grid = element_line(size = 1))











###Figure S14 
###ARchs consistancy with /wo outliers 
archs.outliers = read.table("archs.outliers.id.1sd.csv", sep =",")
archs.gene250$gsm = str_split_fixed(archs.gene250$SID, "-", 2)[,1]
archs.gene250.no.outliers = archs.gene250[! archs.gene250$gsm %in% archs.outliers$V1,]
archs.gene250.no.outliers$gsm = NULL 
archs.gene250$gsm = NULL
archs.gene250.sp = split(archs.gene250, f= archs.gene250$`Tissue` )
archs.gene250.no.outliers.sp = split(archs.gene250.no.outliers, f= archs.gene250.no.outliers$`Tissue` )
archs.df.250 = compute_jaccard(archs.gene250.sp)
archs.df.250.no.outliers = compute_jaccard(archs.gene250.no.outliers.sp )
archs.df.250$type = rep("ARCHS4 with outliers", nrow(archs.df.250))
gtex.gene250.sp = split(gtex.gene250, f= gtex.gene250$`Tissue` )

archs.df.250 = compute_jaccard(archs.gene250.sp)
archs.df.no.outlier.250 = compute_jaccard(archs.gene250.no.outliers.sp)
archs.df.no.outlier.250$type = rep("ARCHS4 without outliers", nrow(archs.df.no.outlier.250))
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

gtex.x.archs = jc_gtex_X_archs(gtex.gene250.sp , archs.gene250.sp)
gtex.x.archs.no.outliers = jc_gtex_X_archs(gtex.gene250.sp , archs.gene250.no.outliers.sp)
gtex.x.archs$type = rep("GTEX X ACRCHS4 with outliers", nrow(gtex.x.archs))
gtex.x.archs.no.outliers$type = rep("GTEX X ACRCHS4 without outliers", nrow(gtex.x.archs.no.outliers))
gtex.x.archs$iter = NULL
gtex.x.archs.no.outliers$iter = NULL


###Figure 5E 
##consistancy of overlap between archs and gtex before and after outlier removal 

archs.outliers = read.table("archs.outliers.id.2sd.csv", sep =",")

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



###Table 1 and 2 
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

df.signatures = rbind(signature.gene,signature.domain,  signature.fm,  signature.spr, signature.fold )
write.table(df.signatures, "tissue-specific-signatures-GES-sGES.csv",sep =",", quote =F, , row.names =F)