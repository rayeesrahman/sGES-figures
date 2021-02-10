rm(list = ls())
setwd("sGES-figures/")
library(ggplot2)
library(magrittr)
library(vegan)
library(tidyverse)
library(data.table)
## Figure 3A -------------------------------------------------------------
### Jc null distribution 
##Read data 

gene50 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.50.csv", header = F , stringsAsFactors = F)
gene250 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.250.csv", header = F , stringsAsFactors = F)
gene1000 = fread("figures/data/gtex/gtex.ranked.genelist.top.gene.1000.csv", header = T , stringsAsFactors = F)
colnames(gene50) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(gene250) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")
colnames(gene1000) =c("SID", "Sub-Tissue", "Tissue", "Gene", "Rank")

gene50$size = rep(50, nrow(gene50))
gene250$size = rep(250, nrow(gene250))
gene1000$size = rep(1000, nrow(gene1000))

gene50.sp = split(gene50, f= gene50$`Sub-Tissue`)
gene250.sp = split(gene250, f= gene250$`Sub-Tissue` )
gene1000.sp = split(gene1000, f= gene1000$`Sub-Tissue` )
random_comparisions = function(splitfile, bootstraps )
{
  df.final = data.frame() 
  for ( i in 1:bootstraps) 
  {
    num.samp = 1:length(splitfile)
    print(paste0("bootstrap number ", i ))
    tisspos1 = sample(num.samp,1)
    num.samp = num.samp[-tisspos1]
    tisspos2 = sample(num.samp,1)
    tiss.1 = splitfile[[tisspos1]]
    tiss.2 = splitfile[[tisspos2]]
    tiss.1 = split(tiss.1, tiss.1$SID) 
    tiss.2 = split(tiss.2, tiss.2$SID) 
    tiss.combined = rbind(tiss.1[[sample(length(tiss.1),1)]], tiss.2[[sample(length(tiss.2),1)]]) 
    set = tiss.combined$size %>%  as.numeric()  %>% unique()
    tiss.combined[,c(2,3,5,6)] = NULL
    tiss.combined$val = rep(1, nrow(tiss.combined))
    tiss.wide = spread(data = tiss.combined, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
    row.names(tiss.wide) = tiss.wide$SID 
    tiss.wide$SID = NULL
    tiss.dist = vegdist(data.matrix(tiss.wide), method = "jaccard", binary = T ,
                     upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
    distances = tiss.dist[upper.tri(tiss.dist)] %>%  as.numeric() 
    distances = 1 - distances 
    df = data.frame()
    df = cbind(rep(i, length(distances)), distances, set) %>%  as.data.frame()
    colnames(df) = c("bootstrap", "distances", "set")
    df.final = rbind(df.final, df)
  }
  return(df.final)
}

gene50.random = random_comparisions(gene50.sp, 1000)
gene250.random = random_comparisions(gene250.sp, 1000)
gene1000.random = random_comparisions(gene1000.sp, 1000)
random.all = rbind(gene50.random,gene250.random,gene1000.random )


##plot
ggplot(random.all, 
       aes(x = as.factor(random.all$set), 
           y = distances)) +
  geom_violin(aes(alpha = .5), fill = "grey") + 
  scale_fill_manual() + 
  scale_x_discrete("Gene list size") + 
  scale_y_continuous("Jaccard Coeficient", ,limit = c(0,1)) + 
  geom_boxplot(outlier.shape = NA, color = "black", fill = "grey") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 15), 
        panel.grid = element_line(size = 1), 
        axis.text = element_text(size = 30))
          

benchmark = median(gene250.random$distances)


## Figure 3B -------------------------------------------------------------
### Robustness of structural signatures 
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

compute_jaccard  = function(df.sp , type ){
    df.final = data.frame()
    for ( i in df.sp)
    {
        tiss = i$`Sub-Tissue` %>% as.character() %>%  unique  
        overalltissue = i$Tissue %>%  as.character() %>% unique
        print(tiss)
        i[,c(2,3,5)] = NULL
        i$val = rep( 1, nrow(i))
        i.wide = spread(data = i, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
        row.names(i.wide) = i.wide$SID 
        i.wide$SID = NULL
        i.dist = vegdist(data.matrix(i.wide), method = "jaccard", binary = T ,
                         upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
        distances = i.dist[upper.tri(i.dist)] %>%  as.numeric() 
        distances = 1 - distances 
        df = data.frame()
        df = cbind(rep(tiss, length(distances)), distances) %>%  as.data.frame()
        df$type = rep( type , nrow(df))
        df$tissue = rep(overalltissue , nrow(df))
        df.final = rbind(df.final, df)
    }
    df.final$V1 = df.final$V1 %>%  as.factor()
    df.final$distances = df.final$distances %>%  as.character() %>% as.numeric()
    return(df.final)
} 

gene.dist = compute_jaccard(genes.sp, "genes")
domain.dist = compute_jaccard(domain.sp, "domain")
fold.dist = compute_jaccard(fold.sp, "fold")
family.dist = compute_jaccard(family.sp, "family")
sfam.dist = compute_jaccard(sfam.sp, "super")

combined = rbind(gene.dist ,  domain.dist, family.dist,  sfam.dist, fold.dist  ) %>%  as.data.frame()


cl = list( genes = "#51B7A1", domain = "#2274A5" , family = "#F75C03" , super = "#00CB64" , fold = "#F5D34E"    )       

sel.tiss = c("Pancreas",  "Cervix", "Ovary", "Heart", "Spleen", "Muscle")
combined.sel = combined[combined$tissue %in% sel.tiss,] 
##Selected 
###Figure 3B
pl =  ggplot(combined.sel, aes(factor(combined.sel$type), combined.sel$distances, fill = combined.sel$type)) + 
    facet_wrap(~V1,  ncol = 4 ) + 
    geom_violin(aes(alpha = .5)) +
    geom_boxplot() +
    scale_fill_manual(values = cl) + 
    geom_hline(yintercept = benchmark, color ="red", size = .75)  + 
    scale_x_discrete( limits = c("genes", "domain", "family", "super", "fold")) + 
    #scale_y_continuous("Jaccard Coeficient") +
    geom_boxplot(outlier.shape = NA, color = "black") + 
    theme_bw() + 
    theme(axis.title = element_blank(), 
          panel.grid = element_line(size = 1),
          axis.text = element_text(size = 20), 
          legend.position ="none")
##Full 
###Figure S4
ggplot(combined, aes(factor(combined$type), combined$distances, fill = combined$type)) + 
    facet_wrap(~V1,  ncol = 4 ) + 
    geom_violin(aes(alpha = .5)) +
    geom_boxplot() +
    scale_fill_manual(values = cl) + 
    geom_hline(yintercept = .25 , color ="red", size = .75)  + 
    scale_x_discrete( limits = c("genes", "domain", "family", "super", "fold")) + 
    #scale_y_continuous("Jaccard Coeficient") +
    geom_boxplot(outlier.shape = NA, color = "black") + 
    theme_bw() + 
    theme(axis.title = element_blank(), 
          panel.grid = element_line(size = 1),
          axis.text.x = element_blank())



####Figure 3C
### generate distrubutions of control comparisons, such as across tissues, and randomly selecting structures for sGES
random_comparisions = function(splitfile, bootstraps, type  )
{
    df.final = data.frame() 
    for ( i in 1:bootstraps)  
    {
        num.samp = 1:length(splitfile)
        print(paste0("bootstrap number ", i ))
        tisspos1 = sample(num.samp,1)
        num.samp = num.samp[-tisspos1]
        tisspos2 = sample(num.samp,1)
        tiss.1 = splitfile[[tisspos1]]
        tiss.2 = splitfile[[tisspos2]]
        tiss.1$SID = factor(tiss.1$SID)
        tiss.2$SID = factor(tiss.2$SID)
        tiss.1 = split(tiss.1, tiss.1$SID) 
        tiss.2 = split(tiss.2, tiss.2$SID) 
        tiss.combined = rbind(tiss.1[[sample(length(tiss.1),1)]], tiss.2[[sample(length(tiss.2),1)]]) 
        tiss.combined[,c(2,3,5,6)] = NULL
        tiss.combined$val = rep(1, nrow(tiss.combined))
        tiss.wide = spread(data = tiss.combined, key = Gene ,  value = val,  fill = 0 ) %>%  as.data.frame()    
        row.names(tiss.wide) = tiss.wide$SID 
        tiss.wide$SID = NULL
        tiss.dist = vegdist(data.matrix(tiss.wide), method = "jaccard", binary = T ,
                            upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
        distances = tiss.dist[upper.tri(tiss.dist)] %>%  as.numeric() 
        distances = 1 - distances 
        df = data.frame()
        df = cbind(rep(i, length(distances)), distances, type) %>%  as.data.frame()
        colnames(df) = c("bootstrap", "distances", "type")
        df.final = rbind(df.final, df)
    }
    return(df.final)
}
#gene.random = random_comparisions(genes.sp, 1000, "genes") 
fold.random = random_comparisions(fold.sp, 1000, "fold") 
family.random = random_comparisions(family.sp, 1000, "family") 
sfam.random = random_comparisions(sfam.sp, 1000, "super") 
domain.random = random_comparisions(domain.sp, 1000, "domain") 
random.combined = rbind(domain.random,sfam.random,family.random,fold.random)
random.combined$group = rep("across-tissues", nrow(random.combined))
combined.2 = combined[,c(1,2,3)]
names(combined.2)  = c("bootstrap", "distances", "type")
combined.2$group = rep("within-tissues", nrow(combined.2))
within.across.combined = rbind(combined.2, random.combined)

### take counts from dataframe and repeat structure n times, convert into vector
domain.cnt = fread("files/backgrounds/human_proteome.background.ipr.domain.cnt", header = F, sep =",")
fold.cnt = fread("files/backgrounds/human_proteome.background.scop.fold.cnt", header = F, sep =",")
family.cnt = fread("files/backgrounds/human_proteome.background.scop.family.cnt", header = F, sep =",")
sfam.cnt = fread("files/backgrounds/human_proteome.background.scop.superfam.cnt", header = F, sep =",")
domain.sampl = apply(domain.cnt, 1, 
                    function(x)
                        {
                            z = rep(as.character(x[1]),as.numeric(x[2]) )
                            z
                        } ) %>% unlist 
fold.sampl = apply(fold.cnt, 1, 
                    function(x)
                        {
                            z = rep(as.character(x[1]),as.numeric(x[2]) )
                            z
                        } ) %>% unlist 
fam.sampl = apply(family.cnt, 1, 
                    function(x)
                        {
                            z = rep(as.character(x[1]),as.numeric(x[2]) )
                            z
                        } ) %>% unlist 
sfam.sampl = apply(sfam.cnt, 1, 
                    function(x)
                        {
                            z = rep(as.character(x[1]),as.numeric(x[2]) )
                            z
                        } ) %>% unlist 
generate_jaccard_dist = function(sampl, bootstraps, nstruct , type  )
{
    sampl = sampl[sample(1:length(sampl), length(sampl) )] ##randomly shuffle input sample 
    distances = c()
    for (i in 1:bootstraps)
    {
        print(paste0("bootstrap number: ", i))
        data.1 = sampl[sample(1:length(sampl), nstruct )] %>% unique
        data.2 = sampl[sample(1:length(sampl), nstruct )] %>% unique 
        i = intersect(data.1, data.2) %>% length
        u = union(data.1, data.2) %>% length
        data.1 = cbind(rep("sample1", length(data.1)), data.1)   
        data.2 = cbind(rep("sample2", length(data.2)), data.2)   
        data.combined = rbind(data.1, data.2) %>% as.data.table
        names(data.combined) = c("sample", "structure")
        data.combined$presense = rep(1, nrow(data.combined))
        data.wide = spread(data.combined, key = structure, value = presense , fill = 0)
        data.wide$sample = NULL
        d =  vegdist(data.wide, method="jaccard" ,  binary = T ) %>% as.vector %>% as.numeric
        d = 1 - d 
        distances = c(d, distances)
    }
    dist.1 = cbind(distances, rep(type, length(distances)))
    dist.2 = cbind(1:length(distances), dist.1) %>% as.data.frame
    dist.2$group = rep("random_structures", length(distances))
    names(dist.2) = c("bootstrap", "distances", "type", "group")
    return(dist.2)
}
domain.dist = generate_jaccard_dist(domain.sampl, 1000, 250, "domain" )
fam.dist = generate_jaccard_dist(fam.sampl, 1000, 150, "family"  )
sfam.dist = generate_jaccard_dist(sfam.sampl, 1000, 150, "super")
fold.dist = generate_jaccard_dist(fold.sampl, 1000, 150 , "fold")
random.struct.combined = rbind(domain.dist , fam.dist, sfam.dist, fold.dist)
within.across.struct.combined = rbind(within.across.combined, random.combined, domain.dist , fam.dist, sfam.dist, fold.dist)
   
within.across.struct.combined = within.across.struct.combined[ which(within.across.struct.combined$type != "genes"  ) ,] 

ggplot(within.across.struct.combined, 
        aes(x = factor(within.across.struct.combined$group), 
            y = as.numeric(as.character(within.across.struct.combined$distances)), 
            fill = factor(within.across.struct.combined$type)))  + 
    facet_wrap(~type,  ncol = 2 ) + 
    geom_violin(aes(alpha = .5)) +
    geom_boxplot() +
    scale_fill_manual(values = cl) + 
    scale_x_discrete("", 
        limits = c("random_structures", "across-tissues", "within-tissues")) + 
    theme_bw() + 
    geom_hline(yintercept = benchmark , color ="red", size = .75)  + 
    scale_y_continuous("Jaccard Coeficient") +
    geom_boxplot(outlier.shape = NA, color = "black") + 
    theme(axis.title = element_blank(), 
          panel.grid = element_line(size = 1), 
          axis.text = element_text(size = 20), 
          legend.position = "none")


## Figure 3D -------------------------------------------------------------
### Benchmarking of sGES on other datasets
#### read data 

###read data 
#########select tissue types from below 
### [1] "esophagus"      "heart"          "kidney"         "liver"         
### [5] "pancreas"       "smallintestine" "stomach 

tiss = list(Esophagus = "esophagus", Kidney = "kidney", Liver = "liver", Pancreas = "pancreas",  Stomach = "stomach")


compute_jaccard_ss = function(msigdb.domain.sp, domain.sp, tissue , tissue2 , method , type )
{
    df.final = data.frame()
    msigdb.domain.sel  = msigdb.domain.sp[[tissue2]]
    msigdb.sel.sp = split(msigdb.domain.sel , f = msigdb.domain.sel $SID)
    total = length(msigdb.sel.sp)
    cnt = 0 
    for ( i in msigdb.sel.sp )
    {
        cnt = cnt + 1
        print(paste("working on" , cnt , "of", total , sep = " "))
        gtex.tiss = domain.sp[[tissue]]
        gtex.tiss.sp =   split(gtex.tiss, f = gtex.tiss$SID)
        #i$method = NULL
        i$Tissue = NULL 
        i$Value = 1 
        cmpname = i$SID %>% unique 
        for ( ii in gtex.tiss.sp )
        {
            ii$`Sub-Tissue` = NULL 
            ii$`Tissue` = NULL 
            ii$Value =1 
            gtex.sid = ii$SID %>% unique
            #print(paste(cmpname, gtex.sid, sep=  "|")) 
            cmp = rbind(ii, i)       
            i.wide = spread(data = cmp, key = Gene ,  value = Value,  fill = 0 ) %>%  as.data.frame()    
            row.names(i.wide) = i.wide$SID 
            i.wide$SID = NULL
            i.dist = vegdist(data.matrix(i.wide), method = "jaccard", binary = T ,
                            upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
            distances = i.dist[upper.tri(i.dist)] %>%  as.numeric() 
            distances = 1 - distances 
            df = data.frame(method = method, tissue = tissue, distance = distances , type = type, compare = paste(cmpname, gtex.sid, sep=  "|"))
            df.final = rbind(df.final, df)
        }
    }
    return(df.final)
}


compute_jaccard_genes = function(msigdb.sp, genes.sp, tissue , method , type = "Gene")
{
    df.final = data.frame()
    msigdb.sp.sel  = msigdb.sp[[tissue2]]
    msigdb.gl.sp = split(msigdb.sp.sel , f = msigdb.sp.sel$SID)
    total = length(msigdb.gl.sp)
    cnt = 0 
    for ( i in msigdb.gl.sp )
    {
        cnt = cnt + 1
        print(paste("working on" , cnt , "of", total , sep = " "))
        ngenes = nrow(i)
        gtex.tiss = genes.sp[[tissue]]
        gtex.tiss.sp =   split(gtex.tiss, f = gtex.tiss$SID)
        i$method = NULL
        i$Tissue = NULL 
        i$Value = 1 
        cmpname = i$SID %>% unique 
        for ( ii in gtex.tiss.sp )
        {
            gtex.topn = head(ii, ngenes)
            gtex.topn$`Sub-Tissue` = NULL 
            gtex.topn$`Tissue` = NULL 
            gtex.topn$`Rank` = NULL 
            gtex.topn$Value =1 
            gtex.sid = gtex.topn$SID %>% unique
            #print(paste(cmpname, gtex.sid, sep=  "|")) 
            cmp = rbind(gtex.topn, i)       
            i.wide = spread(data = cmp, key = Gene ,  value = Value,  fill = 0 ) %>%  as.data.frame()    
            row.names(i.wide) = i.wide$SID 
            i.wide$SID = NULL
            i.dist = vegdist(data.matrix(i.wide), method = "jaccard", binary = T ,
                            upper = T, diag = T ) %>% data.matrix() %>%  as.data.frame() 
            distances = i.dist[upper.tri(i.dist)] %>%  as.numeric() 
            distances = 1 - distances 
            df = data.frame(method = method, tissue = tissue, distance = distances , type = type , compare = paste(cmpname, gtex.sid, sep=  "|"))
            df.final = rbind(df.final, df)
        }
    }
    return(df.final)
}



###read data 
genes = fread("gtex.ranked.genelist.top.gene.1000.csv", header = F , stringsAsFactors = F, data.table = F )
domain = fread("allcombined.gtex.250.domain.csv", sep = "," , stringsAsFactors = F,  header = F, data.table = F )
fold = fread("allcombined.gtex.250.fold.csv", sep = "," , stringsAsFactors = F,  header = F, data.table = F )
family = fread("allcombined.gtex.250.family.csv", sep = "," , stringsAsFactors = F,  header = F, data.table = F )
sfam = fread("allcombined.gtex.250.superfam.csv", sep = "," , stringsAsFactors = F, header = F, data.table = F  )
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

genes.sp = split(x = genes, f = genes$`Tissue`)
domain.sp = split(x = domain, f = domain$`Tissue`)
fold.sp = split(x = fold, f = fold$`Tissue`)
family.sp = split(x = family, f = family$`Tissue`)
sfam.sp = split(x = sfam, f = sfam$`Tissue`)


tissue = "Heart" 
tissue2 = "heart"  
print(paste("working on" , tissue, sep = " "))
msigdb.gl = fread("figures/data/msigdb/allcombined.msigdb.genes.2.csv", sep = ",", header = F)
names(msigdb.gl) = c("method", "Tissue", "SID", "Gene" )
msigdb.gl.sp = split(msigdb.gl, f = msigdb.gl$Tissue) 
method = "msigdb"

msigdb.heart.gene.jc = compute_jaccard_genes(msigdb.gl.sp , genes.sp , tissue , method )

### structures domain

msigdb.domain = fread("figures/data/msigdb/allcombined.msigdb.domain.2.csv", sep = ",", header = F)
msigdb.domain = msigdb.domain[,c(1,10,11)]
names(msigdb.domain) = c("Gene", "SID", "Tissue"  )

msigdb.domain.sp = split(msigdb.domain, f = msigdb.domain$Tissue)



msigdb.heart.domain.jc  = compute_jaccard_ss(msigdb.domain.sp, domain.sp, tissue , tissue2 , method , type = "Domain" )
head(msigdb.heart.domain.jc)

###family 
msigdb.domain = fread("figures/data/msigdb/allcombined.msigdb.family.2.csv", sep = ",", header = F)
msigdb.domain = msigdb.domain[,c(1,10,11)]
names(msigdb.domain) = c("Gene", "SID", "Tissue"  )
msigdb.domain.sp = split(msigdb.domain, f = msigdb.domain$Tissue)
msigdb.heart.family.jc  = compute_jaccard_ss(msigdb.domain.sp, family.sp, tissue , tissue2 , method , type = "Family" )
head(msigdb.heart.family.jc)

###family 
msigdb.domain = fread("figures/data/msigdb/allcombined.msigdb.superfamily.2.csv", sep = ",", header = F)
msigdb.domain = msigdb.domain[,c(1,10,11)]
names(msigdb.domain) = c("Gene", "SID", "Tissue"  )
msigdb.domain.sp = split(msigdb.domain, f = msigdb.domain$Tissue)
msigdb.heart.superfamily.jc  = compute_jaccard_ss(msigdb.domain.sp, sfam.sp, tissue , tissue2 , method , type = "Superfamily" )
head(msigdb.heart.superfamily.jc)

###fold 
msigdb.domain = fread("figures/data/msigdb/allcombined.msigdb.fold.2.csv", sep = ",", header = F)
msigdb.domain = msigdb.domain[,c(1,10,11)]
names(msigdb.domain) = c("Gene", "SID", "Tissue"  )
msigdb.domain.sp = split(msigdb.domain, f = msigdb.domain$Tissue)
msigdb.heart.fold.jc  = compute_jaccard_ss(msigdb.domain.sp, fold.sp, tissue , tissue2 , method , type = "fold" )
head(msigdb.heart.fold.jc )


###human protein atlas 
####analyze genes 

hpa.gl = fread("figures/data/humanprotatlas/allcombined.hpa.2.gl", sep = ",", header = F)
names(hpa.gl) = c("method", "Tissue", "SID", "Gene" )
hpa.gl.sp = split(hpa.gl, f = hpa.gl$Tissue) 
method = "hpa"

hpa.heart.gene.jc = compute_jaccard_genes(hpa.gl.sp , genes.sp , tissue , method )
head(hpa.heart.gene.jc )
### structures domain

hpa.domain = fread("figures/data/humanprotatlas/results/allcombined.hpa.domain.csv", sep = ",", header = F)
hpa.domain = hpa.domain[,c(1,10,11)]
names(hpa.domain) = c("Gene", "SID", "Tissue"  )
hpa.domain.sp = split(hpa.domain, f = hpa.domain$Tissue)
type = "Domain"
hpa.heart.domain.jc  = compute_jaccard_ss(hpa.domain.sp, domain.sp, tissue , tissue2 , method , type )
head(hpa.heart.domain.jc )

###family 
hpa.domain = fread("figures/data/humanprotatlas/results/allcombined.hpa.family.csv", sep = ",", header = F)
hpa.domain = hpa.domain[,c(1,10,11)]
names(hpa.domain) = c("Gene", "SID", "Tissue"  )
hpa.domain.sp = split(hpa.domain, f = hpa.domain$Tissue)
hpa.heart.family.jc  = compute_jaccard_ss(hpa.domain.sp, family.sp, tissue , tissue2 , method , type = "Family" )
head(hpa.heart.family.jc)
###superfamily 
hpa.domain = fread("figures/data/humanprotatlas/results/allcombined.hpa.superfamily.csv", sep = ",", header = F)
hpa.domain = hpa.domain[,c(1,10,11)]
names(hpa.domain) = c("Gene", "SID", "Tissue"  )
hpa.domain.sp = split(hpa.domain, f = hpa.domain$Tissue)
hpa.heart.superfamily.jc  = compute_jaccard_ss(hpa.domain.sp, sfam.sp, tissue , tissue2 , method , type = "Superfamily" )
head(hpa.heart.superfamily.jc)
###fold 
hpa.domain = fread("figures/data/humanprotatlas/results/allcombined.hpa.fold.csv", sep = ",", header = F)
hpa.domain = hpa.domain[,c(1,10,11)]
names(hpa.domain) = c("Gene", "SID", "Tissue"  )
hpa.domain.sp = split(hpa.domain, f = hpa.domain$Tissue)
hpa.heart.fold.jc  = compute_jaccard_ss(hpa.domain.sp, fold.sp, tissue , tissue2 , method , type = "fold" )
head(hpa.heart.fold.jc )


##############TISSUES 


tissues.gl = fread("figures/data/TISSUES-harmonizome/allcombined.tissues.inte.genes.2.csv", sep = ",", header = F)
names(tissues.gl) = c("method", "Tissue", "SID", "Gene" )
tissues.gl.sp = split(tissues.gl, f = tissues.gl$Tissue) 
method = "tissues"
type = "Domain"

tissues.heart.gene.jc = compute_jaccard_genes(tissues.gl.sp , genes.sp , tissue , method )

### structures domain

tissues.domain = fread("figures/data/TISSUES-harmonizome/allcombined.tissue-inte.domain.2.csv", sep = ",", header = F)
tissues.domain = tissues.domain[,c(1,10,11)]
names(tissues.domain) = c("Gene", "SID", "Tissue"  )

tissues.domain.sp = split(tissues.domain, f = tissues.domain$Tissue)



tissues.heart.domain.jc  = compute_jaccard_ss(tissues.domain.sp, domain.sp, tissue , tissue2 , method , type )
head(tissues.heart.domain.jc)

###family 
tissues.domain = fread("figures/data/TISSUES-harmonizome/allcombined.tissue-inte.family.2.csv", sep = ",", header = F)
tissues.domain = tissues.domain[,c(1,10,11)]
names(tissues.domain) = c("Gene", "SID", "Tissue"  )
tissues.domain.sp = split(tissues.domain, f = tissues.domain$Tissue)
tissues.heart.family.jc  = compute_jaccard_ss(tissues.domain.sp, family.sp, tissue , tissue2 , method , type = "Family" )
head(tissues.heart.family.jc)

###family 
tissues.domain = fread("figures/data/TISSUES-harmonizome/allcombined.tissue-inte.superfam.2.csv", sep = ",", header = F)
tissues.domain = tissues.domain[,c(1,10,11)]
names(tissues.domain) = c("Gene", "SID", "Tissue"  )
tissues.domain.sp = split(tissues.domain, f = tissues.domain$Tissue)
tissues.heart.superfamily.jc  = compute_jaccard_ss(tissues.domain.sp, sfam.sp, tissue , tissue2 , method , type = "Superfamily" )
head(tissues.heart.superfamily.jc)

###fold 
tissues.domain = fread("figures/data/TISSUES-harmonizome/allcombined.tissue-inte.fold.2.csv", sep = ",", header = F)
tissues.domain = tissues.domain[,c(1,10,11)]
names(tissues.domain) = c("Gene", "SID", "Tissue"  )
tissues.domain.sp = split(tissues.domain, f = tissues.domain$Tissue)
tissues.heart.fold.jc  = compute_jaccard_ss(tissues.domain.sp, fold.sp, tissue , tissue2 , method , type = "fold" )
head(tissues.heart.fold.jc )


###combine 
complete = rbind(tissues.heart.superfamily.jc, tissues.heart.family.jc , tissues.heart.domain.jc , tissues.heart.fold.jc ,tissues.heart.gene.jc , msigdb.heart.gene.jc, msigdb.heart.domain.jc, msigdb.heart.family.jc, msigdb.heart.superfamily.jc, msigdb.heart.fold.jc,  
                hpa.heart.gene.jc, hpa.heart.domain.jc, hpa.heart.family.jc, hpa.heart.superfamily.jc, hpa.heart.fold.jc) 
complete$method2 =  str_replace(complete$method, "msigdb", "mSigDB")
complete$method2 = str_replace(complete$method2, "tissues", "TISSUES2")
complete$method2 = str_replace(complete$method2, "hpa", "Human Protein Atlas")

cl = list( Gene = "#51B7A1", Domain = "#2274A5" , Family = "#F75C03" , Superfamily = "#00CB64" , fold = "#F5D34E"    )       

ggplot(complete, aes(type, distance, fill = type )) + 
    geom_boxplot() + 
    geom_violin(aes(alpha = .5)) +
    scale_x_discrete( limits = c("Gene", "Domain", "Family", "Superfamily", "fold")) + 
    geom_hline(yintercept = .333 , color ="red", size = .75)  + 
    facet_wrap(~method2) +
    geom_boxplot(outlier.shape = NA, color = "black") + 
    theme_bw() + 
    scale_fill_manual(values = cl) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(axis.title = element_blank(), 
        panel.grid = element_line(size = 1),
        axis.text = element_text(size = 25, colour="black", family = "Arial"), 
        strip.text = element_text(size = 40, colour="black", family = "Arial"),
        strip.background = element_rect(fill="white", colour="black",size=0.5),
        legend.position ="none")