rm(list = ls())
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
convert_to_wide = function(x, value_feat, key_feat)
{
    x = x[,c(key_feat, value_feat, "SID")]
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
    "SID"
)

cl = c(
    "#7e1e9c","#15b01a","#0343df","#ff81c0","#653700", "#e50000","#95d0fc","#f97306","#029386","#c20078",
       "#53fca1","#c04e01","#3f9b0b","#dbb40c","#580f41","#b9a281","#ff474c","#fffe7a","#40a368","#0a888a",
       "#887191","#be6400","#82cafc","#1fa774","#8cffdb","#7bb274","#510ac9","#ff5b00"
    )

domain = fread( "figures/data/dtox/top250_down_fold_combined.csv", sep = ",", header = F)
names(domain) = header 

domain.2 = domain[,c("Structure", "logfoldchange", "SID")]
x.wide = spread(data = domain.2, key = "Structure" ,  value = "logfoldchange",  fill = 0 ) %>%  as.data.frame()  
x.wide.meta = str_split_fixed(x.wide$SID, '\\.', 5) %>% as.data.frame
x.wide$cell = x.wide.meta$V2
x.wide$plate = x.wide.meta$V4
x.wide$nib = x.wide.meta$V5

library(iheatmapr)

selecteddrugs = c(
    "DAB",
    "NIL",
    "PAZ",
    "REG",
    "RUX",
    "SOR",
    "TRA",
    "VEM",
    "TGF",
    "DAU",
    "DOX",
    "EPI",
    "IDA"
)
drug_cate = list(
    "TRS" = "Protein kinase inhibitors", 
    "CTX" = "Protein kinase inhibitors",
    "SOR" = "Protein kinase inhibitors",
    "SUN" = "Protein kinase inhibitors",
    "IMA" = "Protein kinase inhibitors",
    "ERL" = "Protein kinase inhibitors", 
    "DAS" = "Protein kinase inhibitors",
    "TOF" = "Protein kinase inhibitors",
    "CRI" = "Protein kinase inhibitors",
    "GEF" = "Protein kinase inhibitors",
    "BOS" = "Protein kinase inhibitors",
    "VAN" = "Protein kinase inhibitors",
    "LAP" = "Protein kinase inhibitors",
    "NIL" = "Protein kinase inhibitors",
    "AXI" = "Protein kinase inhibitors",
    "PAZ" = "Protein kinase inhibitors",
    "RUX" = "Protein kinase inhibitors",
    "AFA" = "Protein kinase inhibitors",
    "REG" = "Protein kinase inhibitors",
    "PON" = "Protein kinase inhibitors",
    "DAB" = "Protein kinase inhibitors",
    "VEM" = "Protein kinase inhibitors",
    "CAB" = "Protein kinase inhibitors",
    "TRA" = "Protein kinase inhibitors",
    "CER" = "Protein kinase inhibitors",
    "EPI" = "Anthracyclines",
    "IDA" = "Anthracyclines",
    "DOX" = "Anthracyclines",
    "DAU" = "Anthracyclines",
    "AMI" = "Antiarrhythmics",
    "FLE" = "Antiarrhythmics",
    "PHP" = "CARDIAC STIMULANTS EXCL. CARDIAC GLYCOSIDES", 
    "ISO" = "CARDIAC STIMULANTS EXCL. CARDIAC GLYCOSIDES",
    "PRE" = "GENITO URINARY SYSTEM AND SEX HORMONES", 
    "TNF" = "Tumor necrosis factor alpha inhibitors" , 
    "MIL" = "CARDIAC STIMULANTS EXCL. CARDIAC GLYCOSIDES", 
    "EST" = "GENITO URINARY SYSTEM AND SEX HORMONES", 
    "IGF" = "GENITO URINARY SYSTEM AND SEX HORMONES", 
    "END" = "GENITO URINARY SYSTEM AND SEX HORMONES" , 
    "OLM" = "Angiotensin II receptor blockers (ARBs)" , 
    "DOB" = "CARDIAC STIMULANTS EXCL. CARDIAC GLYCOSIDES", 
    "VRP" = "Phenylalkylamine derivatives", 
    "ALE" =  "Bisphosphonates", 
    "DOM" =  "ALIMENTARY TRACT AND METABOLISM", 
    "DTP" =  "ALIMENTARY TRACT AND METABOLISM", 
    "ENT" =  "Nucleoside and nucleotide reverse transcriptase inhibitors", 
    "LOP" =  "ALIMENTARY TRACT AND METABOLISM", 
    "TGF" =  "Transforming growth factor inhibitors", 
    "CYC" =  "Calcineurin inhibitors", 
    "URS" =  "ALIMENTARY TRACT AND METABOLISM", 
    "PAR" =  "Selective serotonin reuptake inhibitors"
)

x.wide= x.wide[which(x.wide$nib %in% selecteddrugs) , ]
x.wide$nib = factor(x.wide$nib)
maxcol = ncol(x.wide) - 3
x.wide.2 = x.wide[,which(apply(x.wide[,2:maxcol], 2, function(x) {sum(x)}) %>% as.numeric() !=  0 ) + 1 ]
x.wide.dat = x.wide.2[,which(apply(x.wide.2[,2:ncol(x.wide.2)], 2, function(x) {min(x)}) %>% as.numeric() >=  0 ) + 1 ]
maxcol2 = ncol(x.wide.dat)
x.wide.dat$nib = x.wide$nib 
x.wide.dat$plate = x.wide$plate
x.wide.dat$cell = x.wide$cell 
x.wide.dat$drug_cat = as.factor(as.character(unlist(drug_cate[as.character(x.wide.dat$nib)]  ) ))
main_heatmap(data.matrix(x.wide.dat[,2:maxcol2]), name = "Log Fold Change") %>%
  add_col_clustering() %>%
  add_row_clustering() %>%
  add_row_title("Folds") %>%
  add_col_title("Sample") %>%
  add_row_annotation(data.frame(
      "Drug Category" = x.wide.dat$drug_cat, 
      "Drug" =  x.wide.dat$nib) )  #%>%
  add_row_labels(x.wide.dat$nib) #%>%
  add_row_annotation() #%>%
  add_row_annotation(x.wide.dat$cell) #%>%
  
  
  add_main_heatmap(Indometh_matrix,
                   name = "Indometacin<br>Concentration") %>%
  
  add_col_title("Time") %>%
  add_col_summary()



x.tsne  = Rtsne(x.wide, dims = 3, perplexity = 15, verbose = T )
x.tsne.y = x.tsne$Y %>% as.data.frame
x.tsne.y = cbind(x.tsne.y , x.wide.meta)

names(x.tsne.y) = c("TSNE1", "TSNE2", "TSNE3", "DRUG", "CELL", "ID", "TYPE")
plot_ly(x.tsne.y, 
        x = as.numeric(x.tsne.y$TSNE1),
        y = as.numeric(x.tsne.y$TSNE2), 
        z = as.numeric(x.tsne.y$TSNE3), 
        color = factor(x.tsne.y$DRUG) , 
        colors = cl ) %>% 
        layout(
            scene = list(
                xaxis = list(title = "TSNE1"),
                yaxis = list(title = "TSNE2"),
                zaxis = list(title = "TSNE3")
            )
        ) %>%  
        add_trace(
            text = paste( "cell:", x.tsne.y$CELL, "\ntype:", x.tsne.y$TYPE,  sep = " "),
            hoverinfo = 'text',
            showlegend = T
        )


######################################################################################
######################################################################################
## ## ## ## ## ## ## ## ## ## ## ## ## 


selecteddrugs = c(
    "regorafenib",
    "sorafenib",
    "Doxorubicin",
    "epirubicin",
    "idarubicin"
)


fold.down = fold.down[which(fold.down$drug %in% selecteddrugs) , ]
#fold.up = fold.up[pvalue< .05]
#fold.down = fold.down[pvalue< .05]
library(gplots)
fold.up = fold.up[,c("SID", "structure", "pvalue", "drug")] 
fold.down = fold.down[,c("SID", "structure", "log_fold_change", "drug")] 
fold.u.w = spread(fold.up, structure, pvalue, fill = 0  ) %>% as.data.frame()
fold.d.w = spread(fold.down, structure, log_fold_change, fill = 0  ) %>% as.data.frame()
#row.names(fold.u.w) = fold.u.w$SID 
#row.names(fold.d.w) = fold.d.w$SID
fold.d.w$SID = NULL 
fold.u.w$SID = NULL
f.u.w.drug = fold.u.w$drug
f.d.w.drug = fold.d.w$drug
fold.u.w$drug = NULL
fold.d.w$drug = NULL
library(RColorBrewer)
Colors=c("#440154FF","#3CBB75FF","#FDE725FF")
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)
col1 <- colorRampPalette(Colors)(15)
library(iheatmapr)
cl3 = c(
    "#8cffdb", "#7bb274", "#510ac9", "#ff5b00",
        "#fffe7a" , "#0a888a",  "#887191" , "#c04e01" , "#95d0fc", "#40a368" ,
        "#53fca1" , "#c04e01" , "#3f9b0b" ,"#580f41" ,  "#b9a281", "#ff474c", 
        "#7e1e9c","#0343df","#95d0fc","#f97306","#029386" ,"#c20078"
        ) 
main_heatmap(data.matrix((fold.u.w)), name = "Log10 Fold Change",) %>%
    add_row_annotation(data.frame("Drugs" = as.character(f.u.w.drug )) ,
                        colors = list("Drugs" = cl3)) %>%
    add_col_clustering()  %>%
    add_row_clustering() %>% 
    add_row_title("Drugs") %>%
    add_col_title("Folds") #%>%
    
heatmap.2(data.matrix((fold.u.w)), trace = "none",col=Colors, labCol =NA,  margins=c(3,15),  
          RowSideColors=col1[as.numeric(as.factor(f.u.w.drug))], labRow = f.u.w.drug, cexRow = 2, keysize =1 )

heatmap.2(data.matrix((fold.d.w)), trace = "none",col=Colors, labCol =NA,  margins=c(3,15),  
          RowSideColors=col1[as.numeric(as.factor(f.d.w.drug))], labRow = f.d.w.drug, cexRow = 2, keysize =1 )
library(ggplot2)
library(ggdendro)
library(plotly)
dd.col <- as.dendrogram(hclust(dist(fold.u.w)))
dd.row <- as.dendrogram(hclust(dist(t(fold.u.w))))
dx <- dendro_data(dd.row)
dy <- dendro_data(dd.col)
ggdend <- function(df) {
    ggplot() +
        geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
        labs(x = "", y = "") + theme_minimal() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank())
}
px <- ggdend(dx$segments)
py <- ggdend(dy$segments) + coord_flip()
col.ord <- order.dendrogram(dd.col)
row.ord <- order.dendrogram(dd.row)
xx = fold.u.w[col.ord, row.ord]
xx$name = row.names(xx) 
xx$name <- with(xx, factor(name, levels=name, ordered=TRUE))
mdf <- reshape2::melt(xx, id.vars="name")


p <- ggplot(mdf, aes(x = variable, y = name)) + geom_tile(aes(fill = value)) 
eaxis <- list(
    showticklabels = FALSE,
    showgrid = FALSE,
    zeroline = FALSE
)
ggplotly(p) %>%
    # note that margin applies to entire plot, so we can
    # add it here to make tick labels more readable
    layout(margin = list(l = 200))
subplot(px, p_empty, p, py, nrows = 2, margin = 0.01)


p <- plot_ly(
    x = c("a", "b", "c"), y = c("d", "e", "f"),
    z = m, type = "heatmap"
)

