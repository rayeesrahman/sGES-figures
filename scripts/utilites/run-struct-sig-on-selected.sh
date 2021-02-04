#!/bin/bash 

for file in selected_gene_lists/*
do 
    name=$( basename $file ".txt") 
    echo $name 
    for i in `seq 10 10 30` 
    do 
        head -n $i $file >  genelist.txt 
        n=$( echo  "$name^$i" ) 
        echo $n
        ../../../Structural_Signatures_2.0/structural-signatures-2.0.sh -i genelist.txt -t both -n gn -o "$name^$i" 
        cat *-domain-enrichments.csv | grep -v "structure" | tr "^" "," >> stable_signature_data/domain-enrichments.csv 
        cat *-family-enrichments.csv | grep -v "structure" | tr "^" "," >> stable_signature_data/family-enrichments.csv
        cat *-fold-enrichments.csv | grep -v "structure" | tr "^" "," >> stable_signature_data/fold-enrichments.csv
        cat *-superfam-enrichments.csv | grep -v "structure" | tr "^" "," >> stable_signature_data/superfam-enrichments.csv
        rm genelist.txt
        rm $n* 
    done 
done  
