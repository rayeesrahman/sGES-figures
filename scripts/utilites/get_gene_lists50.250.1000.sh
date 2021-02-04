#!/bin/bash 
###get top 50, 250 and 1000 genes per GTEX experiment 

while read line 
do 
    echo $line 
    grep "^$line" ../allcombined.txt.2 > $line.tmp 
    cat $line.tmp | head -50 >> top.50.genes.csv 
    cat $line.tmp | head -250 >> top.250.genes.csv
    rm $line.tmp 
done < ../unique.samples.txt

