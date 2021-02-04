#!/bin/bash 
rm selected-id-organ-systems.txt 
while read line 
do 
    grep $line id-organ-system.txt | head -10  >> selected-id-organ-systems.txt 
done < organsystems.txt 

while read line 
do 
    id=$(echo $line | cut -f1 -d",") 
    organ=$(echo $line | cut -f2 -d",") 
    echo $line 
    grep "^$id" ../allcombined.txt.2 | cut -f 4 -d"," > selected_gene_lists/$id^$organ.txt
done < selected-id-organ-systems.txt 
