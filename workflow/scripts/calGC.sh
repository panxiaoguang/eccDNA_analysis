#!/usr/bin/bash
inputFile=$1
sample=$2
outdir=$3
chromSize=$4
Genome=$5

awk -F '\t' '{printf "%s\t%s\t%s\t%s:%s-%s\n",$1,$2,$3,$1,$2,$3}' ${inputFile} | \
bedtools flank -i - -l 1.0 -r 0 -pct -g $chromSize | \
bedtools nuc -fi $Genome -bed - | \
awk -F '\t' '{print $4,$6} ' OFS='\t' > ${outdir}/${sample}.upstream.gc.txt

awk -F '\t' '{printf "%s\t%s\t%s\t%s:%s-%s\n",$1,$2,$3,$1,$2,$3}' ${inputFile} | \
bedtools nuc -fi $Genome -bed - | \
awk -F '\t' '{print $4,$6} ' OFS='\t' > ${outdir}/${sample}.ecc.gc.txt

awk -F '\t' '{printf "%s\t%s\t%s\t%s:%s-%s\n",$1,$2,$3,$1,$2,$3}' ${inputFile} | \
bedtools flank -i - -l 0 -r 1.0 -pct -g $chromSize | \
bedtools nuc -fi $Genome -bed - | \
awk -F '\t' '{print $4,$6} ' OFS='\t' > ${outdir}/${sample}.downstream.gc.txt


