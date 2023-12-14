sort -k1,1 -k2,2n hg38.Gene2kbU.bed|bedtools merge -i - -c 4 -o distinct > hg38.Gene2kbU.merge.bed
sort -k1,1 -k2,2n hg38.Gene2kbD.bed|bedtools merge -i - -c 4 -o distinct > hg38.Gene2kbD.merge.bed
sort -k1,1 -k2,2n hg38.UTR3.bed|bedtools merge -i - -c 4 -o distinct > hg38.UTR3.merge.bed
sort -k1,1 -k2,2n hg38.UTR5.bed|bedtools merge -i - -c 4 -o distinct > hg38.UTR5.merge.bed
sort -k1,1 -k2,2n hg38.Exon.bed|bedtools merge -i - -c 4 -o distinct > hg38.Exon.merge.bed
sort -k1,1 -k2,2n hg38.Intron.bed|bedtools merge -i - -c 4 -o distinct > hg38.Intron.merge.bed
sort -k1,1 -k2,2n hg38.CpG.bed|bedtools merge -i - -c 4 -o distinct > hg38.CpG.merge.bed
