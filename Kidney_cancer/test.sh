source /home/panxiaoguang/miniconda3/bin/activate biosoft
### calculate coverage "deeptools"
bamCoverage -b Kca_45T_sorted_reheader.bam --region chr6:36433634:44493758 -of bedgraph -o RNA.cov.bdg
bamCoverage -b 45T.cs.rmdup.sort.bam --region chr6:36433634:44493758 -of bedgraph -o WGS.cov.bdg
bamCoverage -b sorted_cKca_45T_circle.bam --region chr6:36433634:44493758 -of bedgraph -o CIRCLE.cov.bdg

source /home/panxiaoguang/miniconda3/bin/activate allelecount
### calculate allele count "allelecounter"
alleleCounter -l snpPOS.txt -b Kca_45T_sorted_reheader.bam -o RNA.allele.txt -r /home/panxiaoguang/Project/DataBase/hg38.fa.fai -f 1 -F 1548 -m 0 -q 0
alleleCounter -l snpPOS.txt -b 45T.cs.rmdup.sort.bam -o WGS.allele.txt -r /home/panxiaoguang/Project/DataBase/hg38.fa.fai
alleleCounter -l snpPOS.txt -b sorted_cKca_45T_circle.bam -o CIRCLE.allele2.txt -r /home/panxiaoguang/Project/DataBase/hg38.fa.fai -f 1 -F 1548 -m 20 -q 10