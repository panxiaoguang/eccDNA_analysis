minimap2 -ax map-hifi -t 10 /home/panxiaoguang/Project/DataBase/hg38.fa CCS/16N.ccs.fastq.gz > alignData/16N.align.sam
samtools sort -@ 10 -o alignData/16N.align.sort.bam alignData/16N.align.sam
samtools index -@ 10 alignData/16N.align.sort.bam
