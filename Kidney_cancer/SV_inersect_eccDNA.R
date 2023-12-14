## load eccDNA to gRanges
library(readr)
library(dplyr)
library(GenomicRanges)

get_inters<-function(x){
  print(x)
  ##make eccDNA into gR
  eccDNA<-read_tsv(stringr::str_glue("filter_bed/cKca_{x}T.ecc.bed"),col_names = F)
  names(eccDNA)<-c("seqnames","start","end","eccDNA")
  ecc_gr<-makeGRangesFromDataFrame(eccDNA,keep.extra.columns = T)
  ##make SV effect inro gR from the same sample
  SVs<-read_tsv(stringr::str_glue("../WGS/SV/merged/{x}.merged.allsv.filter.txt"))
  SVs<-SVs%>%
    filter(svclass!="TRA")%>%
    select(chrom1,start1,start2,svclass)%>%
    setNames(c("seqnames","start","end","svclass"))
  SV_gr<-makeGRangesFromDataFrame(SVs,keep.extra.columns = T)
  itrs<-findOverlaps(SV_gr,ecc_gr)
  jiaoji<-SV_gr[queryHits(itrs)]
  jiaoji$eccDNA<-ecc_gr[subjectHits(itrs)]
  jiaoji_df<-as.data.frame(jiaoji)
  jiaoji_df<-as_tibble(jiaoji_df)
  stat<-jiaoji_df%>%
    mutate(haha=paste(seqnames,start,end))%>%
    distinct(haha,.keep_all = T)%>%
    count(svclass)
  stat$sample<-x
  stat
}

infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2


fin<-do.call("rbind",lapply(as.character(need_samples),get_inters))
names(fin)<-c("type","numbers","samples")
fin<-fin%>%group_by(type)%>%summarise(totaln=sum(numbers))

write_tsv(fin,"merged_SV.eccDNA_inters.statistic.tsv")
#write_tsv(fin,"Non-merged_SV.eccDNA_inters.statistic.tsv")


## find diff
# compare amp with sv no-sv

get_diff<-function(x){
  need_chrom<-c(paste0("chr",seq(1,22)),"chrX","chrY")
  total_chrom_length <- 3088269832
  eccDNA<-read_tsv(stringr::str_glue("filter_bed/cKca_{x}T.ecc.bed"),col_names = F)
  names(eccDNA)<-c("seqnames","start","end","eccDNA")
  eccDNA<-eccDNA%>%
    filter(seqnames%in%need_chrom)
  ecc_gr<-makeGRangesFromDataFrame(eccDNA,keep.extra.columns = T)
  ##make SV effect inro gR from the same sample
  SVs<-read_tsv(stringr::str_glue("../WGS/SV/merged/{x}.merged.allsv.txt"))
  SVs<-SVs%>%
    filter(svclass!="TRA")%>%
    filter(chrom1%in%need_chrom)%>%
    select(chrom1,start1,start2,svclass)%>%
    setNames(c("seqnames","start","end","svclass"))
  SV_gr<-makeGRangesFromDataFrame(SVs,keep.extra.columns = T)
  itrs<-findOverlaps(SV_gr,ecc_gr)
  numbers_eccDNA_withSV<-length(unique(subjectHits(itrs)))
  numbers_eccDNA_noSV<-length(ecc_gr)-numbers_eccDNA_withSV
  ### normalize
  total_SV_length<-sum(width(SV_gr))
  fin<-tibble(type=c("inSV","outSV"),
              numbers=c(numbers_eccDNA_withSV,numbers_eccDNA_noSV),
              length=c(total_SV_length,total_chrom_length-total_SV_length))%>%
    mutate(ratio=numbers/length,ratio2=ratio/sum(ratio))
  fin$sample<-x
  fin
}

infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2

wocao<-do.call("bind_rows",lapply(need_samples, get_diff))


plotData<-wocao%>%
  select(type,ratio2,sample)%>%
  tidyr::pivot_wider(names_from="sample",values_from = "ratio2")

write_tsv(plotData,"Non-merged_SV.eccDNA_numbers.tsv")
