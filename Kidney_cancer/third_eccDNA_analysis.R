## calculate intersect numbers between third and next

library(readr)
library(dplyr)
samples<-c("cKca-14N", "cKca-14T", "cKca-15N", "cKca-15T", "cKca-18N", "cKca-18T", "cKca-19N", "cKca-19T", "cKca-1N", "cKca-1T", "cKca-23N", "cKca-23T", "cKca-25N", "cKca-25T", "cKca-2N", "cKca-2T", "cKca-45N", "cKca-45T", "cKca-8N", "cKca-8T")


# filter ------------------------------------------------------------------

ft<-function(x){
  test<-read_tsv(stringr::str_glue("eccDNAs/{x}.info.txt"))
  test<-test%>%
    filter(Nfullpass>=2)%>%
    distinct(fragments,.keep_all = T)
  write_tsv(test,stringr::str_glue("eccDNAs/{x}.info.ft.txt"))
}

lapply(samples, ft)


# expand eccDNA 5bp -------------------------------------------------------
bedFormat<-function(x){
  df<-read_tsv(stringr::str_glue("eccDNAs/{x}.info.txt"))
  haha<-df%>%
    filter(Nfullpass>=2)%>%
    mutate(eccDNAs=fragments)%>%
    filter(!stringr::str_detect(fragments,"\\|"))%>%
    select(fragments,eccDNAs)%>%
    mutate(fragments=stringr::str_remove(fragments,"\\([+-]\\)"))%>%
    distinct(fragments,.keep_all=T)%>%
    tidyr::separate(fragments,into=c("chrom","haha"),sep=":")%>%
    tidyr::separate(haha,into=c("start","end"),sep="-")%>%
    arrange(chrom,start)
  
  haha$start<-as.numeric(haha$start)
  haha$end<-as.numeric(haha$end)
  haha<-haha%>%mutate(start=start-6,end=end+5)
  haha<-haha%>%
    mutate(start=if_else(start<0,0,start))
  write_tsv(haha,stringr::str_glue("eccDNAs/{x}.eccDNA_expand.bed"),col_names=F)
}

lapply(samples,bedFormat)

## bedtools intersect -a ../CircleSeq/filter_bed/cKca_1N.ecc.bed -b  eccDNAs/cKca-1N.eccDNA_expand.bed

# remove duplicated results from intersetcs -------------------------------
de_dup <- function(x){
  df<-read_tsv(stringr::str_glue("eccDNAs/{x}.intersect.bed"),col_names = F)
  df<-df%>%
    arrange(X4,desc(X9))%>%
    group_by(X4)%>%
    slice_head(n=1)%>%
    ungroup()
  write_tsv(df,stringr::str_glue("eccDNAs/{x}.intersect.ft.bed"),col_names = F)
}
lapply(samples,de_dup)
# intersect ---------------------------------------------------------------

cal_nums<-function(x){
  df<-read_tsv(stringr::str_glue("eccDNAs/{x}.intersect.ft.bed"),col_names = F)
  eccDNA_short<-nrow(df)
  eccDNA_long<-read_tsv(stringr::str_glue("eccDNAs/{x}.eccDNA_expand.bed"),col_names = F)%>%nrow()
  eccDNA_inters<-df%>%
    filter(X5!=".")%>%
    distinct(X8,.keep_all = T)%>%
    nrow()
  tibble(sample=x,eccDNALong=eccDNA_long,eccDNAShort=eccDNA_short,eccDNAInters=eccDNA_inters)
}

fin<-do.call("bind_rows",lapply(samples, cal_nums))
write_tsv(fin,"eccDNA_intersect.statistic.txt")


# length ------------------------------------------------------------------

getLength_long<-function(x){
  df<-read_tsv(stringr::str_glue("eccDNAs/cKca-{x}.info.txt"))
  tmp<-df%>%
    filter(Nfullpass>=2)%>%
    distinct(fragments,.keep_all = T)%>%
    select(seqLength)
  tmp$sample<-x
  tmp
}

chang_eccDNA<-do.call("bind_rows",lapply(stringr::str_remove(samples,"cKca-"), getLength_long))

getLength_duan<-function(x){
  df<-read_tsv(stringr::str_glue("../CircleSeq/filter/cKca_{x}_circle_site.filter.tsv"))
  tibble(seqLength=df$length,sample=x)
}

duan_eccDNA<-do.call("bind_rows",lapply(stringr::str_remove(samples,"cKca-"),getLength_duan))
duan_eccDNA$type<-"short_read"
chang_eccDNA$type<-"long_read"
fin<-bind_rows(duan_eccDNA,chang_eccDNA)
ggplot(fin,aes(x=seqLength))+
  geom_density(aes(fill=type),color="black",alpha=0.6)+
  theme_pubr()+
  xlab("The length distribution of eccDNA")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,2000))+
  #scale_y_continuous(limits = c(0,20000))+
  scale_fill_manual(values = c("#4475A7","#C5362C"))
ggsave("length.pdf",width =5.17 ,height =2.97)

wocao<-fin%>%
  group_by(sample,type)%>%
  summarise(meanSize=mean(seqLength))%>%
  ungroup()%>%
  tidyr::pivot_wider(names_from = "type",
                     values_from = "meanSize")

write_tsv(wocao,"ecc_average_length.txt")


# plot juction ----------------------------------------------------------------
need_chroms<-c(paste0("chr",seq(1,22)),"chrX","chrY")
joint_links<-read_tsv("eccDNAs/cKca-14N.breakpoints.txt",col_names = F)
names(joint_links)<-c("chrm1","pos1","chrm2","pos2")
joint_links<-joint_links%>%
  filter(chrm1%in%need_chroms)%>%
  filter(chrm2%in%need_chroms)

expands_joint_links<-joint_links%>%
  mutate(start1=pos1-1,end1=pos1,start2=pos2-1,end2=pos2)%>%
  select(chrm1,start1,end1,chrm2,start2,end2)
bed1<-expands_joint_links%>%
  select(1:3)
bed2<-expands_joint_links%>%
  select(4:6)
names(bed1)<-c("chrm","start","end")
names(bed2)<-c("chrm","start","end")
pdf("test.pdf",width = 12,height =6.57)
circos.initializeWithIdeogram(plotType = c("axis", "labels"),species = "hg38")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
}, track.height = 0.05, bg.border = NA)
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
                   border = NA)
dev.off()


# events ------------------------------------------------------------------

get_events<-function(x){
  df<-read_tsv(stringr::str_glue("eccDNAs/{x}.info.txt"))%>%
    mutate(frags=stringr::str_remove(fragments,"\\(+\\)"))%>%
    mutate(frags=stringr::str_remove(fragments,"\\(-\\)"))%>%
    filter(Nfullpass>=2)%>%
    group_by(frags)%>%
    summarise(count=n())%>%
    mutate(count2=if_else(count>10,">10",as.character(count)))%>%
    mutate(count2=if_else(count>10,">10",as.character(count)))%>%
    group_by(count2)%>%
    summarise(number=n())
  names(df)<-c("Events","numbers")
  df$sample<-x
  df
}

fin<-do.call("rbind",lapply(samples,get_events))

fin2<-fin%>%tidyr::pivot_wider(names_from = "sample",values_from = "numbers",values_fill = 0)
fin3<-fin%>%group_by(Events)%>%summarise(totals=sum(numbers))

cal_totals<-function(x){
  df<-read_tsv(stringr::str_glue("eccDNAs/{x}.info.txt"))%>%
    mutate(frags=stringr::str_remove(fragments,"\\(+\\)"))%>%
    mutate(frags=stringr::str_remove(fragments,"\\(-\\)"))%>%
    filter(Nfullpass>=2)%>%
    group_by(frags)%>%
    summarise(count=n())%>%
    nrow()
  df
  
}
totals<-tibble(sample=samples,numbers = sapply(samples, cal_totals))

# qianheti ----------------------------------------------------------------
get_nsegment<-function(x){
  test<-read_tsv(stringr::str_glue("eccDNAs/{x}.info.txt"))
  test<-test%>%
    filter(Nfullpass>=2)%>%
    distinct(fragments,.keep_all = T)
  test$label<-seq(1,nrow(test))
  haha<-test%>%
    select(label,fragments)%>%
    tidyr::separate_rows(fragments,sep="\\|")%>%
    mutate(orign=sapply(stringr::str_split(fragments,":"),function(x) x[[1]]))%>%
    group_by(label)%>%
    summarise(count=n())%>%
    group_by(count)%>%
    summarise(num=n())
  haha$sample<-x
  haha
}

plotData<-do.call("bind_rows",lapply(samples, get_nsegment))
plotData<-plotData%>%tidyr::pivot_wider(names_from = count,values_from = num,values_fill = 0)
write_tsv(plotData,"nsegemnt.chinmeric.txt")
