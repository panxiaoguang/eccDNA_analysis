library(GenomicRanges)
library(readr)
library(dplyr)

df<-read_tsv("~/Project/gastric/WGS/merged/J1.merged.allsv.filter.txt")
df<-df%>%
  select(1,2,5,9,10,11)%>%
  filter(svclass!="TRA")%>%
  mutate(strand=paste0(strand1,"/",strand2))%>%
  select(chrom1,start1,start2,strand)%>%
  setNames(c("seqnames","start","end","type"))
  
grSV<-makeGRangesFromDataFrame(df,keep.extra.columns = T)

hits<-findOverlaps(grSV,GRanges(seqnames = "chr10",IRanges(start = 22538311,end = 124834538)))
regionGR<-grSV[queryHits(hits)]
#regionGR$value<-width(regionGR)/(129940547-1042086)

######################################################
cna<-read_tsv("~/Project/gastric/WGS/CNVkit/J1.cs.rmdup.sort.call.cns")

cn_chr <- purrr::map(split(cna$cn,cna$chromosome),~Rle(.x))
idx_chr <- do.call("rbind",purrr::map(names(cn_chr),~tibble(chroms=.x,starts=start(cn_chr[[.x]]),ends=end(cn_chr[[.x]]))))

mg <- function(x,y,z){
  tmp <- cna%>%
    filter(chromosome==z)
  df <- tmp[x:y,]
  chr <- df$chromosome[1]
  starts <- df$start
  ends <- df$end
  start <- starts[1]
  end <- ends[length(starts)]
  cn <- df$cn[1]
  tibble(chromsome=chr,start=start,end=end,cn=cn)
}
cn_fin<-purrr::pmap_dfr(list(idx_chr$starts,idx_chr$ends,idx_chr$chroms),mg)
write_tsv(cn_fin,"~/Project/gastric/WGS/CNVkit/J12.cs.rmdup.sort.call.merge.cns")


names(cn_fin)[1]<-"seqnames"

grCN<-makeGRangesFromDataFrame(cn_fin,keep.extra.columns = T)

hits<-findOverlaps(grCN,GRanges(seqnames = "chr10",IRanges(start = 22538311,end = 124834538)))
regionGR2<-grCN[queryHits(hits)]




p2<-ggplot(regionGR)+
  geom_arch(aes(color=type))+
  theme_pack_panels()
p3<-ggplot(regionGR2)+
  geom_rect(aes(y=cn),stat="identity")+
  theme_pack_panels()


tracks(sv=p2,cn=p3)+
  scale_x_continuous(labels = function(x){paste0(x/1000/1000,"Mb")},n.breaks = 6)



