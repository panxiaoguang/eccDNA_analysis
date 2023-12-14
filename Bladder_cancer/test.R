### try to use ggbio to plot CNA result from gistic2
library(readr)
library(dplyr)
library(biovizBase)
library(ggbio)
library(ggplot2)
library(GenomicRanges)
### first to load data 
scores<-read_tsv("CNV/gistic2_results/scores.gistic")
### transfer the gistic2 score to GRanges obj.

scores<-scores%>%
  select(Chromosome,Start,End,Type,`G-score`)%>%
  mutate(`G-score`=if_else(Type=="Del",0-`G-score`,`G-score`))

names(scores)[1:3]<-c("seqnames","start","end")  

gr<-makeGRangesFromDataFrame(scores,keep.extra.columns = T)
seqlevelsStyle(gr)<-"UCSC"
rm(scores)
## get seqinfo 
sqinfos<-seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
sqinfos<-keepSeqlevels(sqinfos,paste0("chr",1:22))
seqinfo(gr)<-sqinfos
## get the coord::genome using BioVizBase
gr.dt<-transformToGenome(gr,space.skip=0)
bks<-getXScale(gr.dt)

ggplot(gr.dt)+
  geom_area(aes(x=.start,y=gr.dt$`G-score`,group=Type,fill=Type))+
  scale_x_continuous(breaks = bks$breaks)+
  scale_fill_manual(values = c("Amp"="#C5362C","Del"="#4475A7"))+
  theme_clear(axis.ticks.x = TRUE,axis.line.color = "black")
