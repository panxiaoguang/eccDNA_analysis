
allGenes<-read_tsv("../../Bladder/dbs/hg38.coding.bed",col_names = F)
names(allGenes)<-c("seqnames","start","end","genes")
genes<-makeGRangesFromDataFrame(allGenes,keep.extra.columns = T)
rm(allGenes)

df<-read_tsv("eccDNAs/cKca-14T.info.ft.txt")
df$label<-seq(1,nrow(df))

df%>%
  select(label,fragments)%>%
  tidyr::separate_rows(fragments,sep="\\|")%>%
  tidyr::separate(fragments,into=c("seqnames","start","end","strand"))
