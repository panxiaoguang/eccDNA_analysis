library(readr)
library(dplyr)
TPM<-read_tsv("count_TPM.txt")
names(TPM)[1]<-"Geneid"
goodGenes<-TPM%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  mutate(type=if_else(TPM>3,"good","bad"))%>%
  group_by(Geneid,type)%>%
  summarise(count=n())%>%
  tidyr::pivot_wider(names_from = "type",values_from = count,values_fill = 0)%>%
  filter(good==230)%>%
  pull(Geneid)
TPM3<-TPM%>%
  filter(Geneid%in%goodGenes)%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  filter(stringr::str_ends(sample,"T"))

stst<-TPM3%>%
  group_by(Geneid)%>%
  summarise(mean=mean(TPM),sd=sd(TPM))
## zscore
zscores<-TPM3%>%
  left_join(stst,by="Geneid")%>%
  mutate(zscore=(TPM-mean)/sd)

zscores<-zscores%>%
  dplyr::select(Geneid,sample,zscore)%>%
  dplyr::rename(TPM=zscore)%>%
  mutate(sample=stringr::str_extract(sample,"\\d+"))%>%
  mutate(sample=as.numeric(sample))%>%
  mutate(sample=as.character(sample))

### define which gene been effected by eccDNA(overlap 20%)

## readin genes
allGenes<-read_tsv("../../Bladder/dbs/hg38.coding.bed",col_names = F)
names(allGenes)<-c("seqnames","start","end","genes")
genes<-makeGRangesFromDataFrame(allGenes,keep.extra.columns = T)
rm(allGenes)

find_ecc_genes<-function(x){
  df<-read_tsv(stringr::str_glue("../CircleSeq/filter_bed/cKca_{x}T.ecc.bed"),col_names = F)
  names(df)<-c("seqnames","start","end","eccDNA")
  gr<-makeGRangesFromDataFrame(df,keep.extra.columns = T)
  rm(df)
  hits <- findOverlaps(gr,genes,ignore.strand=T)
  ints<-pintersect(gr[queryHits(hits)],genes[subjectHits(hits)])
  int_p<-width(ints)/width(genes[subjectHits(hits)])
  hits <- hits[int_p>=0.2]
  hitGenes<-(genes[subjectHits(hits)])$genes
  tibble(sample=x,Geneid=goodGenes)%>%
    mutate(eccDNA=if_else(Geneid%in%hitGenes,"T","F"))
}

fin<-do.call("rbind",lapply(unique(zscores$sample),find_ecc_genes))

zscores<-zscores%>%left_join(fin,by=c("sample","Geneid"))

ggplot(zscores,aes(x=TPM,color=eccDNA))+
  geom_density()+
  scale_color_manual(values = c("T"="#9E3735","F"="#48436D"))+
  theme_pubr()+
  xlab("RNA expression(Z-score)")

ggsave("Non_vs_eccDNA.pdf",width = 5.52,height = 2.97)

ggplot(zscores,aes(x=eccDNA,y=TPM,color=eccDNA))+
  geom_violin()+
  scale_color_manual(values = c("T"="#9E3735","F"="#48436D"))+
  theme_pubr()+
  xlab("")
ggsave("Non_vs_eccDNA_violin.pdf",width = 3.98,height = 3.47)

saveRDS(zscores,"Non_vs_eccDNA.RDS")

############################new type of plot

p1<-ggplot(ecc_exp,aes(x=TPM))+geom_density(n=80)+
  scale_x_continuous(limits = c(-3,6))
data2<-layer_data(p1,1)

ggplot(ecc_exp,aes(x=TPM))+
  geom_bar(aes(x=x,y=density),
           color="black",
           fill="white",
           linewidth=0.2,
           stat = "identity",
           position = position_dodge(0),
           data = data2)+
  geom_density(aes(color=eccDNA),n=80)+
  scale_x_continuous(limits = c(-3,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values = c("T"="#DA60DE","F"="#ECBF5D"))+
  theme_pubr()

ggsave("Non_vs_eccDNA.pdf",width = 7.52,height =3.31 )

ecc_exp<-readRDS("Non_vs_eccDNA.RDS")
