library(GenomicRanges)
library(readr)
library(dplyr)

makeGranges<-function(x,sample=NULL){
  df<-read_tsv(x,col_names=F)%>%
    select(1:3)%>%
    setNames(c("seqname","start","end"))%>%
    mutate(eccName=paste0(sample,"_",1:nrow(.)))
  gr <- makeGRangesFromDataFrame(df,keep.extra.columns = T,ignore.strand = T)
  gr <- keepStandardChromosomes(gr,species = "Homo_sapiens",pruning.mode = "coarse")
  tmpinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  tmpinfo <- keepStandardChromosomes(tmpinfo,species = "Homo_sapiens")
  seqinfo(gr)<-tmpinfo
  gr
}

### test
##################make protein coding gene from txdb#####################################
gr <- makeGranges("bedFile/JXM-1.ecc.bed",sample = "JXM-1")
database<-genes(TxDb.Hsapiens.UCSC.hg38.knownGene,single.strand.genes.only=F)
genetypes<-select(Homo.sapiens,keys = keys(Homo.sapiens,keytype = "SYMBOL"),keytype = "SYMBOL",columns = c("GENEID","GENETYPE","SYMBOL"))
genetypes<-genetypes%>%
  filter(GENETYPE=="protein-coding")%>%
  filter(is.na(SYMBOL))
protein_Codings<-database[names(database)%in%(genetypes$GENEID)]
protein_Codings$GENEID=names(protein_Codings)
shuju<-as_tibble(values(protein_Codings))%>%
  left_join(genetypes,by="GENEID")%>%
  pull(SYMBOL)
protein_Codings$SYMBOL<-shuju
save(protein_Codings,file = "../DataBase/proteinCoding.Rdata")
##########################################################################################
anno_start<-function(x,y){
  x_tmp<-narrow(x,start=1,width = 1)
  res<-findOverlaps(x_tmp,y,ignore.strand=T)
  jiaoji<-x[queryHits(res)]
  mcols(jiaoji)$anno<-mcols(y[subjectHits(res)])[[2]]
  jiaoji
}


eccDNA<-makeGranges("qd-ECC4/S/ECC_report/FinallyData/bedFile/S10.ecc.bed",meta.name = "eccDNA")

jiaoji<-anno_start(gr,protein_Codings)

get_overlap_percentage<-function(x,y,query.restrict=T,pct=0.2){
  hits <- findOverlaps(x,y,ignore.strand=T)
  ints<-pintersect(x[queryHits(hits)],y[subjectHits(hits)])
  if(query.restrict==TRUE){
    int_p<-width(ints)/width(x[queryHits(hits)])
  }else{
    int_p<-width(ints)/width(y[subjectHits(hits)])
  }
  hits <- hits[int_p>=pct]
  jiaoji<-x[queryHits(hits)]
  mcols(jiaoji)$anno<-mcols(y[subjectHits(hits)])[[1]]
  jiaoji
}

get_overlap_percentage(eccDNA,genes)


