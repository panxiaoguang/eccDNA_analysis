
## make gene tracks######################################################################
library(Homo.sapiens)
TxDb(Homo.sapiens)<-TxDb.Hsapiens.UCSC.hg38.knownGene
Homo.sapiens
get_gene_id<-function(gene.name="FGFR2"){
  AnnotationDbi::select(Homo.sapiens,keys = "FGFR2",keytype = "SYMBOL",columns = c("SYMBOL","GENEID"))$GENEID
}

fgfr2<-GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene,filter=c("gene_id"=get_gene_id("FGFR2")))
#geneTK<-autoplot(Homo.sapiens, which = fgfr2, stat = "reduce",fill="#3288bd",color="black")
########################################################################################
#####make eccDNA links

make_links<-function(sample,gene.name="FGFR2"){
  df<-readr::read_tsv(stringr::str_glue("genedrops/{sample}.startAnno.bed"),col_names=F)
  test<-df%>%
    filter(X8==gene.name)%>%
    select(1:4)%>%
    setNames(c("seqnames","start","end","name"))%>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  test$value<-runif(length(test),1,5)
  test
}

### to build eccDNA tracks
gct<-make_links("JXM-7")
nat<-make_links("JXM-20")

link1<-ggplot(gct)+geom_arch(aes(height=value),color="#31a354")
link2<-ggplot(nat)+geom_arch(aes(height=value),color="#878787")
geneTK<-ggplot(fgfr2)+geom_arrowrect(fill="#3288bd",color="black")+xlim(fgfr2+20000)
pdf("test.plot.pdf",width =8.89,height = 2.88)
tracks(gene=geneTK,GCT=link1,NAT=link2,heights = c(1,4,4))+theme_pack_panels()
dev.off()
