genes<-read_tsv("../../Bladder/dbs/hg38.coding.bed",col_names = F)
names(genes)<-c("chr","Start","End","gene")
genes<-genes%>%
  mutate(length=End-Start)

get_genes<-function(x,dbs=genes){
  df<-read_tsv(stringr::str_glue("../CircleSeq/start_anno/{x}.startAnno.bed"),col_names = F)
  df<-df%>%
    dplyr::select(1:4,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    group_by(X8)%>%
    distinct(ecc,.keep_all = T)%>%
    summarise(count=n())%>%
    arrange(desc(count))%>%
    rename(gene=X8)%>%
    left_join(dbs,by="gene")%>%
    mutate(pct=count/length)%>%
    mutate(pct2=pct/sum(pct)*(10^6))%>%
    dplyr::select(1,8)
  names(df)[2]<-"start_ratio"
  df
}

get_genes2<-function(x,dbs=genes){
  df<-read_tsv(stringr::str_glue("../CircleSeq/start_anno/{x}.endAnno.bed"),col_names = F)
  df<-df%>%
    dplyr::select(1:4,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    group_by(X8)%>%
    distinct(ecc,.keep_all = T)%>%
    summarise(count=n())%>%
    arrange(desc(count))%>%
    rename(gene=X8)%>%
    left_join(dbs,by="gene")%>%
    mutate(pct=count/length)%>%
    mutate(pct2=pct/sum(pct)*(10^6))%>%
    dplyr::select(1,8)
  names(df)[2]<-"end_ratio"
  df
}


df3<-get_genes("cKca_11T")
df4<-get_genes2("cKca_11T")
fin2<-df3%>%inner_join(df4,by="gene")

p2<-ggscatter(fin2, x = "logs", y = "loge",
              color = "#5181B0", shape = 20, size = 1.5, 
              xlab = "Start breakpoints in genes",
              ylab = "End breakpoints in genes",
              conf.int = TRUE, 
              cor.coef = TRUE, 
              cor.coeff.args = list(method = "pearson", 
                                    label.x = 3, 
                                    label.sep = "\n")
)
ggarrange(p1,p2)

ggsave("start_vs_end.corr.pdf",width = 8,height = 3.51)