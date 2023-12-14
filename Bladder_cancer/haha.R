new<-read_tsv("RNAseq/Protein_coding.tpms.txt")
old<-openxlsx::read.xlsx("下机信息（含barcode）/#Table S4_guang.xlsx",sheet = 3)
colnames(old)[1]<-"Symbol"
need_genes1<-openxlsx::read.xlsx("Cancer_driver genes_Bladder cancer.xlsx",sheet = 1)
oncogenes<-need_genes1$Gene

sel_samples<-names(new)%>%sample(size = 16)

jiancha<-function(x){
  tmp1<-new%>%select(Symbol,all_of(x))
  tmp2<-old%>%select(Symbol,all_of((paste0("CCGA-UBC-",stringr::str_remove(x,"^0+")))))
  tmp1%>%
    left_join(tmp2,by="Symbol")%>%
    mutate(col=if_else(Symbol%in%oncogenes,"Oncogenes","others"))%>%
    mutate(col=forcats::fct_relevel(col,c("others","Oncogenes")))%>%
    setNames(c("Symbol","xvar","yvar","col"))%>%
    plot_ly(x=~xvar,y=~yvar)%>%
    add_trace(text=~Symbol,color=~col,colors=c("#d9d9d9","#e41a1c"),hovertemplate ="Gene:%{text}<br>new:%{x}</br>old:%{y}")%>%
    layout(xaxis=list(title="New result"),yaxis=list(title="Old calculating"),showlegend=F)
}

jiancha("008T")

fig_lst<-lapply(sel_samples,jiancha)
subplot(fig_lst,nrows = 2)
