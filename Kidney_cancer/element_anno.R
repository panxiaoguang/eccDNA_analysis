# Figure3I------data ------------------------------------------------------
#proteinGenes<-read_tsv("~/Project/Bladder/dbs/hg38.coding.bed",col_names = F)%>%
#  pull(X4)
#newdb<-dbs%>%
#  tidyr::separate(X4,into = c("gene","type"),sep = ":")%>%
#  filter(gene%in%proteinGenes)%>%
#  tidyr::unite("X4",gene,type,sep = ":")
cal_element_length<-function(x){
  lgth <- read_tsv(stringr::str_glue("~/Project/Bladder/dbs/hg38.{x}.merge.bed"),col_names = F)%>%
    mutate(length=X3-X2)%>%
    summarise(total_l=sum(length))%>%
    pull(total_l)
  lgth
}
elements<-c("Gene2kbU","Gene2kbD","UTR5","UTR3","Exon","Intron","CpG")
length_corr<-tibble(type=elements,total_l=sapply(elements,cal_element_length))
total_chrom_length<-read_tsv("~/Project/甲状腺癌/old/plasma/dbs/hg38.chromo.size",col_names = F)
total_chrom_length<-sum((total_chrom_length[1:24,])$X2)
length_corr<-length_corr%>%
  mutate(pct=total_l/total_chrom_length)
length_corr<-length_corr%>%
  mutate(type=case_when(type=="UTR5" ~ "5'UTR",
                        type=="UTR3" ~ "3'UTR",
                        TRUE ~ type))
cal_element<-function(x){
  fs<-read_tsv(stringr::str_glue("element_anno/cKca_{x}.startAnno.bed"),col_names = F)
  tongji<-fs%>%
    dplyr::select(1,2,3,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    tidyr::separate_rows(X8,sep =",")%>%
    tidyr::separate(X8,into=c("gene","type"),sep=":")%>%
    dplyr::select(ecc,type)%>%
    group_by(type)%>%
    summarise(ecc_c = n_distinct(ecc))
  tongji$samples<-x
  tongji
}

infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
samples<-infos$Sample2
normals<-paste0(samples,"N")
cases<-paste0(samples,"T")

fin<-do.call('rbind',lapply(c(normals,cases),cal_element))
saveRDS(fin,"gene_elements.RDS")
### pie plot data
pie_plot_data<-fin%>%
  group_by(type)%>%
  summarise(counts=sum(ecc_c))%>%
  mutate(ratio=counts/sum(counts))

## correct by element length for barplot
## should calculate ecc numbers per sample
eccNumbers<-openxlsx::read.xlsx("Kidney_circle_epms.xlsx",sheet=1)
eccNumbers<-eccNumbers%>%
  select(1:2)%>%
  setNames(c("totalNumbers","samples"))

barplot_data<-fin%>%
  left_join(eccNumbers,by="samples")%>%
  mutate(ratio=ecc_c/totalNumbers)%>%
  left_join(length_corr,by="type")%>%
  mutate(enrichment=ratio/pct)

data2 <- barplot_data%>%
  select(samples,type,enrichment)

openxlsx::write.xlsx(list(pie_plot_data,data2),"ecc_element_enrich.xlsx")


barplotData<-data2%>%
  mutate(gp=if_else(samples%in%cases,"Tumor","NAT"))%>%
  group_by(gp,type)%>%
  summarise(md=mean_sdl(enrichment))%>%
  tidyr::unnest()

pvalue<-data2%>%
  mutate(gp=if_else(samples%in%cases,"Tumor","NAT"))%>%
  group_by(type)%>%
  rstatix::t_test(enrichment~gp,paired = T)

ggplot(barplotData,aes(x=type,y=y))+
  geom_col(aes(fill=gp),position = position_dodge2(width = 0.8,padding = 0.2),color="black")+
  geom_errorbar(aes(group=gp,ymin=ymin,ymax=ymax),
                position = position_dodge2(width = 0.1,padding = 0.2),
                linewidth=0.3)+
  scale_y_continuous(expand = c(0,0))+
  theme_pubr()+
  geom_text(aes(label=p,y=2.5),data=pvalue)

ggsave("element_diff.pdf",width =6.34 ,height = 3.16)  
