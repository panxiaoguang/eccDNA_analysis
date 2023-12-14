infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2

cases<-paste0(need_samples,"T")
normals<-paste0(need_samples,"N")

eccNumbers<-openxlsx::read.xlsx("Kidney_circle_epms.xlsx",sheet=1)
eccNumbers<-eccNumbers%>%
  select(1:2)%>%
  setNames(c("totalNumbers","sample"))

cal_ingenes<-function(x){
  read_tsv(stringr::str_glue("start_anno/cKca_{x}.startAnno.bed"),col_names = F)%>%
    select(X4)%>%
    distinct(X4)%>%
    nrow()
}
genenumbers<-tibble(sample=c(cases,normals),eccs=sapply(c(cases,normals),cal_ingenes))
haha<-genenumbers%>%
  left_join(eccNumbers,by="sample")

haha<-haha%>%
  mutate(geneP=eccs/totalNumbers*100)

openxlsx::write.xlsx(haha,"ecc_in_genes.xlsx")
