infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2

cases<-paste0(need_samples,"T")
normals<-paste0(need_samples,"N")

cases<-c("141T","142T")
normals<-c("141N","142N")
readGC<-function(x){
  df<-read_tsv(stringr::str_glue("filter_bed/GCs/cKca_{x}.gcContents.txt"))%>%
    tidyr::gather(type,value,-ecc)%>%
    mutate(type2=case_when(type=="self" ~ "eccDNA",
                           type=="downstream" ~ "Downstream",
                           type=="upstream" ~ "Upstream"))
  df$sample<-x
  df
}

plotData<-do.call("bind_rows",lapply(c(cases,normals),readGC))

saveRDS(plotData,"ecc.gc.Rds")

eccGC<-plotData%>%
  filter(type2=="eccDNA")


avgGC<-eccGC%>%
  group_by(sample)%>%
  summarise(avgGC=mean(value))

write_tsv(avgGC,"ecc.averageGC.tsv")
