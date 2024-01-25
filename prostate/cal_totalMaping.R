

## calculate total mappings

getTotalMapReads<-function(sample){
  df<-read_tsv(stringr::str_glue("Total.Mapping.Reads/{sample}_mapping.flag.txt"),col_names = F)
  df%>%
    filter(X1!="*" )%>%
    filter(X1!="chrM")%>%
    summarise(totalMaps=sum(X3))%>%
    pull(totalMaps)
}

newTb<-tibble(sample=c(plasma_samples,normal_samples),
              totalmappings=as.numeric(sapply(c(plasma_samples,normal_samples),getTotalMapReads)))

write_tsv(newTb,"FinallyData/urine/eccmaps.tsv")


### find miss files

NormalInfo<-read.xlsx("sampleInfo.xlsx",sheet = 2)
NormalInfo<-NormalInfo%>%
  as_tibble()
normal_sanples<-NormalInfo%>%
  filter(!is.na(SampleID))%>%
  pull(SampleName)

nowhave<-Sys.glob("Total.Mapping.Reads/*_mapping.flag.txt")
nowhave<-nowhave%>%
  stringr::str_remove("Total.Mapping.Reads/")%>%
  stringr::str_remove("_mapping.flag.txt")

nd<-setdiff(normal_sanples,nowhave)
write_tsv(tibble(samples=nd),"lost_mapping_data.tsv")
