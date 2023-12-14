library(readr)
library(dplyr)


bedFormat<-function(x){
    df<-read_tsv(stringr::str_glue("filter/{x}_circle_site.filter.tsv"),col_names = T)
    df<-df%>%
      select(1:3)%>%
      arrange(chrom,start)%>%
      mutate(X4=seq(1,nrow(.)))%>%
      mutate(X4=paste(x,X4,sep="_"))
    write_tsv(df,stringr::str_glue("filter_bed/{x}.ecc.bed"),col_names = F)
}

files <- Sys.glob("filter/*_circle_site.filter.tsv")
samples <- stringr::str_remove(files,"_circle_site.filter.tsv")
samples <- stringr::str_remove(samples,"filter/")

lapply(samples, bedFormat)
