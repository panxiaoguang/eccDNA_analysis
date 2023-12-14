library(readr)
library(dplyr)

hebing<-function(sample){
    df1<-read_tsv(stringr::str_glue("format_svaba/format.{sample}.sv.svaba.txt"))
    df2<-read_tsv(stringr::str_glue("format_svaba/format.{sample}.smallIndel.svaba.txt"))
    df2$sv_id<-as.character(df2$sv_id)
    fin<-bind_rows(df1,df2)%>%
        arrange(chrom1,start1)
    write_tsv(fin,stringr::str_glue("format_svaba/format.{sample}.allsvs.svaba.txt"))
}

samples<-read_tsv("samples",col_names=F,col_types=cols("c"))%>%pull(X1)
tk <- lapply(samples,hebing)