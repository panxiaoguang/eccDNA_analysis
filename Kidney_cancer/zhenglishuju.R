preprocess<-function(x){
  df<-read_tsv(stringr::str_glue("CNV/anno.cnv/{x}T.anno.cns"))
  df%>%
    select(gene,log2)%>%
    tidyr::separate_rows(gene,sep = ",")%>%
    filter(gene!="-")%>%
    arrange(gene,desc(log2))%>%
    setNames(c("Geneid",paste0(x,"T")))%>%
    group_by(Geneid)%>%
    slice_head(n=1)%>%
    distinct(Geneid,.keep_all = T)
}



samples<-c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 11L, 
           12L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 
           21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 
           29L, 30L, 31L, 32L, 33L, 34L, 35L, 36L, 
           38L, 39L, 45L, 46L, 47L, 48L, 49L, 50L, 
           51L, 52L, 53L, 54L, 55L, 56L, 57L, 58L, 
           59L, 60L, 61L, 62L, 63L, 64L, 65L, 66L, 
           68L, 70L, 71L, 72L, 73L, 74L, 75L, 76L, 
           78L, 79L, 80L, 81L, 82L, 84L, 85L, 86L, 
           87L, 88L, 89L, 90L, 91L, 92L, 93L, 94L, 
           95L, 96L, 97L, 98L, 99L, 100L, 102L, 103L, 
           104L, 105L, 106L, 107L, 109L, 110L, 111L, 
           112L, 113L, 114L, 115L, 116L, 117L, 118L, 
           119L, 120L, 121L, 122L, 123L, 124L, 125L, 
           126L, 127L, 128L, 129L, 130L, 131L, 132L, 
           133L, 134L, 136L, 138L)

samples<-as.character(samples)

CNVs<-Reduce(function(x,y){full_join(x,y,by="Geneid")},lapply(samples,preprocess))

names(CNVs)<-c("Geneid",paste0("CCGA-RCC-",stringr::str_pad(names(CNVs)[-1],width = 4,side = "left",pad = "0")))

openxlsx::write.xlsx(CNVs,"CNVs.xlsx")

### SV


process_SV<-function(x){
  tmp<-read_tsv(stringr::str_glue("SV/merged/{x}.merged.allsv.filter.txt"))%>%
    select(chrom1,start1,chrom2,start2,strand1,strand2,svclass,svmethod)%>%
    mutate(strand=paste0(strand1,"/",strand2))%>%
    select(1:4,9,7,8)
  tmp$sample<-paste0("CCGA-RCC-",stringr::str_pad(x,width = 4,side = "left",pad = "0"))
  tmp
}

samples<-read_tsv("SV/samples",col_names = F)

fin<-do.call("rbind",lapply(samples$X1, process_SV))

openxlsx::write.xlsx(fin,"SVs.xlsx")
