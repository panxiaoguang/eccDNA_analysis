
df<-read_tsv("../../Bladder/dbs/hg38.genetic_elements.bed",col_names = F)
## split into different type

tp1 <- df%>%
  filter(stringr::str_detect(X4,"Gene2kbU"))
tp2 <- df%>%
  filter(stringr::str_detect(X4,"Gene2kbD"))
tp3 <- df%>%
  filter(stringr::str_detect(X4,"3'UTR"))
tp4 <- df%>%
  filter(stringr::str_detect(X4,"5'UTR"))
tp5 <- df%>%
  filter(stringr::str_detect(X4,"Exon"))
tp6 <- df%>%
  filter(stringr::str_detect(X4,"Intron"))
tp7 <- df%>%
  filter(stringr::str_detect(X4,"CpG"))

write_tsv(tp1,"../../Bladder/dbs/hg38.Gene2kbU.bed",col_names = F)
write_tsv(tp2,"../../Bladder/dbs/hg38.Gene2kbD.bed",col_names = F)
write_tsv(tp3,"../../Bladder/dbs/hg38.3'UTR.bed",col_names = F)
write_tsv(tp4,"../../Bladder/dbs/hg38.5'UTR.bed",col_names = F)
write_tsv(tp5,"../../Bladder/dbs/hg38.Exon.bed",col_names = F)
write_tsv(tp6,"../../Bladder/dbs/hg38.Intron.bed",col_names = F)
write_tsv(tp7,"../../Bladder/dbs/hg38.CpG.bed",col_names = F)

### use bedtools merge to merge beds

## then integrate
tp1<-read_tsv("../../Bladder/dbs/hg38.Gene2kbU.merge.bed",col_names = F)
tp2<-read_tsv("../../Bladder/dbs/hg38.Gene2kbD.merge.bed",col_names = F)
tp3<-read_tsv("../../Bladder/dbs/hg38.UTR3.merge.bed",col_names = F)
tp4<-read_tsv("../../Bladder/dbs/hg38.UTR5.merge.bed",col_names = F)
tp5<-read_tsv("../../Bladder/dbs/hg38.Exon.merge.bed",col_names = F)
tp6<-read_tsv("../../Bladder/dbs/hg38.Intron.merge.bed",col_names = F)
tp7<-read_tsv("../../Bladder/dbs/hg38.CpG.merge.bed",col_names = F)

fin<-rbind(tp1,tp2,tp3,tp4,tp5,tp6,tp7)

fin<-fin%>%
  arrange(X1,X2)

write_tsv(fin,"../../Bladder/dbs/hg38.genetic_elements.merge.bed",col_names = F)
