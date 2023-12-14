### calculate TPM(https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
counts<-read_tsv("RNAseq/allSample.count.txt")
## gene id transform
txname<-read_tsv("RNAseq/Allgenes_meta.txt")
## remove id version
counts$Geneid<-stringr::str_remove(counts$Geneid,"\\.\\d+")
## transfer
txname<-txname%>%
  dplyr::select(1,5)%>%
  setNames(c("Geneid","Symbol"))%>%
  na.omit()
counts<-counts%>%left_join(txname,by="Geneid")%>%
  dplyr::select(1,129,2:128)%>%
  na.omit()

##keep only protein coding genes [19137]
load("StableProteinCodings.Rdata")
counts<-counts%>%
  filter(Symbol%in%(StableProteinGenes))
## random drop duplicated genes["TBCE","HERC3","POLR2J3","PINX1","LINC02203","GOLGA8M","SIGLEC5"],rest=19135
counts<-counts%>%
  distinct(Symbol,.keep_all = T)
### get length 
geneLength<-read_tsv("RNAseq/Gene_length.txt")

geneLength$Geneid<-stringr::str_remove(geneLength$Geneid,"\\.\\d+")

TPMs<-counts%>%
  left_join(geneLength,by="Geneid")%>%
  select(1:2,130,3:129)%>%
  tidyr::gather(sample,counts,-Geneid,-Symbol,-Length)%>%
  mutate(cpm=(counts*1000/Length))%>%
  group_by(sample)%>%
  mutate(totals=sum(cpm))%>%
  ungroup()%>%
  mutate(TPM=cpm*1000000/totals)
  
newTPM<-TPMs%>%
  select(Symbol,sample,TPM)%>%
  tidyr::pivot_wider(names_from = "sample",values_from = "TPM")

#write_tsv(newTPM,"RNAseq/Protein_coding.tpms.txt")
###fix samples#############################################################
case<-newTPM%>%
  select(1:72)
ctrl<-newTPM%>%
  select(73:128)
case<-case%>%
  rename("019N"="019T",
         "020N"="020T",
         "039N"="039T")
ctrl<-ctrl%>%
  rename("019T"="019N",
         "020T"="020N",
         "039T"="039N")

TPMs2<-bind_cols(case,ctrl)
TPMs3<-TPMs2%>%
  select(1:13,89,90,16:32,105,34:88,14,15,91:104,33,106:128)
write_tsv(TPMs3,"RNAseq/Protein_coding.tpms.txt")

#######################################################################################
#########################################################################################