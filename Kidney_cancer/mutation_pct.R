Variants <- c("Missense_Mutation",
              "Nonsense_Mutation", 
              "Nonstop_Mutation",
              "Frame_Shift_Ins",
              "Frame_Shift_Del",
              "In_Frame_Del",
              "In_Frame_Ins",
              "Splice_Site",
              "Translation_Start_Site")

df2<-read_tsv("KCA_Combine.maf")

df2<-df2%>%
  filter(Variant_Classification%in%Variants)
df2$Tumor_Sample_Barcode<-stringr::str_remove(df2$Tumor_Sample_Barcode,"CCGA-RCC-")
df2<-df2%>%
  filter(!Tumor_Sample_Barcode%in%c("68","70"))

trans<-df2%>%
  filter(Variant_Type=="SNP")%>%
  select(11:13)%>%
  mutate(type=paste0(Reference_Allele,">",Tumor_Seq_Allele2))%>%
  filter(type%in%c("C>T","C>G","C>A","T>A","T>C","T>G"))%>%
  group_by(type)%>%
  summarise(count=n())%>%
  mutate(pct=count/sum(count))





