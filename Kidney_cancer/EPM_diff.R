###first cal EPM
infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2
normals<-paste0(need_samples,"N")
cases<-paste0(need_samples,"T")
buchong<-c("141Bone","141M","142M","142Vein","142PBMC")

getNums<-function(x){
  df<-read_tsv(stringr::str_glue("filter/cKca_{x}_circle_site.filter.tsv"))
  nrow(df)
}

fin<-tibble(eccnumbers=as.character(sapply(c(normals,cases),getNums)),sample=c(normals,cases))
fin<-fin%>%
  mutate(group=if_else(sample %in% cases,"Tumour","Normal"))

mappings<-openxlsx::read.xlsx("total_mappings.xlsx")
mappings$sample<-stringr::str_remove(mappings$sample,"cKca_")
fin<-fin%>%
  left_join(mappings,by="sample")
fin$totalmappings<-as.numeric(fin$totalmappings)
fin$eccnumbers<-as.numeric(fin$eccnumbers)
fin<-fin%>%
  mutate(EPM=eccnumbers/totalmappings*10^6)
openxlsx::write.xlsx(fin,"Kidney_circle_epms.xlsx")

### compare_with_UBC
RCC<-fin%>%
  dplyr::select(EPM,group)%>%
  mutate(group=paste0("RCC-",group))

fin2<-openxlsx::read.xlsx("../../Bladder/Bladder_circle_numbers.xlsx")

UBC<-fin2%>%
  dplyr::select(EPM,group)%>%
  mutate(group=paste0("UBC-",group))

merges<-rbind(RCC,UBC)

tongji<-merges%>%rstatix::t_test(EPM~group,p.adjust.method = "BH")

RCC%>%
  rstatix::t_test(EPM~group,paired = T)
UBC%>%
  rstatix::t_test(EPM~group,paired = T)

save(fin,file = "ecc_element_numbers.Rdata")


### plot EPM
ggplot(plotdata,aes(x=label))+
  geom_point(aes(y=EPM),color="#34578B")+
  geom_smooth(aes(x=label,y=var),method = "loess",color="#9A3735",se = F)+
  geom_hline(yintercept = avg,linetype="dashed")+
  theme_classic(base_size = 14)+
  theme(panel.grid.major.y = element_line(color="grey"))

ggsave("epm.plot.pdf",width = 5.57,height =3.06 )

## EPM > average 
naxie <- plotdata%>%filter(EPM>avg)%>%pull(sample)%>%unique()
naxie <-stringr::str_remove(naxie,"T")
# load anno
infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
infos$Grade<-as.character(infos$Grade)

new_info <-infos%>%
  filter(Sample2%in%naxie)%>%
  select(Sample2,Histology,Grade)
new_info$Sample2<-as.character(new_info$Sample2)
# load MAF

sigGenes<-c("MSI-H","POLE","PMS2","PMS1","MSH6","MSH3","MSH2","MLH3","MLH1")
Variants <- c("Missense_Mutation",
              "Nonsense_Mutation", 
              "Nonstop_Mutation",
              "Frame_Shift_Ins",
              "Frame_Shift_Del",
              "In_Frame_Del",
              "In_Frame_Ins",
              "Splice_Site",
              "Translation_Start_Site")

df2<-read_tsv("../WGS/KCA.WGS.maf")
df2<-df2%>%
  filter(Variant_Classification%in%Variants)
df2<-df2%>%
  filter(Tumor_Sample_Barcode!="CCGA-RCC-141M")%>%
  filter(Tumor_Sample_Barcode!="CCGA-RCC-142M")
df2$Tumor_Sample_Barcode<-stringr::str_remove(df2$Tumor_Sample_Barcode,"CCGA-RCC-")
df2$Tumor_Sample_Barcode<-stringr::str_remove(df2$Tumor_Sample_Barcode,"^0+")
df2$Tumor_Sample_Barcode<-stringr::str_remove(df2$Tumor_Sample_Barcode,"T")
#write_tsv(df2,"../WGS/KCA.WGS.modify.maf")
df2<-df2%>%
  filter(Tumor_Sample_Barcode%in%naxie)
df2<-df2%>%
  filter(Hugo_Symbol%in%sigGenes)
maf <- read.maf(df2)
#oncoplot(maf,top = 30,writeMatrix = T)

plotData<-read.table("onco_matrix.txt",sep="\t",check.names = F)
otherSamples<-setdiff(naxie,colnames(plotData))
plotData_buchong<-matrix(rep("",66),nrow =2 )
colnames(plotData_buchong)<-otherSamples
rownames(plotData_buchong)<-rownames(plotData)
plotData<-cbind(plotData,plotData_buchong)
plotData<-plotData[,new_info$Sample2]

col = c("Missense_Mutation" = "#a6bddb", 
        "Nonsense_Mutation" = "#fdc086", 
        "Frame_Shift_Ins"="#e78ac3",
        "Frame_Shift_Del"="#a6d854",
        "In_Frame_Del"="#96CED5",
        "In_Frame_Ins"="#CF4A31",
        "Nonstop_Mutation"="#484125",
        "Splice_Site"="#beaed4",
        "Translation_Start_Site"="#F88D51",
        "Multi_Hit"="#DC494C")

alter_fun = list(
  background = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = "#f0f0f0"),   
  Missense_Mutation = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Missense_Mutation"]),
  Nonsense_Mutation = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Nonsense_Mutation"]),
  Frame_Shift_Ins = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Frame_Shift_Ins"]),
  Frame_Shift_Del = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Frame_Shift_Del"]),
  In_Frame_Del = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["In_Frame_Del"]),
  In_Frame_Ins = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["In_Frame_Ins"]),
  Nonstop_Mutation = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Nonstop_Mutation"]),
  Splice_Site = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Splice_Site"]),
  Translation_Start_Site = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Translation_Start_Site"]),
  Multi_Hit = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Multi_Hit"])
)

ht1<-oncoPrint(plotData,
          alter_fun = alter_fun, 
          col = col,
          alter_fun_is_vectorized = FALSE,
          show_pct = F,
          row_names_side = "left",
          row_names_gp = gpar(fontface="italic",fontsize=8),
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(axis = T,axis_param = list(side="top")),
                                           show_annotation_name=F),
          top_annotation=HeatmapAnnotation(Histology=new_info$Histology,
                                           Grade=new_info$Grade,
                                           simple_anno_size = unit(0.3, "cm"),
                                           col = list(Histology=c("ccRCC"="#EBABAC","CDC"="#7662A5","chRCC"="#A68575","pRCC"="#73A99B","sRCC"="#747070"),
                                                      Grade=c("0"="#747070","1"="#F5F6F7","2"="#E1E2D2","3"="#8A8A8A","4"="#EBD4E4")))
                                                      
)

pdf("test_info.pdf",width = 6.03,height = 1)
draw(ht1)
dev.off()
