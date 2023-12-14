data<-openxlsx::read.xlsx("mapping_ratio.xlsx")
data<-as_tibble(data)

data<-data%>%
  mutate(category=case_when(stringr::str_ends(sample,"T") ~ "Tumor",
                            stringr::str_ends(sample,"N") ~ "Normal",
                            stringr::str_starts(sample,"C-") ~ "cellline",
                            TRUE ~ "others"))

infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2
normals<-paste0("cKca_",need_samples,"N")
cases<-paste0("cKca_",need_samples,"T")
cl <-c("C-769P","C-CaKi","C-ACHV","C-HK-2")
data<-data%>%
  filter(sample%in%c(cases,normals,cl))

plotData<-data%>%
  select(sample,category,ratio_chr,ratio_mt,retio_unmap)%>%
  dplyr::rename(chromosome=ratio_chr,chrMT=ratio_mt,others=retio_unmap)%>%
  tidyr::gather(type,value,-sample,-category)

plotData$category<-factor(plotData$category,levels = c("Tumor","Normal","cellline"))
plotData$type<-factor(plotData$type,levels = c("others","chrMT","chromosome"))
plotData<-plotData%>%
  arrange(category,sample)
plotData$sample<-factor(plotData$sample,levels = unique(plotData$sample))

ggplot(plotData,aes(x=sample,y=value,fill=type))+
  geom_col()+
  scale_fill_manual(values = c(chromosome="#9FBED9",others="#DC8F8E",chrMT="#E8AF8D"))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(size = 6,angle = 90,hjust = 1,vjust = 1))

ggsave("test.pdf",width =15.52 ,height =2.13 )
