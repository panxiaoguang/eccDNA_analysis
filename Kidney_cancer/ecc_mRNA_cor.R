fin<-readRDS("../CircleSeq/gene_drops_122.RDS")
fin<-fin%>%
  dplyr::select(1:123)%>%
  tidyr::gather(sample,value,-gene)%>%
  dplyr::rename(Geneid=gene,EPM=value)%>%
  tidyr::replace_na(replace = list(EPM=0))
TPM<-read_tsv("count_TPM.txt")
names(TPM)[1]<-"Geneid"
goodGenes<-TPM%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  mutate(type=if_else(TPM>3,"good","bad"))%>%
  group_by(Geneid,type)%>%
  summarise(count=n())%>%
  tidyr::pivot_wider(names_from = "type",values_from = count,values_fill = 0)%>%
  filter(good==230)%>%
  pull(Geneid)

TPMs<-TPM%>%
  filter(Geneid%in%goodGenes)%>%
  dplyr::select(1:116)%>%
  tidyr::gather(sample,TPM,-Geneid)

fin$sample<-stringr::str_remove(fin$sample,"T")
TPMs$sample<-as.character(as.numeric(stringr::str_extract(TPMs$sample,"\\d+")))
fin<-inner_join(fin,TPMs,by=c("Geneid","sample"))
rst<-fin%>%
  mutate(logT=log2(TPM+1),
         logE=log2(EPM+1))%>%
  group_by(Geneid)%>%
  rstatix::cor_test(vars = "logE",vars2 = "logT",method = "spearman")

save(rst,file="spearCor.tumor.Rdata")

ggplot(rst,aes(x=cor,y=after_stat(count)))+
  geom_histogram(fill="#BD5251",color="black",bins = 100,size=0.2)+
  scale_x_continuous(limits = c(-0.6,0.6))+
  theme(legend.position = "none")+
  theme_pubr()+ylab("Gene count")+
  geom_vline(xintercept = 0.037,linetype="dashed",color="black")+
  xlab("Correlation")
ggsave("ecc_RNA_cor_tumor.pdf",width = 5.88,height = 2.59)

# for normal --------------------------------------------------------------

fin<-readRDS("../CircleSeq/gene_drops_122.RDS")
fin<-fin%>%
  dplyr::select(1,124:245)%>%
  tidyr::gather(sample,value,-gene)%>%
  dplyr::rename(Geneid=gene,EPM=value)%>%
  tidyr::replace_na(replace = list(EPM=0))
TPM<-read_tsv("count_TPM.txt")
names(TPM)[1]<-"Geneid"
goodGenes<-TPM%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  mutate(type=if_else(TPM>3,"good","bad"))%>%
  group_by(Geneid,type)%>%
  summarise(count=n())%>%
  tidyr::pivot_wider(names_from = "type",values_from = count,values_fill = 0)%>%
  filter(good==230)%>%
  pull(Geneid)

TPMs<-TPM%>%
  filter(Geneid%in%goodGenes)%>%
  dplyr::select(1,117:231)%>%
  tidyr::gather(sample,TPM,-Geneid)

fin$sample<-stringr::str_remove(fin$sample,"N")
TPMs$sample<-as.character(as.numeric(stringr::str_extract(TPMs$sample,"\\d+")))
fin<-inner_join(fin,TPMs,by=c("Geneid","sample"))
rst2<-fin%>%
  mutate(logT=log2(TPM+1),
         logE=log2(EPM+1))%>%
  group_by(Geneid)%>%
  rstatix::cor_test(vars = "logE",vars2 = "logT",method = "spearman")
save(rst2,file="spearCor.normal.Rdata")

ggplot(rst2,aes(x=cor))+
  geom_histogram(fill="#4D4FAE",color="black",bins = 100,size=0.2)+
  scale_x_continuous(limits = c(-0.6,0.6))+
  theme(legend.position = "none")+
  theme_pubr()+ylab("Gene count")+
  geom_vline(xintercept = 0.026,linetype="dashed",color="black")+
  xlab("Correlation")
ggsave("ecc_RNA_cor_normal.pdf",width = 5.88,height = 2.59)


p1<-ggplot(rst,aes(x=cor))+geom_density(n=80)

data2<-layer_data(p1,1)


ggplot(rst,aes(x=cor))+
  geom_bar(aes(x=x,y=density),
           color="black",
           fill="white",
           linewidth=0.2,
           stat = "identity",
           position = position_dodge(0),data=data2)+
  geom_density(n=80,color="#BD5251")+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values = c("Tumor"="#DA60DE","Normal"="#ECBF5D"))+
  theme_pubr()

ggsave("ecc_RNA_cor_tumor.pdf",width = 5.88,height = 2.59)

p2<-ggplot(rst2,aes(x=cor))+geom_density(n=80)

data3<-layer_data(p2,1)
ggplot(rst2,aes(x=cor))+
  geom_bar(aes(x=x,y=density),
           color="black",
           fill="white",
           linewidth=0.2,
           stat = "identity",
           position = position_dodge(0),data=data3)+
  geom_density(n=80,color="#4D4FAE")+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values = c("Tumor"="#DA60DE","Normal"="#ECBF5D"))+
  theme_pubr()
ggsave("ecc_RNA_cor_normal.pdf",width = 5.88,height = 2.59)


part1<-rst%>%select(Geneid,cor)
part2<-rst2%>%select(Geneid,cor)
part1$type<-"Tumor"
part2$type<-"NAT"
alls<-inner_join(part1,part2,by="Geneid")
nima<-alls%>%select(Geneid,cor.x,cor.y)%>%setNames(c("Geneid","Tumor","Normal"))%>%tidyr::gather(value,Group,-Geneid)
ggplot(nima,aes(x=value,y=Group))+
  geom_boxplot(aes(fill=value),outlier.shape = NA,width=0.5)+
  scale_fill_manual(values = c("Tumor"="#9E3735","Normal"="#48436D"))+
  theme_pubr()

ggsave("cor.diff.pdf",width =3.04 ,height =4.19 )

nima%>%rstatix::t_test(Group~value,paired = F)
