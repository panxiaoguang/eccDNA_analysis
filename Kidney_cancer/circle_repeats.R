infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2
normals<-paste0(need_samples,"N")
cases<-paste0(need_samples,"T")

calrepeat<-function(x){
  df<-read_tsv(stringr::str_glue("repeatStat/repStats.{x}.tsv"))
  sum(df$ratio)
}
repnumbers<-tibble(sample=c(cases,normals),eccs=sapply(c(cases,normals),calrepeat))



dbs<-read_tsv("../../CKDs/repeats/UCSC_repeatMarsker.region.bed",col_names = F)
stst<-dbs%>%
  filter(X1 !="chrM")%>%
  mutate(length=X3-X2)%>%
  group_by(X4)%>%
  summarise(total_l=sum(length))

chrmSize<-read_tsv("../../甲状腺癌/old/plasma/dbs/hg38.chromo.size",col_names = F)
chrmSize<-chrmSize%>%
  filter(X1!="chrM")%>%
  summarise(total_l=sum(X2))

##normalized repeats
stst<-stst%>%
  mutate(ratio=total_l/(chrmSize$total_l))
names(stst)[1]<-"name"
stst<-stst%>%
  select(1,3)%>%
  rename(lgthR=ratio)

getReps<-function(x,dd=stst){
  df<-read_tsv(stringr::str_glue("repeatStat/repStats.{x}.tsv"))
  df<-df%>%
    left_join(dd,by="name")%>%
    mutate(fin=ratio/lgthR)%>%
    select(name,fin)%>%
    setNames(c("name",x))
  df
}
fin<-Reduce(\(x,y){full_join(x,y,by="name")},lapply(c(cases,normals), getReps))
need_repeats<-c("SINE","LINE","LTR","DNA","Simple_repeat","Low_complexity","Satellite","RNA")
fin<-fin%>%
  filter(name%in%need_repeats)

barplotData<-fin%>%
  tidyr::gather(samples,value,-name)%>%
  mutate(gp=if_else(samples%in%cases,"Tumor","NAT"))%>%
  #mutate(logr=log10(value))%>%
  group_by(gp,name)%>%
  summarise(md=mean_sdl(value))%>%
  tidyr::unnest()

pvalue<-fin%>%
  tidyr::gather(samples,value,-name)%>%
  mutate(gp=if_else(samples%in%cases,"Tumor","NAT"))%>%
  group_by(name)%>%
  rstatix::t_test(value~gp,paired = T)


ggplot(barplotData,aes(x=name,y=y))+
  geom_col(aes(fill=gp),position = position_dodge2(width = 0.8,padding = 0.2),color="black")+
  geom_errorbar(aes(group=gp,ymin=ymin,ymax=ymax),
                position = position_dodge2(width = 0.1,padding = 0.2),
                linewidth=0.3)+
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label=p,y=8),data=pvalue,size=3,angle=30)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))

ggsave("repeats.diff.pdf",width = 5.50,height =3.80 )
