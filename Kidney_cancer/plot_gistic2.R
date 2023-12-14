df<-read_tsv("CNV/gistic2_results/amp_genes.conf_90.txt",col_names = F)
amp_genes<-df[,2:72]%>%
  as.data.frame()%>%
  t()%>%
  as_tibble()%>%
  rowwise()%>%
  mutate(genes=paste(c_across(5:ncol(.)),collapse = ","))%>%
  select(1:4,ncol(.))%>%
  ungroup()%>%
  mutate(genes=stringr::str_remove_all(genes,",NA"))%>%
  mutate(genes=stringr::str_remove_all(genes,"\\["))%>%
  mutate(genes=stringr::str_remove_all(genes,"\\]"))%>%
  setNames(c("arm","qvalue","qvalue2","loc","genes"))%>%
  tidyr::separate_rows(genes,sep=",")
df2<-read_tsv("CNV/gistic2_results/del_genes.conf_90.txt",col_names = F)
del_genes<-df2[,2:78]%>%
  as.data.frame()%>%
  t()%>%
  as_tibble()%>%
  rowwise()%>%
  mutate(genes=paste(c_across(5:ncol(.)),collapse = ","))%>%
  select(1:4,ncol(.))%>%
  ungroup()%>%
  mutate(genes=stringr::str_remove_all(genes,",NA"))%>%
  mutate(genes=stringr::str_remove_all(genes,"\\["))%>%
  mutate(genes=stringr::str_remove_all(genes,"\\]"))%>%
  setNames(c("arm","qvalue","qvalue2","loc","genes"))%>%
  tidyr::separate_rows(genes,sep=",")
oncogenes<-openxlsx::read.xlsx("../dbs/ccRcc_CAGs.xlsx")
amp_genes_onco<-amp_genes%>%
  filter(genes%in%(oncogenes$Gene.Symbol))
del_genes_onco<-del_genes%>%
  filter(genes%in%(oncogenes$Gene.Symbol))
amp_genes_onco<-amp_genes_onco%>%
  tidyr::separate(loc,into = c("chrom","start","end"),sep = "[:-]")
del_genes_onco<-del_genes_onco%>%
  tidyr::separate(loc,into = c("chrom","start","end"),sep = "[:-]")

###################plot for gistic score###########################
chrom.size<-read_tsv("../../qd-ECC4/S/CircleSeq/ECC_report/DataBase/hg38.chromo.size",col_names = F)
names(chrom.size)<-c("chromName","chromlength")
chrom.size<-chrom.size[1:24,]
chrom.size<-chrom.size%>%
  mutate(chromName=forcats::fct_relevel(chromName,c(paste0("chr",seq(1,22)),"chrX","chrY")))%>%
  arrange(chromName)

chrom.size<-chrom.size%>%
  mutate(Chromosome=stringr::str_remove(chromName,"chr"),
         chromlengthCumsum=cumsum(chromlength),
         chromStartFrom0=c(0,cumsum(chromlength)[-24]),
         chromMiddlePosFrom0=chromStartFrom0+chromlength/2)
chrom.size<-chrom.size[1:22,]
chrom.size$Chromosome<-as.double(chrom.size$Chromosome)
chrom.size$ypos<-rep(c(1.5,1.75),11)



scores<-read_tsv("CNV/gistic2_results/scores.gistic")
scores<-scores%>%
  #mutate(qvalue=10**(0-`-log10(q-value)`))%>%
  left_join(chrom.size,by="Chromosome")%>%
  mutate(Start=Start+chromStartFrom0,
         End=End+chromStartFrom0,
         `G-score`=if_else(Type=="Del",`-log10(q-value)`,`-log10(q-value)`))
###arm regions
arms.db<-read_tsv("dbs/arm_region.bed",col_names = F)
tmp<-chrom.size%>%select(chromName,chromStartFrom0)
names(arms.db)<-c("chromName","start","end","armName")
arms.db<-arms.db%>%
  left_join(tmp,by="chromName")
arms.db<-arms.db%>%
  mutate(start=start+chromStartFrom0,
         end=end+chromStartFrom0)
rm(tmp)
arms.db<-arms.db%>%
  select(armName,start,end)
arms.db$armName<-stringr::str_to_lower(arms.db$armName)
arm_levels<-read_tsv("CNV/gistic2_results/broad_significance_results.txt")
amps<-arm_levels%>%
  select(Arm,`Amp frequency`)%>%
  setNames(c("armName","ampF"))%>%
  left_join(arms.db,by="armName")
dels<-arm_levels%>%
  select(Arm,`Del frequency`)%>%
  setNames(c("armName","delF"))%>%
  mutate(delF=0-delF)%>%
  left_join(arms.db,by="armName")

ggplot(scores,aes(x=Start,y=`G-score`))+
  geom_area(aes(group=Type,fill=factor(Type,levels = c("Del","Amp"))))+
  geom_vline(aes(xintercept=chromlengthCumsum),linetype=3,size=0.25)+
  geom_segment(aes(x=start,y=ampF,xend=end,yend=ampF),data = amps,size=0.25)+
  geom_segment(aes(x=start,y=delF,xend=end,yend=delF),data = dels,size=0.25)+
  scale_fill_manual(values = c("Del"="#4779A1","Amp"="#DA6967"))+
  scale_y_continuous(breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1))+
  scale_x_continuous(expand = c(0.01,0.01))+
  theme_classic()+
  theme(axis.text.x = element_blank())

ggsave("gistic_results.png",width = 8.03,height = 2.79,dpi = 300)
###################################################################
tmp<-chrom.size%>%select(chromName,chromStartFrom0)
names(amp_genes_onco)[4]<-"chromName"
names(del_genes_onco)[4]<-"chromName"
amp_genes_onco$start<-as.numeric(amp_genes_onco$start)
del_genes_onco$start<-as.numeric(del_genes_onco$start)
amp_genes_onco<-amp_genes_onco%>%
  left_join(tmp,by="chromName")%>%
  mutate(start=start+chromStartFrom0)
del_genes_onco<-del_genes_onco%>%
  left_join(tmp,by="chromName")%>%
  mutate(start=start+chromStartFrom0)

amp_genes_onco$qvalue<-as.numeric(amp_genes_onco$qvalue)
del_genes_onco$qvalue<-as.numeric(del_genes_onco$qvalue)



scores<-read_tsv("CNV/gistic2_results/scores.gistic")
scores<-scores%>%
  left_join(chrom.size,by="Chromosome")%>%
  mutate(Start=Start+chromStartFrom0,
         End=End+chromStartFrom0)
names(scores)[5]<-"score"

ampData<-scores%>%filter(Type=="Amp")
delData<-scores%>%filter(Type=="Del")      
ggplot(ampData,aes(x=Start,y=score))+
  geom_line(color="#DA6967")+
  geom_vline(aes(xintercept=chromlengthCumsum),linetype=3,size=0.25)+
  geom_text_repel(aes(x=start,y=-log10(qvalue),label=genes),data=amp_genes_onco,arrow = arrow(length = unit(0.02, "npc")),
                  box.padding = 1)+
  #scale_fill_manual(values = c("Del"="#4779A1","Amp"="#DA6967"))+
  #scale_y_continuous(breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1))+
  scale_x_continuous(expand = c(0.01,0.01))+
  theme_classic()+
  ylab("-log10(qvalue)")+
  theme(axis.text.x = element_blank())
ggsave("focal_amp.pdf",width = 7.38,height = 2.01)
ggplot(delData,aes(x=Start,y=score))+
  geom_line(color="#4779A1")+
  geom_vline(aes(xintercept=chromlengthCumsum),linetype=3,size=0.25)+
  geom_text_repel(aes(x=start,y=-log10(qvalue),label=genes),
                  data=del_genes_onco,
                  angle=90,
                  size=3,
                  arrow = arrow(length = unit(0.02, "npc")),
                  box.padding = 1)+
  #scale_fill_manual(values = c("Del"="#4779A1","Amp"="#DA6967"))+
  #scale_y_continuous(breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1))+
  scale_y_continuous(limits = c(0,60))+
  scale_x_continuous(expand = c(0.01,0.01))+
  ylab("-log10(qvalue)")+
  theme_classic()+
  theme(axis.text.x = element_blank())
ggsave("focal_del.pdf",width = 7.38,height = 2.01)
