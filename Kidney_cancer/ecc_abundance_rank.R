# Figure3G-----plot -------------------------------------------------------
EPMs<-readRDS("gene_drops_122.RDS")
newEPM<-EPMs%>%tidyr::gather(samples,counts,-gene)%>%mutate(gp=if_else(stringr::str_ends(samples,"T"),"Tumor","Normal"))
newEPM<-newEPM%>%
  tidyr::replace_na(replace = list(counts=0))%>%
  mutate(logc=log2(counts+1))%>%
  group_by(samples)%>%
  mutate(rank=dense_rank(desc(logc)))%>%
  ungroup()

ggplot(newEPM,aes(x=rank,y=logc))+
  geom_line(aes(group=samples,color=gp))+
  scale_color_manual(values = c("Tumor"="#C8413B","Normal"="#364BBA"))+
  xlab("Genes ranks")+
  ylab("eccDNA density")+
  facet_wrap(.~gp)+
  theme_pubr()+
  theme(strip.background = element_blank())

ggsave("ecc_abundance_rank.pdf",width = 5.16,height = 4.4)

normal<-newEPM%>%filter(gp=="Normal")
tm<-newEPM%>%filter(gp=="Tumor")

tm%>%
  arrange(rank)%>%
  filter(rank<=10)%>%
  count(gene)%>%
  arrange(desc(n))
### H3C7;C6orf226;H1-3;H2BC17;H3C1
normal%>%
  arrange(rank)%>%
  filter(rank<=10)%>%
  count(gene)%>%
  arrange(desc(n))
### OR2T1;H4C13;OR5AC2;SPACA4;FTMT

### sup

EPMs<-readRDS("gene_drops_122.RDS")
newEPM<-EPMs%>%tidyr::gather(samples,counts,-gene)%>%mutate(gp=if_else(stringr::str_ends(samples,"T"),"Tumor","Normal"))
newEPM<-newEPM%>%
  tidyr::replace_na(replace = list(counts=0))%>%
  mutate(logc=log2(counts+1))%>%
  group_by(samples)%>%
  mutate(rank=dense_rank(desc(logc)))%>%
  ungroup()
ggplot(newEPM,aes(x=logc,color=samples))+
  geom_density(linewidth=0.3)+
  scale_color_manual(values = c(rep("#c5362c",122),rep("#4475a7",122)))+
  theme_pubr()+
  theme(legend.position = "none",axis.title.x = element_blank())+
  ylab("Density")
ggsave("ecc-abundance-density.pdf",width = 4.5,height = 2.72)
