library(readr)
library(dplyr)
df<-read_tsv("FinallyData/plasma/overlaps/Finally.out.bed",col_names = F)
names(df)<-c("chrom","start","end","hits","group","cluster_id")

## filter

newTb <- df%>%
  select(4:6)%>%
  dplyr::group_by(cluster_id, group) %>%
  dplyr::summarise(hits=round(mean(hits)), .groups = "drop")%>%
  tidyr::pivot_wider(names_from = "group",values_from = "hits",values_fill = 0)%>%
  mutate(Case_ot=length(plasma_samples)-Case,
         Control_ot=length(normal_sanples)-Control)

fisherTable<-newTb%>%
  rowwise()%>%
  mutate(pvalue=fisher.test(matrix(c(Case,Control,Case_ot,Control_ot),nrow=2))$p)%>%
  ungroup()

fisherTable<-fisherTable%>%
  rstatix::adjust_pvalue(p.col = "pvalue",method = "BH")

saveRDS(fisherTable,"FinallyData/plasma/overlaps/fisher_stat.RDS")
## find sig regions
sigPoint<-fisherTable%>%filter(pvalue.adj<=0.05)

## find raw regions in sig 

sigRegions<-df%>%
  filter(cluster_id %in% (sigPoint$cluster_id))

## remove super longer regions

keptID<-sigRegions%>%
  mutate(length=end-start)%>%
  select(cluster_id,length)%>%
  group_by(cluster_id)%>%
  filter(max(length)<=1000)%>%
  ungroup()%>%
  pull(cluster_id)

fina_region<- sigRegions%>%
  filter(cluster_id %in% keptID)

write_tsv(fina_region,"FinallyData/plasma/overlaps/sig.regions.tsv")
### ronghe changdu fenbu

merge_case<-read_tsv("FinallyData/plasma/overlaps/Case.merged_ecc.bed",col_names = F)
merge_case<-merge_case%>%
  mutate(length=X3-X2)
ggplot(merge_case,aes(x=length))+
  geom_density()+
  scale_x_log10()

merge_control<-read_tsv("FinallyData/plasma/overlaps/Control.merged_ecc.bed",col_names = F)
merge_control<-merge_control%>%
  mutate(length=X3-X2)
ggplot(merge_control,aes(x=length))+
  geom_density()+
  scale_x_log10()

## plot regions 
fina_region<-fina_region%>%
  filter(hits>1)

col_fun<-colorRamp2(c(8,12),c("#F0B0AF","#C9432F"))

plot_region<-fina_region%>%
  mutate(color=col_fun(hits))%>%
  select(chrom,start,end,cluster_id,color)%>%
  setNames(c("seqnames","start","end","name","color"))

pl<-makeGRangesFromDataFrame(plot_region,keep.extra.columns = T,ignore.strand = T)
seqlevelsStyle(pl)<-"NCBI"

ideogRam(organism = "human") %>%
  set_option(chromosomes = as.character(seq(1,22)),
             chrHeight = 300, chrMargin = 4,
             orientation = "vertical") %>%
  add_track(pl)

#####annotate
fina_region%>%
  filter(group=="Case")%>%
  select(1:4)%>%
  write_tsv("FinallyData/plasma/overlaps/High_in_case_region.bed",col_names = F)

fina_region%>%
  filter(group=="Control")%>%
  select(1:4)%>%
  write_tsv("FinallyData/plasma/overlaps/High_in_control_region.bed",col_names = F)
