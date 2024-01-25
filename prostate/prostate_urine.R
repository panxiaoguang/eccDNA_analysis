library(openxlsx)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
#2.	eccDNA features in tissue biopsie
#2.1 General informations
## load sample infors
sampleInfo<-read.xlsx("sampleInfo.xlsx",sheet = 1)
sampleInfo<-sampleInfo%>%
  as_tibble()
NormalInfo<-read.xlsx("sampleInfo.xlsx",sheet = 2)
NormalInfo<-NormalInfo%>%
  as_tibble()
## plasma_samples
plasma_samples<-sampleInfo%>%
  filter(Origin=="Urine")%>%
  pull(SampleName)

normal_samples<-NormalInfo%>%
  filter(Origin=="Urine")%>%
  filter(!is.na(SampleID))%>%
  pull(SampleName)
## EPM--> numbers

getNums<-function(sample){
  df<-read_tsv(stringr::str_glue("FinallyData/urine/filterData/{sample}_circle_site.filter.tsv"))
  nrow(df)
}

fin<-tibble(eccnumbers=as.character(sapply(c(plasma_samples,normal_samples),getNums)),
            sample=c(plasma_samples,normal_samples))

fin<-read_tsv("FinallyData/urine/eccNumbers.tsv")
mappings<-read_tsv("FinallyData/urine/eccmaps.tsv")

fin<-fin%>%
  left_join(mappings,by="sample")
fin$eccnumbers<-as.numeric(fin$eccnumbers)

fin<-fin%>%mutate(group=if_else(sample %in% plasma_samples,"PRAD","NAT"))

ggplot(fin,aes(x=log2(totalmappings),y=log2(eccnumbers)))+
  geom_point(aes(color=group))
ggsave("urine.mapping.png",width =5.9 ,height = 4.54)
## length
getLength<-function(sample){
  df<-read_tsv(stringr::str_glue("FinallyData/urine/filterData/{sample}_circle_site.filter.tsv"))
  tibble(sample=sample,length=df$length)
}

fin<-do.call("bind_rows",lapply(c(plasma_samples,normal_samples),getLength))

fin<-fin%>%
  mutate(group=if_else(sample %in% plasma_samples,"PRAD","NAT"))
## save length 
save(fin,file="FinallyData/urine/allurine_length_total.Rdata")

## length plot distribution

plotData2<-fin%>%
  filter(length<=1000)

ggplot(plotData2,aes(x=length,y=after_stat(count)))+
  geom_density(aes(fill=factor(group,levels = c("PRAD","NAT"))),alpha=0.6)+
  scale_fill_manual(values = c("NAT"="#5575AB","PRAD"="#C9432F"))+
  theme_pubr()

ggsave("FinallyData/urine/length.plot.pdf",width =6.56 ,height =3.42)  

## length category compareL:

haha2<-fin%>%
  mutate(gps=cut(length,
                 breaks = c(-1,250,500,1000,Inf),
                 labels = c("<250","250~500","500~1000",">1000")))%>%
  group_by(sample,gps)%>%
  summarise(number=n())%>%
  mutate(ratio=number/sum(number)*100)%>%
  ungroup()


haha2<-haha2%>%
  mutate(group=if_else(sample %in% plasma_samples,"PRAD","NAT"))

haha<-haha2%>%group_by(gps)%>%rstatix::wilcox_test(ratio~group)

## normalise
haha3<-haha2%>%
  group_by(gps,group)%>%
  summarise(meanR=mean(ratio))%>%
  reframe(group=group,
          newMean=meanR/sum(meanR))

## circle plot
## 0.107, 0.000000356, 0.464, 0.0000483
ggplot(haha3, aes(x = 3,
                  y = newMean,
                  fill = group)) +
  geom_col(width = 1.5, 
           color = 'white') + 
  facet_grid(.~gps)+
  coord_polar(theta = "y")+ 
  xlim(c(0.2, 3.8))+
  scale_fill_manual(values = c("NAT"="#5575AB","PRAD"="#C9432F")) +
  theme_void()+
  theme(
    strip.text.x = element_text(size = 14), 
    legend.title = element_text(size = 15), 
    legend.text = element_text(size = 14) 
  )+
  geom_text(aes(label = paste0(round(newMean,2),'%')),
            position = position_stack(vjust = 0.5),
            size = 4)

ggsave("FinallyData/urine/length_cat.plot.pdf",width =7.49 ,height =3.43) 

##GC content
readGC<-function(x){
  df<-read_tsv(stringr::str_glue("FinallyData/urine/GCs/{x}.gcContents.txt"))%>%
    tidyr::gather(type,value,-ecc)%>%
    mutate(type2=case_when(type=="self" ~ "eccDNA",
                           type=="downstream" ~ "Downstream",
                           type=="upstream" ~ "Upstream"))
  df$sample<-x
  df
}
fin<-do.call('bind_rows',lapply(c(plasma_samples,normal_samples),readGC))

plotData2<-fin%>%
  mutate(type=case_when(type2=="Downstream"  ~ "Downstream ecc",
                        type2=="Upstream" ~ "Upstream ecc",
                        type2=="eccDNA"& (sample %in%plasma_samples) ~ "PRAD-eccDNA",
                        type2=="eccDNA"& (sample %in%normal_samples) ~ "NAT-eccDNA"))

plotData2$type<-factor(plotData2$type,levels = rev(c("PRAD-eccDNA","NAT-eccDNA","Upstream ecc",
                                                     "Downstream ecc")))

saveRDS(plotData2,"FinallyData/urine/urine.gc.Rds")

ggplot(plotData2,aes(x=value,fill=type))+
  geom_density_ridges(aes(y=type),alpha=0.5)+
  cowplot::theme_cowplot()+
  xlab("GC Content(%)")+
  ylab("Density")+
  scale_x_continuous(labels = c("0","25","50","75","100"),breaks = c(0,0.25,0.5,0.75,1))+
  geom_vline(xintercept = 0.39,linetype="dashed")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("PRAD-eccDNA"="#C9432F",
                               "NAT-eccDNA"="#F9E49A",
                               "Upstream ecc"="#9BBDDE",
                               "Downstream ecc"="#5575AB"))

ggsave("FinallyData/urine/GC_ridges.pdf",width = 6.09,height = 4.16)

# ecc per chrom -----------------------------------------------------------

stat_chromo_num<-function(x) {
  data<-read_tsv(stringr::str_glue("FinallyData/urine/filterData/{x}_circle_site.filter.tsv"))
  names(data)[1]<-"Chromosome"
  result<-data%>%
    group_by(Chromosome)%>%
    dplyr::summarise(num=n())
  result$samples<-x
  result
}
chromos<-c(paste0("chr",seq(1,22)),"chrX","chrY")
fin<-do.call('rbind',lapply(c(plasma_samples,normal_samples),stat_chromo_num))
fin<-fin%>%
  filter(Chromosome %in% chromos)
chromosome_size<-read_tsv("FinallyData/urine/dbs/hg38.chromo.size",col_names = F)
names(chromosome_size)<-c("Chromosome","size")
need_size<-chromosome_size%>%filter(Chromosome%in%chromos)
plotData<-fin%>%left_join(need_size,by="Chromosome")%>%
  tidyr::replace_na(replace = list(num=0))%>%
  mutate(type=if_else(samples %in% plasma_samples,"PRAD","NAT"))


plotData<-plotData%>%
  mutate(countP=num/size*10^6)%>%
  group_by(samples)%>%
  mutate(countP=countP/sum(countP)*100)%>%
  ungroup()

plotData$Chromosome<-factor(plotData$Chromosome,levels = chromos)

plotData2<-plotData%>%
  filter(Chromosome!="chrY")

ggplot(plotData2,aes(x=Chromosome,y=countP,fill=type))+
  geom_boxplot(width=0.5,color="black",linewidth=0.25,outlier.size = 0.1)+
  geom_hline(yintercept = 4.3,linetype="dashed")+
  #scale_y_continuous(limits = c(1,7))+
  scale_fill_manual(values = c("NAT"="#5575AB","PRAD"="#C9432F"))+
  theme_pubr(base_size = 12,x.text.angle = 30)+
  xlab("")+
  ylab("Percent eccDNA Per Mb")
##2,6,7,8,15
stat<-plotData2%>%
  group_by(Chromosome)%>%
  rstatix::t_test(formula =countP~type)%>%
  rstatix::adjust_pvalue(p.col = "p",method = "BH")
write_tsv(stat,"FinallyData/urine/ecc_chromosome_pvalue.tsv")
ggsave("FinallyData/urine/ecc_chromosome.pdf",width = 8.61,height = 3.23)

##oncogenes
##ecc abundance
genes<-read_tsv("FinallyData/urine/dbs/hg38.coding.bed",col_names = F)
names(genes)<-c("chr","Start","End","gene")
genes<-genes%>%
  mutate(length=End-Start)
get_genes<-function(x,dbs=genes){
  df<-read_tsv(stringr::str_glue("FinallyData/urine/start_anno/{x}.startAnno.bed"),col_names = F)
  df<-df%>%
    dplyr::select(1:4,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    dplyr::group_by(X8)%>%
    distinct(ecc,.keep_all = T)%>%
    summarise(count=n())%>%
    arrange(desc(count))%>%
    dplyr::rename(gene=X8)%>%
    left_join(dbs,by="gene")%>%
    mutate(pct=count/length)%>%
    mutate(pct2=pct/sum(pct)*(10^6))%>%
    select(1,8)
  df<-df%>%
    group_by(gene)%>%
    summarise(new=max(pct2))
  names(df)[2]<-x
  df
}



hebing<-function(x,y){full_join(x,y,by="gene")}

expres_lst<-lapply(c(plasma_samples,normal_samples), get_genes)
fin<-Reduce(hebing,expres_lst)
### remove 0 value
fin<-fin%>%
  tidyr::gather(sample,value,-gene)%>%
  tidyr::replace_na(replace = list(value=0))%>%
  tidyr::pivot_wider(names_from = "sample",values_from = "value")

saveRDS(fin,"FinallyData/urine/gene_drops.RDS")

###logFC

logfc<-fin%>%
  mutate(mean2=rowMeans(across(2:42)),
         mean1=rowMeans(across(43:101)))%>%
  select(gene,mean1,mean2)%>%
  mutate(logfc=log2(mean1)-log2(mean2))

#fin<-readRDS("gene_drops.RDS")
dffs<-fin%>%
  tidyr::gather(sample,TPM,-gene)%>%
  mutate(gp=if_else(sample %in% plasma_samples,"PRAD","NAT"))%>%
  group_by(gene)%>%
  rstatix::wilcox_test(TPM ~ gp)


saveRDS(dffs,"FinallyData/urine/gene_drops.diff.RDS")

dffs<-dffs%>%left_join(logfc,by="gene")
dffs<-dffs%>%rstatix::adjust_pvalue(p.col = "p",method = "BH")
plotData<-fin%>%tibble::column_to_rownames(var="gene")%>%as.matrix()
nima<-dffs%>%filter(p.adj<0.05,abs(logfc)>0.5)
#nima2<-dffs%>%filter(is.infinite(logfc))
#nima<-bind_rows(nima,nima2)
nima<-nima%>%
  arrange(logfc)

plt<-plotData[nima$gene,]
plt[is.na(plt)]<-0
plt<-log2(plt+1)
plt2<-t(apply(plt,1,function(x){scale(x)}))
colnames(plt2)<-colnames(plt)
col_fun<-colorRamp2(c(-4,0,4),c("#2166ac","#f7f7f7","#b2182b"))
ht<-Heatmap(plt2,
            col = col_fun,
            show_row_names = F,
            show_column_names = F,
            cluster_rows = F,
            cluster_columns = F,
            top_annotation = HeatmapAnnotation(group=c(rep("PRAD",41),rep("NAT",59)),
                                               col = list(group=c("PRAD"="#C9432F","NAT"="#5575AB"))),
            name = "Log normalized \neccDNA count")

pdf("FinallyData/urine/gene_drops.pdf",width = 5.69,height = 5.44)
draw(ht)
dev.off()


### GO and KEGG
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

dffs<-readRDS("FinallyData/urine/gene_drops.diff.RDS")
fin<-readRDS("FinallyData/urine/gene_drops.RDS")
topGenes<-nima%>%
  filter(logfc<0)%>%
  pull(gene)
#lowGenes<-nima%>%
#  filter(logfc>0)%>%
#  pull(gene)

backgrounds<-bitr(genes$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
topIds<-bitr(topGenes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
#lowIDs<-bitr(lowGenes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

eg_top<-enrichGO(topIds$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 universe = backgrounds$ENTREZID,
                 pvalueCutoff = 0.05,
                 readable = T)
#eg_low<-enrichGO(lowIDs$ENTREZID,
#                 OrgDb = org.Hs.eg.db,
#                 ont = "BP",
#                 pvalueCutoff = 0.05,
#                 readable = T)
ek_top<-enrichKEGG(topIds$ENTREZID,
                   organism = "hsa",
                   universe = backgrounds$ENTREZID,
                   pvalueCutoff = 0.05)
ek_top<-setReadable(ek_top,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
#ek_low<-enrichKEGG(lowIDs$ENTREZID,
#                   organism = "hsa",
#                   pvalueCutoff = 0.05)
#ek_low<-setReadable(ek_low,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

result<-list(GO_up=eg_top@result,
             #GO_DOWN=eg_low@result,
             KEGG_UP=ek_top@result)
#KEGG_DOWN=ek_low@result)

openxlsx::write.xlsx(result,"FinallyData/urine/diff_genes.enrich.xlsx")

## find miRNA

allmirna<-read_tsv("../eccDNA_analysis/dbs/miRNA_primary_transcript.bed",col_names = F)
names(allmirna)<-c("seqnames","start","end","strand","id","name")
mirna<-makeGRangesFromDataFrame(allmirna,keep.extra.columns = T)
rm(allmirna)
find_intact<-function(x){
  ## readin ecc
  df<-read_tsv(stringr::str_glue("FinallyData/urine/filter_bed/{x}.ecc.bed"),col_names = F)
  names(df)<-c("seqnames","start","end","eccDNA")
  gr<-makeGRangesFromDataFrame(df,keep.extra.columns = T)
  rm(df)
  ## find overlap with intact mirna
  hits <- findOverlaps(gr,mirna,ignore.strand=T)
  rst<-gr[queryHits(hits)]
  rst$mirna<-mirna[subjectHits(hits)]$id
  ints<-pintersect(gr[queryHits(hits)],mirna[subjectHits(hits)])
  rst$intW=width(ints)
  rst$geneC <- mirna[subjectHits(hits)]
  intact<-rst[width(rst$geneC)==rst$intW]
  if(length(intact)>0){
    return(tibble(sample=x,mirna=unique(intact$mirna)))
  }else{
    return(tibble(sample="NA",mirna="NA"))
  }
}

fin<-do.call('rbind',lapply(c(plasma_samples,normal_samples),find_intact))
fin<-fin%>%
  filter(sample!="NA")
fin<-fin%>%group_by(sample)%>%distinct(mirna,.keep_all = T)
fin$num<-1
## expand into matrix
miRNA_mat<-fin%>%tidyr::pivot_wider(names_from ="sample" ,values_from ="num")
miRNA_mat<-miRNA_mat%>%
  tibble::column_to_rownames(var="mirna")%>%
  as.matrix()
miRNA_mat[is.na(miRNA_mat)]<-0
left<-miRNA_mat[,1:41]
right<-miRNA_mat[,42:99]
middle<-matrix(rep(0,nrow(left)),ncol = 1)
colnames(middle)<-"SP204"
rownames(middle)<-rownames(left)
fin<-cbind(left,middle,right)

case_mat<-fin[,1:41]
normal_mat<-fin[,42:100]

case_numbers<-apply(case_mat,1,function(x){sum(x==1)})
control_numbers<-apply(normal_mat,1,function(x){sum(x==1)})

stat_table<-tibble(miRNA=rownames(fin),Casein=case_numbers,Casenon=41-case_numbers,Normalin=control_numbers,Normalnon=61-control_numbers)

### fisher test

haha<-stat_table%>%
  rowwise()%>%
  mutate(pvalues=fisher.test(matrix(c(Casein,Casenon,Normalin,Normalnon),
                                    nrow = 2,
                                    byrow = T),
                             alternative = "two.sided")$p)
haha<-haha%>%arrange(pvalues,desc(Casein))
haha<-haha%>%rstatix::adjust_pvalue(p.col = "pvalues",method = "BH")
saveRDS(haha,"FinallyData/urine/miRNA_statMat.RDS")
write_tsv(haha,"FinallyData/urine/miRNA_statMat.tsv")

topGenes<-haha[1:13,]$miRNA

plotData<-fin[topGenes,]
lieshu<-apply(plotData,2,function(x) sum(x,na.rm=T))

ht<-Heatmap(plotData,
            show_column_names = F,
            cluster_columns = T,
            cluster_rows = T,
            show_column_dend = F,
            show_row_dend = F,
            show_heatmap_legend = F,
            rect_gp = gpar(type="none"),
            column_split = c(c(rep("PRAD",41),rep("NAT",59))),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.rect(x = x, y = y, width = width, height = height, 
                        gp = gpar(col = "white", fill = "#CCCCCC",lwd=1.4))
              if(plotData[i,j]==1){
                grid.rect(x = x, y = y, width = width, height = height/2, 
                          gp = gpar(fill = "#2E318C",lwd=0,col="#CCCCCC"))
              }
            },
            top_annotation = HeatmapAnnotation(genes=anno_barplot(lieshu,
                                                                  border=F,
                                                                  axis=T,
                                                                  gp=gpar(fill="#2E318C",
                                                                          col="transparent"))),
            bottom_annotation = HeatmapAnnotation(group=c(rep("PRAD",41),rep("NAT",59)),
                                                  simple_anno_size = unit(0.3, "cm"),
                                                  col=list(group=(c("PRAD"="#C9432F","NAT"="#5575AB")))),
            row_names_gp = gpar(fontface="italic",cex=0.7)
            
)

pdf("FinallyData/urine/ecc-miRNA.pdf",width = 8/09,height = 4.17)
draw(ht)
dev.off()