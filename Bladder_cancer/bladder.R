library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
## length
normals<-c("cBca_1N", "cBca_2N", "cBca_3N", "cBca_4N", "cBca_5N", "cBca_6N", "cBca_7N", "cBca_8N", "cBca_9N", "cBca_10N", "cBca_11N", "cBca_12N", "cBca_13N", "cBca_14N", "cBca_15N", "cBca_16N", "cBca_17N", "cBca_18N", "cBca_19T", "cBca_20T", "cBca_21N", "cBca_23N", "cBca_24N", "cBca_25N", "cBca_26N", "cBca_27N", "cBca_28N", "cBca_29N", "cBca_30N", "cBca_31N", "cBca_32N", "cBca_33N", "cBca_34N", "cBca_35N", "cBca_36N", "cBca_37N", "cBca_38N", "cBca_39T", "cBca_40N", "cBca_41N", "cBca_42N", "cBca_43N", "cBca_44N", "cBca_45N", "cBca_46N", "cBca_47N", "cBca_48N", "cBca_50N", "cBca_51N", "cBca_52N", "cBca_54N", "cBca_55N", "cBca_56N", "cBca_57N", "cBca_58N", "cBca_59N", "cBca_60N", "cBca_62N", "cBca_63N", "cBca_64N", "cBca_65N", "cBca_66N", "cBca_67N", "cBca_68N", "cBca_69N", "cBca_70N", "cBca_71N", "cBca_72N", "cBca_74N", "cBca_75N", "cBca_76N", "cBca_77N", "cBca_78N", "cBca_79N", "cBca_80N", "cBca_84N", "cBca_85N", "cBca_86N", "cBca_87N", "cBca_88N")

cases<-c("cBca_1T", "cBca_2T", "cBca_3T", "cBca_4T", "cBca_5T", "cBca_6T", "cBca_7T", "cBca_8T", "cBca_9T", "cBca_10T", "cBca_11T", "cBca_12T", "cBca_13T", "cBca_14T", "cBca_15T", "cBca_16T", "cBca_17T", "cBca_18T", "cBca_19N", "cBca_20N", "cBca_21T", "cBca_23T", "cBca_24T", "cBca_25T", "cBca_26T", "cBca_27T", "cBca_28T", "cBca_29T", "cBca_30T", "cBca_31T", "cBca_32T", "cBca_33T", "cBca_34T", "cBca_35T", "cBca_36T", "cBca_37T", "cBca_38T", "cBca_39N", "cBca_40T", "cBca_41T", "cBca_42T", "cBca_43T", "cBca_44T", "cBca_45T", "cBca_46T", "cBca_47T", "cBca_48T", "cBca_50T", "cBca_51T", "cBca_52T", "cBca_54T", "cBca_55T", "cBca_56T", "cBca_57T", "cBca_58T", "cBca_59T", "cBca_60T", "cBca_62T", "cBca_63T", "cBca_64T", "cBca_65T", "cBca_66T", "cBca_67T", "cBca_68T", "cBca_69T", "cBca_70T", "cBca_71T", "cBca_72T", "cBca_74T", "cBca_75T", "cBca_76T", "cBca_77T", "cBca_78T", "cBca_79T", "cBca_80T", "cBca_84T", "cBca_85T", "cBca_86T", "cBca_87T", "cBca_88T")


getLength<-function(x){
  df<-read_tsv(stringr::str_glue("filter_eccs/{x}_circle_site.filter.tsv"))
  tibble(sample=x,length=df$length)
}

fin<-do.call("bind_rows",lapply(c(normals,cases),getLength))

fin<-fin%>%
  mutate(group=if_else(sample %in% cases,"Tumour","Normal"))
##catogray
haha2<-fin%>%
  mutate(gps=cut(length,
                 breaks = c(-1,2000,10000,Inf),
                 labels = c("<2K","2k~10K",">10K")))%>%
  group_by(sample,gps)%>%
  summarise(number=n())%>%
  mutate(ratio=number/sum(number)*100)%>%
  ungroup()

haha2<-haha2%>%
  mutate(group=if_else(sample %in% cases,"Tumour","Normal"))

openxlsx::write.xlsx(list(haha,haha2),"bladder_length_category.xlsx")
###
plotData<-fin%>%
  filter(length<=2000)

ggdensity(plotData,
          x="length",
          y="..count..",
          color = "group",
          palette = c("#3E68B2","#A5303B"),
          alpha=0.5,
          xlab = "The length distribution of eccDNA",
          ylab = "Count")+
  theme(panel.border = element_rect(color = "black"))

ggplot(plotData,aes(x=length,y=..count..))+
  geom_density(aes(color = group))+
  theme_prism(border = T)+
  xlab("The length distribution of eccDNA")+
  ylab("Count")+
  scale_color_manual(values = c("#3E68B2","#A5303B"))

ggsave("figs/length.density_withBorder.pdf",width = 7.08,height = 2.69)

## length different
fin<-fin%>%
  mutate(length2=length/1000)
ggviolin(fin,
         x="group",
         y="length2",
         fill="group",
         palette = c("#3E68B2","#A5303B"),
         yscale = "log2",
         xlab = "",
         ylab = "Length log2(Kb)")+
  stat_compare_means()

ggsave("figs/length.violin.pdf",width = 3.16,height =3.58 )

###eccnumbers

getNums<-function(x){
  df<-read_tsv(stringr::str_glue("filter_eccs/{x}_circle_site.filter.tsv"))
  nrow(df)
}

fin<-tibble(eccnumbers=as.character(sapply(c(normals,cases),getNums)),sample=c(normals,cases))
fin<-fin%>%
  mutate(group=if_else(sample %in% cases,"Tumour","Normal"))

mappings<-read_tsv("total_mappings.tsv")
fin<-fin%>%
  left_join(mappings,by="sample")
fin$totalmappings<-as.numeric(fin$totalmappings)
fin$eccnumbers<-as.numeric(fin$eccnumbers)
fin<-fin%>%
  mutate(EPM=eccnumbers/totalmappings*10^6)

ggpaired(fin,x="group",
          y='EPM',
          color = "group",
          palette = c("#3E68B2","#A5303B"),
          width = 0.5,
          line.color = "gray", line.size = 0.4,
         xlab = "",
         ylab = "eccDNA numbers per\n million mapping reads")+
  stat_compare_means(paired = TRUE)
ggsave("figs/epm.diff.pdf",width = 3.35,height = 4.31)


df<-read_tsv("filter_eccs/cBca_59T_circle_site.filter.tsv")

### diff genes
genes<-read_tsv("../qd-ECC4/S/CircleSeq/ECC_report/FinallyData/bedFile/dbs/hg38.coding.bed",col_names = F)
names(genes)<-c("chr","Start","End","gene")
genes<-genes%>%
  mutate(length=End-Start)
get_genes<-function(x,dbs=genes){
  df<-read_tsv(stringr::str_glue("gene_anno/{x}.startAnno.bed"),col_names = F)
  df<-df%>%
    select(1:4,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    group_by(X8)%>%
    distinct(ecc,.keep_all = T)%>%
    summarise(count=n())%>%
    arrange(desc(count))%>%
    rename(gene=X8)%>%
    left_join(dbs,by="gene")%>%
    mutate(pct=count/length)%>%
    mutate(pct2=pct/sum(pct)*(10^6))%>%
    select(1,8)
  names(df)[2]<-x
  df
}
hebing<-function(x,y){full_join(x,y,by="gene")}
fin<-Reduce(hebing,lapply(c(cases,normals), get_genes))
fin<-fin%>%tidyr::gather(sample,value,-gene)%>%tidyr::replace_na(replace = list(value=0))%>%tidyr::pivot_wider(names_from = "sample",values_from = "value")
logfc<-fin%>%
  mutate(mean1=rowMeans(across(cBca_1T:cBca_88T)),
         mean2=rowMeans(across(cBca_1N:cBca_88N)))%>%
  select(gene,mean1,mean2)%>%
  mutate(logfc=log2(mean1)-log2(mean2))

dffs<-fin%>%
  tidyr::gather(sample,TPM,-gene)%>%
  tidyr::replace_na(replace = list(TPM=0))%>%
  mutate(gp=if_else(sample %in% cases,"Case","Normal"))%>%
  group_by(gene)%>%
  rstatix::pairwise_wilcox_test(TPM ~ gp,p.adjust.method = "BH")
dffs<-dffs%>%select(-mean1,-mean2,-logfc)%>%left_join(logfc,by="gene")
saveRDS(dffs,"gene_drop_diffs.RDS")
dffs<-readRDS("gene_drop_diffs.RDS")

fin<-readRDS("../Bladder/gene_drops.RDS")
dffs<-readRDS("../Bladder/gene_drop_diffs.RDS")


plotData<-fin%>%tibble::column_to_rownames(var="gene")%>%as.matrix()
nima<-dffs%>%filter(p.adj<0.01,abs(logfc)>0.5)
nima<-nima%>%
  arrange(logfc)
plt<-plotData[nima$gene,]
plt[is.na(plt)]<-0
plt<-log2(plt+1)

fin%>%filter(gene=="ADGRL3")%>%tidyr::gather(sample,value,-gene)%>%mutate(value=log2(value+1))%>%mutate(group=if_else(sample %in% cases,"Tumour","Normal"))%>%ggplot(aes(x=group,y=value,color=group))+geom_violin()
Heatmap(plt,
        show_row_names = F,
        show_column_names = F,
        cluster_rows = T,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(group=c(rep("Case",80),rep("Normal",80)),
                                           col = list(group=c("Case"="#E7B800","Normal"="#00AFBB"))),
        name = "Log normalized \neccDNA count")
write_tsv(tibble(keygenes=nima$gene),"test.genes.tsv")
col_fun<-colorRamp2(c(-2,0,2),c("#2166ac","#f7f7f7","#b2182b"))
ht<-Heatmap(plt2,
        show_row_names = F,
        show_column_names = F,
        cluster_rows = F,
        cluster_columns = F,
        #col = col_fun,
        top_annotation = HeatmapAnnotation(group=c(rep("Case",80),rep("Normal",80)),
                                           col = list(group=c("Case"="#e41a1c","Normal"="#377eb8"))),
        name = "Log normalized \neccDNA count",
        row_split = c(rep("down regular",326),rep("up regular",1019)))

keyGenes1<-nima%>%filter(logfc<0)
keyGenes2<-nima%>%filter(logfc>0)
pdf("gene_drops.pdf",width = 5.58,height = 5.83)
draw(ht)
dev.off()
write_tsv(keyGenes2,"gene_drops.genes_bot.tsv")
write_tsv(keyGenes1,"gene_drops.genes_top.tsv")
## GO enrich

ego<-read_tsv("keygenes_enrich/enrichment_results_wg_result1649823777.txt")
ggbarplot(ego,
          y="enrichmentRatio",
          x="description",
          color = "white",
          orientation = "horiz",
          fill="FDR",
          sort.val="asc",
          width = 0.9,
          xlab = "",
          ylab = "Enrichment Ratio")+
  scale_fill_gradient(low = "#084594",high = "#eff3ff")
ggsave("gene_drops.genesenich.pdf",width = 5.92,height = 4.46)

fin<-readRDS("gene_drops.RDS")
dffs<-readRDS("gene_drop_diffs.RDS")
diffGenes<-dffs%>%
  filter(p.adj<0.01,abs(logfc)>0.5)
pca_d<-fin%>%
  filter(gene%in%(diffGenes$gene))
pca_d<-pca_d%>%tibble::column_to_rownames(var="gene")%>%as.matrix()
pca_d[is.na(pca_d)]<-0
pca_d<-t(pca_d)
fin.pca<-PCA(pca_d)

plotData<-fin.pca[["ind"]][["coord"]]%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="sample")%>%
  as_tibble()%>%
  mutate(group=if_else(sample %in% cases,"Case","Normal"))


ggplot(plotData,aes(x=`Dim.1`,y=`Dim.2`))+
  stat_ellipse(aes(fill=group),geom = "polygon",alpha=0.6)+
  geom_point(aes(color=group))+
  geom_hline(yintercept = 0,linetype="dashed",color="#B0AFB0")+
  geom_vline(xintercept = 0,linetype="dashed",color="#B0AFB0")+
  scale_color_manual(values = c("Case"="#9C362F","Normal"="#393B84"))+
  xlab("PC1(8.216%)")+
  ylab("PC2(3.084%)")+
 
  scale_fill_manual(values = c("Case"="#9C362F","Normal"="#393B84"))+
  ggprism::theme_prism(border = T)

ggsave("gene_drop_PCA.pdf",width =4.85 ,height =3.33)


## element annotation

df<-read_tsv("bedFile/partGenes.anno.bed",col_names = F)
partial<-df%>%
  select(X1,X2,X3,X4)%>%
  distinct(X4,.keep_all = T)
df2<-read_tsv("bedFile/fullGenes.anno.bed",col_names = F)
fulls<-df2%>%
  select(X1,X2,X3,X4)%>%
  distinct(X4,.keep_all = T)
onlypart<-partial%>%
  filter(!X4%in%(fulls$X4))
write_tsv(onlypart,"only_partial.ecc.bed",col_names = F)
write_tsv(fulls,"fulls.ecc.bed",col_names = F)

## total=45434989;fulls=2971;part=22756690;inters=22675328
## ecDNA total=96;fulls=
dbs<-read_tsv("dbs/hg38.genetic_elements_exceptCPG_fix.bed",col_names = F)
new<-dbs%>%
  distinct(X1,X2,X3,X4)
write_tsv(new,"dbs/hg38.genetic_elements_exceptCPG_fix.bed",col_names = F)
new%>%
  tidyr::separate(X4,into = c("gene","type"),sep=":")%>%
  mutate(length=X3-X2)%>%
  group_by(type)%>%
  summarise(count=n(),
            total_l=sum(length))
flst<-Sys.glob("bedFile/x*")

getCT<-function(x){
  read_tsv(x,col_names = F)%>%
    tidyr::separate(X8,into = c("gene","type"),sep = ":")%>%
    group_by(type)%>%
    summarise(count=n())%>%
    setNames(c("type",x))
}
rst<-lapply(flst,getCT)

## total gene region and repeats
cal_total_eccs<-function(x){
  read_tsv(stringr::str_glue("filter_eccs/{x}_circle_site.filter.tsv"))%>%
    nrow()
}

eccnumbers<-tibble(sample=c(cases,normals),eccs=sapply(c(cases,normals),cal_total_eccs))

cal_ingenes<-function(x){
  read_tsv(stringr::str_glue("gene_anno/{x}.startAnno.bed"),col_names = F)%>%
    select(X4)%>%
    distinct(X4)%>%
    nrow()
}
genenumbers<-tibble(sample=c(cases,normals),eccs=sapply(c(cases,normals),cal_ingenes))
haha<-full_join(eccnumbers,genenumbers,by="sample")
haha<-haha%>%
  mutate(geneP=eccs.y/eccs.x*100)
########################################for random############################################
cal_total_eccs<-function(x){
  read_tsv(stringr::str_glue("bedFile/randomBEDs/{x}.random.bed"),col_names = F)%>%
    nrow()
}

eccnumbers<-tibble(sample=normals,eccs=sapply(normals,cal_total_eccs))

cal_ingenes<-function(x){
  read_tsv(stringr::str_glue("bedFile/randomBEDs/gene_anno/{x}.startAnno.bed"),col_names = F)%>%
    dplyr::select(X4)%>%
    distinct(X4)%>%
    nrow()
}
genenumbers<-tibble(sample=normals,eccs=sapply(normals,cal_ingenes))
haha<-full_join(eccnumbers,genenumbers,by="sample")
haha<-haha%>%
  mutate(geneP=eccs.y/eccs.x*100)
###############################################################################################
calrepeat<-function(x){
  df<-read_tsv(stringr::str_glue("repeatStat/repStats.{x}.tsv"))
  sum(df$ratio)
}
repnumbers<-tibble(sample=c(cases,normals),eccs=sapply(c(cases,normals),calrepeat))

all<-haha%>%full_join(repnumbers,by="sample")
write_tsv(all,"gene_vs_repeat.tsv")

haha<-read_tsv("gene_vs_repeat.tsv")
## repeat
dbs<-read_tsv("../CKDs/repeats/UCSC_repeatMarsker.region.bed",col_names = F)
stst<-dbs%>%
  filter(X1 !="chrM")%>%
  mutate(length=X3-X2)%>%
  group_by(X4)%>%
  summarise(total_l=sum(length))

chrmSize<-read_tsv("../甲状腺癌/old/plasma/dbs/hg38.chromo.size",col_names = F)
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
write_tsv(fin,"repeats.tsv")

## for CCR4
df<-read_tsv("bedFile/CCR4_inter_eccs.bed",col_names = F)
genes<-readRDS("gene_drops.RDS")

## gene20 anno 

df<-read_tsv("allgene20.bed",col_names = F)

df<-df%>%
  filter(stringr::str_starts(X8,"APOBE"))
write_tsv(df,"APOBE.eccDNA_list.tsv")
mt<-df%>%
  mutate(sample=stringr::str_extract(X4,"cBca_\\d+[NT]"))%>%
  mutate(length=X3-X2)%>%
  filter(length<=1000000)%>%
  select(sample,X8)%>%
  distinct(sample,X8,.keep_all = T)%>%
  mutate(num=1)%>%
  tidyr::pivot_wider(names_from = "sample",
                     values_from = "num",
                     values_fill = 0)%>%
  tibble::column_to_rownames(var="X8")%>%
  as.matrix()
tongji<-df%>%
  mutate(sample=stringr::str_extract(X4,"cBca_\\d+[NT]"))%>%
  mutate(length=X3-X2)%>%
  filter(length<=1000000)%>%
  select(sample,X8)%>%
  distinct(sample,X8,.keep_all = T)%>%
  mutate(cat=if_else(sample%in%cases,"Cancer","Normal"))%>%
  group_by(X8,cat)%>%
  summarise(num=n())%>%
  tidyr::pivot_wider(names_from = cat,values_from = num,values_fill = 0)

onlyincancer<-tongji%>%
  filter(Normal==0)%>%
  arrange(desc(Cancer))%>%
  filter(X8%in%(oncogenes))

onlyinnormal<-tongji%>%
  filter(Cancer==0)%>%
  arrange(desc(Normal))%>%
  filter(X8%in%(oncogenes$prevSymbols))

oncogenes1<-openxlsx::read.xlsx("../WGS/Cancer_driver genes_Bladder cancer.xlsx")$Gene
oncogenes2<-read_tsv("../WGS/allOnco_June2021.tsv")$prevSymbols
oncogenes<-unique(c(oncogenes1,oncogenes2))
plotData<-mt[onlyincancer[1:37,]$X8,cases]
lieshu<-apply(plotData,2,function(x) sum(x,na.rm=T))
anno_table<-openxlsx::read.xlsx("../WGS/BC_Sample infornation(2)_最新.xlsx",sheet = 6)
anno_table<-anno_table%>%
  as_tibble()%>%
  select(2,5:13)%>%
  rename(Gender=`Gender*`,Grade=`Grade#`)%>%
  mutate(sampleID=stringr::str_extract(Sample,"\\d+"))
anno_table<-anno_table%>%
  mutate(Smoke=if_else(Smoke=="yes","Yes","No"),
         MPMT=if_else(MPMT=="NO","No",MPMT))
ecDNA<-openxlsx::read.xlsx("../WGS/all_detective_gene_list.xlsx")
ecDNA_samples<-ecDNA%>%mutate(sample_name=paste0("Bca_",sample_name))%>%filter(stringr::str_detect(feature,"ecDNA"))%>%pull(sample_name)%>%unique()
anno_table<-anno_table%>%
  mutate(ecDNA=if_else(Sample%in%ecDNA_samples,"Yes","No"))
ht<-Heatmap(plotData,
        show_column_names = F,
        cluster_columns = T,
        cluster_rows = F,
        show_column_dend = F,
        show_row_dend = F,
        show_heatmap_legend = F,
        rect_gp = gpar(type="none"),
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
        bottom_annotation = HeatmapAnnotation(`NMIBC/MIBC`=anno_table$`NMIBC/MIBC`,
                                              Grade=anno_table$Grade,
                                              ecDNA=anno_table$ecDNA,
                                              col = list(`NMIBC/MIBC`=c("MIBC"="#848CF6","NMIBC"="#9B69F6"),
                                                         Grade=c("High"="#000000","Low"="#E4E3E4"),
                                                         ecDNA=c("Yes"="#000000","No"="#E4E3E4")
                                              ),
                                              simple_anno_size = unit(0.3, "cm"),
                                              annotation_name_gp = gpar(fontsize=8.5)),
        row_names_gp = gpar(fontface="italic",cex=0.7)
        
        )

pdf("heatmap_plot.pdf",width =5.84 ,height =6.57)
draw(ht)
dev.off()

##load R packages
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
## preprocess data
EPM<-read_tsv("EPMs.tsv")
EPM<-EPM%>%mutate(paired=str_extract(sample,"cBca_\\d+"))
EPM$paired<-factor(EPM$paired,levels = unique(str_sort(EPM$paired,numeric = T)))
NAT<-EPM%>%filter(group=="Normal")  
TUM<-EPM%>%filter(group=="Tumour")
NAT<-NAT%>%
  arrange(paired)%>%
  mutate(label=seq(1,80))
TUM<-TUM%>%
  arrange(paired)%>%
  mutate(label=seq(1,80))
plotData<-bind_rows(NAT,TUM)
## one line to plot figures
ggplot(plotData,aes(x=label,y=EPM))+
  geom_line(aes(group=paired),
            color="#bdbdbd",
            size=0.2)+
  geom_point(aes(color=group))+
  geom_smooth(aes(color=group,fill=group),
              linetype="longdash",
              size=0.7,
              method = "loess",
              alpha=0.5)+
  scale_color_manual(values = c("Normal"="#263272",
                                "Tumour"="#B83C3E"))+
  scale_fill_manual(values = c("Normal"="#263272",
                                "Tumour"="#B83C3E"
                               ))+
  theme_pubr()+
  xlab("Samples")+
  ylab("Number of eccDNAs per \nmillion mapped reads")

## heatmap for multiomics
TPMs<-readRDS("gene_drops.RDS")
ecDNA<-openxlsx::read.xlsx("../WGS/all_detective_gene_list.xlsx")

## heatmap
## heatmap?
mat<-read_tsv("symbol.txt")
col_info<-read_tsv("cluster.txt",col_names = F)
mat<-mat[c(-1,-2),]
plotData<-mat%>%tibble::column_to_rownames(var="id")%>%
  as.matrix()
col_info<-col_info%>%arrange(X2)
plotData<-plotData[,col_info$X1]
col_info$X2<-as.character(col_info$X2)
ht<-Heatmap(plotData,
        show_row_names = F,
        show_column_names = F,
        cluster_rows = T,
        cluster_columns = F,
        show_row_dend =F,
        row_km =3,
        top_annotation = HeatmapAnnotation(group=col_info$X2,
                                           col = list(group=c("1"="#282F66","2"="#E2A164","3"="#C74744"))),
        name = "Log normalized \neccDNA count")

shunxu<-list(rownames(plotData)[row_order(ht)$`1`],
             rownames(plotData)[row_order(ht)$`2`],
             rownames(plotData)[row_order(ht)$`3`])
openxlsx::write.xlsx(shunxu,"KM_3.genes.xlsx")
pdf("KM_3.pdf",width =5.27 ,height = 4.43)
draw(ht)
dev.off()

###############################survival
library(survival)
rt<-read.table("clinicalExp.txt",header = T,sep="\t",check.names = F)
rt$futime<-rt$futime/365
outTab<-data.frame()

for(gene in colnames(rt[,4:ncol(rt)])){
  a=rt[,gene]<=median(rt[,gene])
  cox<-coxph(Surv(futime,fustat) ~ a ,data=rt)
  coxSummary<-summary(cox)
  outTab<-rbind(outTab,cbind(gene=gene,HR=coxSummary$coefficients[,"exp(coef)"],
                             z=coxSummary$coefficients[,"z"],
                             pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

##################all tumor read coverage over gene regions
get_info<-function(x){
  df<-read_tsv(stringr::str_glue("{x}.txt"),col_names = F,skip = 2)
  df<-df%>%
    select(1,3:224)%>%
    tidyr::gather(pos,value,-X1)%>%
    mutate(pos=1:222)%>%
    dplyr::rename(sample=X1)
  df
}
allsamples<-c("cBca_1T","cBca_2T","cBca_3T","cBca_4T","cBca_5T","cBca_6T","cBca_7T","cBca_8T",
              "cBca_9T","cBca_10T","cBca_11T","cBca_12T","cBca_13T","cBca_14T","cBca_15T",
              "cBca_16T","cBca_17T","cBca_18T","cBca_19T","cBca_20T","cBca_21T","cBca_23T",
              "cBca_24T","cBca_25T","cBca_26T","cBca_27T","cBca_28T","cBca_29T","cBca_30T",
              "cBca_31T","cBca_32T","cBca_33T","cBca_34T","cBca_35T","cBca_36T","cBca_37T",
              "cBca_38T","cBca_39T","cBca_40T","cBca_41T","cBca_42T","cBca_43T","cBca_44T",
              "cBca_45T","cBca_46T","cBca_47T","cBca_48T","cBca_50T","cBca_51T","cBca_52T",
              "cBca_54T","cBca_55T","cBca_56T","cBca_57T","cBca_58T","cBca_59T","cBca_60T",
              "cBca_62T","cBca_63T","cBca_64T","cBca_65T","cBca_66T","cBca_67T","cBca_68T",
              "cBca_69T","cBca_70T","cBca_71T","cBca_72T","cBca_74T","cBca_75T","cBca_76T",
              "cBca_77T","cBca_78T","cBca_79T","cBca_80T","cBca_84T","cBca_85T","cBca_86T",
              "cBca_87T","cBca_88T")

allsamples2<-c("cBca_1N","cBca_2N","cBca_3N","cBca_4N","cBca_5N","cBca_6N","cBca_7N",
               "cBca_8N","cBca_9N","cBca_10N","cBca_11N","cBca_12N","cBca_13N",
               "cBca_14N","cBca_15N","cBca_16N","cBca_17N","cBca_18N","cBca_19N",
               "cBca_20N","cBca_21N","cBca_23N","cBca_24N","cBca_25N","cBca_26N",
               "cBca_27N","cBca_28N","cBca_29N","cBca_30N","cBca_31N","cBca_32N",
               "cBca_33N","cBca_34N","cBca_35N","cBca_36N","cBca_37N","cBca_38N",
               "cBca_39N","cBca_40N","cBca_41N","cBca_42N","cBca_43N","cBca_44N",
               "cBca_45N","cBca_46N","cBca_47N","cBca_48N","cBca_50N","cBca_51N",
               "cBca_52N","cBca_54N","cBca_55N","cBca_56N","cBca_57N","cBca_58N",
               "cBca_59N","cBca_60N","cBca_62N","cBca_63N","cBca_64N","cBca_65N",
               "cBca_66N","cBca_67N","cBca_68N","cBca_69N","cBca_70N","cBca_71N",
               "cBca_72N","cBca_74N","cBca_75N","cBca_76N","cBca_77N","cBca_78N",
               "cBca_79N","cBca_80N","cBca_84N","cBca_85N","cBca_86N","cBca_87N",
               "cBca_88N")

df1<-do.call("rbind",lapply(allsamples, get_info))
df2<-do.call("rbind",lapply(allsamples2, get_info))

fin<-df2
fin<-fin%>%group_by(sample)%>%mutate(value2=as.numeric(scale(value,center = T,scale = T)))
plotData<-fin%>%select(sample,pos,value2)%>%tidyr::pivot_wider(names_from = "pos",values_from = "value2")%>%tibble::column_to_rownames(var="sample")%>%as.matrix()
plotData2<-fin%>%group_by(pos)%>%summarise(mean=mean(value2))
ht2<-Heatmap(plotData,cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,top_annotation = HeatmapAnnotation(average = anno_lines(plotData2$mean,ylim = c(-2,2))))
ht1+ht2
########################################################
df<-openxlsx::read.xlsx("../Bladder/protein-coding.xlsx")%>%
  as_tibble()
df%>%
  mutate(type=if_else(stringr::str_ends(sample,"T"),"Tumor","Normal"))%>%
  ggplot(aes(x=Protein.ratio,y=EPM,color=type))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("Normal"="#263272",
                                "Tumor"="#B83C3E"))+
  theme_pubr()

ggsave("../Bladder/corplot.epm.pdf",width =4.89 ,height =4.45 )

df%>%
  mutate(type=if_else(stringr::str_ends(sample,"T"),"Tumor","Normal"))%>%
  group_by(type)%>%
  rstatix::cor_test(vars = "Protein.ratio",vars2 = "EPM")

###########################groupheatmap###########################

df<-read_tsv("EPM and clinical variables.tsv")
mat<-df%>%
  dplyr::select(1,2)%>%
  tibble::column_to_rownames(var="CaseID")%>%
  as.matrix()%>%
  t()
ht<-Heatmap(mat,
        cluster_rows = F,
        cluster_columns = T,
        show_column_names = F,
        top_annotation = HeatmapAnnotation(EPMgroup=df$`EPM Group`,
                                           Age=df$Age,
                                           Gender=df$Gender,
                                           MIBC=df$`NMIBC/ MIBC`,
                                           Grade=df$Grade,
                                           ecDNA=df$ecDNA_status,
                                           Msig=df$MSig,
                                           PR=df$`P/R`,
                                           MPMT=df$MPMT,
                                           chromoth=df$Chromothripsis_status,
                                           gp=gpar(col="white",lwd=0.5),
                                           col = list(Age=c("<= 65"="#B9DFFB","> 65"="#68C84D"),
                                                      EPMgroup=c("Low"="#4D4FAE","High"="#BD5251"),
                                                      Gender=c("Female"="#E93420","Male"="#316DBB"),
                                                      MIBC=c("MIBC"="#AE2417","NMIBC"="#F3F3F4"),
                                                      Grade=c("High"="#4FADEB","Low"="#F3F3F4"),
                                                      ecDNA=c("Positive"="#E5CFAB","Negative"="#A5C5CE"),
                                                      Msig=c("MSig1"="#9C362F","MSig2"="#E5CFAB","MSig3"="#A5C5CE"),
                                                      PR=c("Primary"="#B0AFB0","Relapesd"="#1B1819"),
                                                      MPMT=c("Yes"="#1B1819","No"="#B0AFB0"),
                                                      chromoth=c("Positive"="#E5CFAB","Negative"="#A5C5CE")
                                                      )
                                                       
        ),name = "Abundance")


pdf("testPDF.pdf",width = 8.2,height = 2.67)
draw(ht)
dev.off()
#############################long reads###############################################
library(readr)
library(dplyr)
samples<-c("5N", "5T", "13N", "13T", "15N", "15T", "16N", "16T", "17N", "17T", "21N", "21T", "29N", "29T", "37N", "37T", "50N", "50T")

##
ft<-function(x){
  test<-read_tsv(stringr::str_glue("../Bladder/longreads/bakups/{x}.info.txt"))
  test<-test%>%
    filter(Nfullpass>=2)%>%
    distinct(fragments,.keep_all = T)
  write_tsv(test,stringr::str_glue("../Bladder/longreads/bakups/{x}.info.ft.txt"))
}
lapply(samples, ft)
# jiaoji ------------------------------------------------------------------


cal_nums<-function(x){
  df<-read_tsv(stringr::str_glue("{x}.intersect.bed"),col_names = F)
  eccDNA_short<-nrow(df)
  eccDNA_long<-read_tsv(stringr::str_glue("{x}.eccDNA_expand.bed"),col_names = F)%>%nrow()
  eccDNA_inters<-df%>%
    filter(X5!=".")%>%
    nrow()
  tibble(sample=x,eccDNALong=eccDNA_long,eccDNAShort=eccDNA_short,eccDNAInters=eccDNA_inters)
}

fin<-do.call("bind_rows",lapply(samples, cal_nums))
#####################################################

# qianheti ----------------------------------------------------------------
get_nsegment<-function(x){
  test<-read_tsv(stringr::str_glue("{x}.info.txt"))
  test<-test%>%
    filter(Nfullpass>=2)%>%
    distinct(fragments,.keep_all = T)
  test$label<-seq(1,nrow(test))
  haha<-test%>%
    select(label,fragments)%>%
    tidyr::separate_rows(fragments,sep="\\|")%>%
    mutate(orign=sapply(stringr::str_split(fragments,":"),function(x) x[[1]]))%>%
    group_by(label)%>%
    summarise(count=n())%>%
    group_by(count)%>%
    summarise(num=n())
  haha$sample<-x
  haha
}

plotData<-do.call("bind_rows",lapply(samples, get_nsegment))
plotData<-plotData%>%tidyr::pivot_wider(names_from = count,values_from = num,values_fill = 0)
write_tsv(plotData,"nsegemnt.chinmeric.txt")
################################################

# changdu -----------------------------------------------------------------
getLength_long<-function(x){
  df<-read_tsv(stringr::str_glue("{x}.info.txt"))
  df%>%
    filter(Nfullpass>=2)%>%
    distinct(fragments,.keep_all = T)%>%
    select(seqLength)
}

chang_eccDNA<-do.call("bind_rows",lapply(samples, getLength_long))

getLength_duan<-function(x){
  df<-read_tsv(stringr::str_glue("../../filter_eccs/cBca_{x}_circle_site.filter.tsv"))
  tibble(seqLength=df$length,sample=x)
}

duan_eccDNA<-do.call("bind_rows",lapply(samples,getLength_duan))
duan_eccDNA$type<-"short_read"
chang_eccDNA$type<-"long_read"
fin<-bind_rows(duan_eccDNA,chang_eccDNA)
ggplot(fin,aes(x=seqLength))+
  geom_density(aes(fill=type),color="black",alpha=0.6)+
  theme_pubr()+
  xlab("The length distribution of eccDNA")+
  ylab("Density")+
  scale_x_continuous(limits = c(0,2000))+
  #scale_y_continuous(limits = c(0,20000))+
  scale_fill_manual(values = c("#4475A7","#C5362C"))
ggsave("length.pdf",width =5.17 ,height =2.97)
##############################plot junctions#################################################
samples<-c("5N", "5T", "13N", "13T", "15N", "15T", "16N", "16T", "17N", "17T", "21N", "21T", "29N", "29T", "37N", "37T", "50N", "50T")
need_chroms<-c(paste0("chr",seq(1,22)),"chrX","chrY")
joint_links<-read_tsv("../Bladder/longreads/bakups/13N.breakpoints.txt",col_names = F)
names(joint_links)<-c("chrm1","pos1","chrm2","pos2")
joint_links<-joint_links%>%
  filter(chrm1%in%need_chroms)%>%
  filter(chrm2%in%need_chroms)

expands_joint_links<-joint_links%>%
  mutate(start1=pos1-1,end1=pos1,start2=pos2-1,end2=pos2)%>%
  select(chrm1,start1,end1,chrm2,start2,end2)
bed1<-expands_joint_links%>%
  select(1:3)
bed2<-expands_joint_links%>%
  select(4:6)
names(bed1)<-c("chrm","start","end")
names(bed2)<-c("chrm","start","end")
pdf("test.pdf",width = 12,height =6.57)
circos.initializeWithIdeogram(plotType = c("axis", "labels"),species = "hg38")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
}, track.height = 0.05, bg.border = NA)
circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5), 
                   border = NA)
dev.off()
############################################################################################
#####################################events#######################################################
get_events<-function(x){
  df<-read_tsv(stringr::str_glue("longreads/eccDNAs/{x}.info.txt"))%>%
    mutate(frags=stringr::str_remove(fragments,"\\(+\\)"))%>%
    mutate(frags=stringr::str_remove(fragments,"\\(-\\)"))%>%
    filter(Nfullpass>=2)%>%
    group_by(frags)%>%
    summarise(count=n())%>%
    mutate(count2=if_else(count>10,">10",as.character(count)))%>%
    mutate(count2=if_else(count>10,">10",as.character(count)))%>%
    group_by(count2)%>%
    summarise(number=n())
  names(df)<-c("Events","numbers")
  df$sample<-x
  df
}

fin<-do.call("rbind",lapply(samples,get_events))

fin2<-fin%>%tidyr::pivot_wider(names_from = "sample",values_from = "numbers",values_fill = 0)
fin3<-fin%>%group_by(Events)%>%summarise(totals=sum(numbers))

cal_totals<-function(x){
    df<-read_tsv(stringr::str_glue("longreads/eccDNAs/{x}.info.txt"))%>%
      mutate(frags=stringr::str_remove(fragments,"\\(+\\)"))%>%
      mutate(frags=stringr::str_remove(fragments,"\\(-\\)"))%>%
      filter(Nfullpass>=2)%>%
      group_by(frags)%>%
      summarise(count=n())%>%
      nrow()
    df

}
totals<-tibble(sample=samples,numbers = sapply(samples, cal_totals))

############################################################################################
readGC<-function(x){
  df<-read_tsv(stringr::str_glue("bedFile/randomBEDs/GCs/{x}.random.gcContents.txt"))%>%
    tidyr::gather(type,value,-ecc)%>%
    mutate(type2=case_when(type=="self" ~ "eccDNA",
                           type=="downstream" ~ "Downstream",
                           type=="upstream" ~ "Upstream"))
  df$sample<-x
  df
}

plotData<-plotData%>%
  mutate(group=case_when(sample %in% cases ~ "Tm",
                         sample %in% normals ~ "Normal"))
plotData<-do.call("bind_rows",lapply(c(cases,normals),readGC))
plotData2<-do.call("bind_rows",lapply(normals,readGC))
##########################################GC##########################################
genes<-read_tsv("../qd-ECC4/S/ECC_report/FinallyData/bedFile/dbs/hg38.coding.bed",col_names = F)
names(genes)<-c("chr","Start","End","gene")
genes<-genes%>%
  mutate(length=End-Start)

get_genes<-function(x,dbs=genes){
  df<-read_tsv(stringr::str_glue("gene_anno/{x}.startAnno.bed"),col_names = F)
  df<-df%>%
    dplyr::select(1:4,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    group_by(X8)%>%
    distinct(ecc,.keep_all = T)%>%
    summarise(count=n())%>%
    arrange(desc(count))%>%
    rename(gene=X8)%>%
    left_join(dbs,by="gene")%>%
    mutate(pct=count/length)%>%
    mutate(pct2=pct/sum(pct)*(10^6))%>%
    dplyr::select(1,8)
  names(df)[2]<-"start_ratio"
  df
}

get_genes2<-function(x,dbs=genes){
  df<-read_tsv(stringr::str_glue("gene_anno/{x}.endAnno.bed"),col_names = F)
  df<-df%>%
    dplyr::select(1:4,8)%>%
    mutate(ecc=paste(paste(X1,X2,sep=":"),X3,sep="-"))%>%
    group_by(X8)%>%
    distinct(ecc,.keep_all = T)%>%
    summarise(count=n())%>%
    arrange(desc(count))%>%
    rename(gene=X8)%>%
    left_join(dbs,by="gene")%>%
    mutate(pct=count/length)%>%
    mutate(pct2=pct/sum(pct)*(10^6))%>%
    dplyr::select(1,8)
  names(df)[2]<-"end_ratio"
  df
}


fin<-df1%>%inner_join(df2,by="gene")

p1<-ggscatter(fin, x = "start_ratio", y = "end_ratio",
          color = "#D9914C", shape = 20, size = 3, 
          xlab = "Start breakpoints in genes",
          ylab = "End breakpoints in genes",
          add = "reg.line", 
          add.params = list(color = "#B0AFB0"), 
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)

df3<-get_genes("cBca_2T")
df4<-get_genes2("cBca_2T")
fin2<-df3%>%inner_join(df4,by="gene")

p2<-ggscatter(fin2, x = "start_ratio", y = "end_ratio",
              color = "#D9914C", shape = 20, size = 3, 
              xlab = "Start breakpoints in genes",
              ylab = "End breakpoints in genes",
              add = "reg.line", 
              add.params = list(color = "#B0AFB0"), 
              conf.int = TRUE, 
              cor.coef = TRUE, 
              cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)
ggarrange(p1,p2)

ggsave("start_vs_end.corr.pdf",width = 8,height = 3.51)

#############################################################################

abundance<-readRDS("gene_drops.RDS")
only_tm<-abundance%>%
  dplyr::select(1:81)
only_nnor<-abundance%>%
  dplyr::select(1,82:161)

wokankan2<-only_tm%>%
  tidyr::gather(sample,abundance,-gene)%>%
  tidyr::replace_na(replace = list(abundance=0))%>%
  mutate(ifo=if_else(abundance==0,"yes","no"))%>%
  group_by(gene,ifo)%>%
  summarise(count=n())

heihei2<-wokankan2%>%
  filter(ifo=="yes")%>%
  arrange(desc(count))

genes<-heihei2%>%filter(count==80)%>%pull(gene)

wocao<-heihei%>%filter(gene%in%genes)%>%mutate(pct=count/160)

wocao<-wocao%>%
  dplyr::select(gene,pct)
openxlsx::write.xlsx(wocao,"zero_genes.xlsx")
