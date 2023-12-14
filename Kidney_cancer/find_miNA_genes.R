## read miRNA from miRNAbaseV22
allmirna<-read_tsv("miRNA_primary_transcript.bed",col_names = F)
names(allmirna)<-c("seqnames","start","end","strand","id","name")
mirna<-makeGRangesFromDataFrame(allmirna,keep.extra.columns = T)
rm(allmirna)
find_intact<-function(x){
  print(x)
  ## readin ecc
  df<-read_tsv(stringr::str_glue("filter_bed/cKca_{x}.ecc.bed"),col_names = F)
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

infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2
cases<-paste0(need_samples,"T")
normals<-paste0(need_samples,"N")

fin<-do.call('rbind',lapply(c(cases,normals),find_intact))
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
left<-miRNA_mat[,1:68]
right<-miRNA_mat[,69:243]
middle<-matrix(rep(0,nrow(left)),ncol = 1)
colnames(middle)<-"82T"
rownames(middle)<-rownames(left)
fin<-cbind(left,middle,right)

### calculate numbers
case_mat<-fin[,1:122]
normal_mat<-fin[,123:244]

case_numbers<-apply(case_mat,1,function(x){sum(x==1)})
control_numbers<-apply(normal_mat,1,function(x){sum(x==1)})

stat_table<-tibble(miRNA=rownames(fin),Casein=case_numbers,Casenon=122-case_numbers,Normalin=control_numbers,Normalnon=122-control_numbers)

### fisher test

haha<-stat_table%>%
  rowwise()%>%
  mutate(pvalues=fisher.test(matrix(c(Casein,Casenon,Normalin,Normalnon),
                                    nrow = 2,
                                    byrow = T),
                             alternative = "two.sided")$p)
haha<-haha%>%arrange(pvalues,desc(Casein))

topGenes<-haha[1:30,]$miRNA

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
        column_split = c(c(rep("RCC",122),rep("NAT",122))),
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
        bottom_annotation = HeatmapAnnotation(group=c(rep("RCC",122),rep("NAT",122)),
                                              simple_anno_size = unit(0.3, "cm"),
                                              col=list(group=(c("RCC"="#9A3735","NAT"="#34578B")))),
        row_names_gp = gpar(fontface="italic",cex=0.7)
        
)

pdf("ecc-miRNA.pdf",width = 7.42,height = 5.24)
draw(ht)
dev.off()

### usefull
find_intact_ecc<-function(x,mirna){
  print(x)
  ## readin ecc
  df<-read_tsv(stringr::str_glue("filter_bed/cKca_{x}.ecc.bed"),col_names = F)
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
  intact
}

find_intact_ecc("18T",mirna[mirna$id=="MI0003160"])



mirna[mirna$id=="MI0000114"]



## filter by pvalue<0.01
cutData<-haha%>%
  filter(pvalues<0.01)

## load names

miRNA_name<-read_tsv("miRNA_duiying.txt",col_names = F)
names(miRNA_name)<-c("miRNA","name")

diffTable<-tibble(miRNA=cutData$miRNA)%>%
  left_join(miRNA_name,by="miRNA")

write_tsv(diffTable%>%select(name),"diff_miRNA.names.tsv")

tmp<-cutData%>%filter(Normalin==0)

diffTable<-tibble(miRNA=tmp$miRNA)%>%
  left_join(miRNA_name,by="miRNA")

write_tsv(diffTable%>%select(name),"diff_miRNA.names.tsv")


saveRDS(fin,"ecc_miRNA_mat.RDS")
saveRDS(haha,"ecc_miRNA_dff.RDS")


fin<-readRDS("ecc_miRNA_mat.RDS")
haha<-readRDS("ecc_miRNA_dff.RDS")

wocao$qvalues<-p.adjust(wocao$pvalues,method = "BH")
wocao2<-wocao%>%filter(qvalues<0.1)
wocao3<-wocao2%>%arrange(pvalues)%>%.[1:25,]

wocao3<-wocao3%>%
  mutate(log2=if_else(log2==Inf,5,log2))

wocao3<-wocao3%>%arrange(desc(log2))
wocao3$miRNA<-factor(wocao3$miRNA,levels = wocao3$miRNA)
p1<-ggplot(wocao3,aes(x=ratio1,y=miRNA))+
  geom_col(fill="#9a3735",width = 0.9)+
  scale_x_continuous(limits = c(0,0.4))+
  theme_pubr()
ggsave("p1.pdf",width = 3.27,height = 4.62)
ggplot(wocao3,aes(x=ratio2,y=miRNA))+
  geom_col(fill="#34578b",width = 0.9)+
  scale_x_continuous(limits = c(0,0.4))+
  theme_pubr()
ggsave("p2.pdf",width = 3.27,height = 4.62)

wocao4$qvalues<-round(wocao4$qvalues,digits = 3)
p2<-ggtexttable(wocao4, rows = NULL, theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)

ggsave("p3.pdf",width = 4.75,height = 6.86)

ggarrange(p1,p2)

wocao3%>%left_join(miRNA_name,by="miRNA")

topmiRNAs<-wocao3%>%
  select(miRNA)

topmiRNAs<-topmiRNAs%>%left_join(miRNA_name,by="miRNA")

write_tsv(topmiRNAs,"top_miRNA_names.txt")


sam1<-names(fin["MI0003160",][fin["MI0003160",]==1])
sam2<-names(fin["MI0003135",][fin["MI0003135",]==1])
sam3<-names(fin["MI0000114",][fin["MI0000114",]==1])
sam4<-names(fin["MI0000238",][fin["MI0000238",]==1])
sam5<-names(fin["MI0003182",][fin["MI0003182",]==1])


needMi<-c("MI0003160","MI0003135","MI0000114","MI0000238","MI0003182")


haha<-lapply(table$sample,function(x){find_intact_ecc(x,mirna[mirna$id%in%needMi])})


