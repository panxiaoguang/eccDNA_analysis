## readin genes
allGenes<-read_tsv("../../Bladder/dbs/hg38.coding.bed",col_names = F)
names(allGenes)<-c("seqnames","start","end","genes")
genes<-makeGRangesFromDataFrame(allGenes,keep.extra.columns = T)
rm(allGenes)

find_intact<-function(x){
  print(x)
  ## readin ecc
  df<-read_tsv(stringr::str_glue("filter_bed/cKca_{x}T.ecc.bed"),col_names = F)
  names(df)<-c("seqnames","start","end","eccDNA")
  gr<-makeGRangesFromDataFrame(df,keep.extra.columns = T)
  rm(df)
  ## find overlap with intact genes
  hits <- findOverlaps(gr,genes,ignore.strand=T)
  rst<-gr[queryHits(hits)]
  rst$genes<-genes[subjectHits(hits)]$genes
  ints<-pintersect(gr[queryHits(hits)],genes[subjectHits(hits)])
  rst$intW=width(ints)
  rst$geneC <- genes[subjectHits(hits)]
  intact<-rst[width(rst$geneC)==rst$intW]
  if(length(intact)>0){
    return(tibble(sample=x,genes=unique(intact$genes)))
  }else{
    return(tibble(sample="NA",genes="NA"))
  }
}

infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
need_samples<-infos$Sample2

fin<-do.call('rbind',lapply(need_samples,find_intact))
fin<-fin%>%
  filter(sample!="NA")
fin$num<-1
topGenes<-fin%>%
  group_by(genes)%>%
  summarise(count=n())%>%
  filter(count>5)%>%
  arrange(desc(count))%>%
  pull(genes)

#oncogenes<-openxlsx::read.xlsx("../dbs/ccRcc_CAGs.xlsx")
#intersect(topGenes,oncogenes$Gene.Symbol)

topGenes<-topGenes[-c(which(stringr::str_starts(topGenes,"AC")),which(stringr::str_starts(topGenes,"AL")))]
mat<-fin%>%
  tidyr::pivot_wider(names_from = "sample",values_from = "num")
others<-setdiff(infos$Sample2,colnames(mat))

mat<-mat%>%
  tibble::column_to_rownames(var="genes")%>%
  as.data.frame()
plotData<-mat[topGenes[1:25],]
buchong<-matrix(rep(NA,nrow(plotData)*length(others)),ncol=length(others))
buchong<-as.data.frame(buchong)
colnames(buchong)<-others
rownames(buchong)<-rownames(plotData)
plotData<-cbind(plotData,buchong)

plotData<-plotData[,as.character(infos$Sample2)]
##plot heatmap
lieshu<-apply(plotData,2,function(x) sum(x,na.rm=T))
plotData<-as.matrix(plotData)
plotData[is.na(plotData)]<-0

ht<-Heatmap(plotData,
        show_column_names = F,
        cluster_columns = F,
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
        bottom_annotation = HeatmapAnnotation(Grade=infos$Grade,
                                              Stage = infos$stage,
                                              Histology=infos$Histology,
                                              col = list(Histology=c("ccRCC"="#EBABAC","CDC"="#7662A5","chRCC"="#A68575","pRCC"="#73A99B","sRCC"="#747070"),
                                                        Grade=c("0"="#747070","1"="#F5F6F7","2"="#E1E2D2","3"="#8A8A8A","4"="#EBD4E4"),
                                                        Stage=c("1"="#F5F6F7","2"="#E1E2D2","3"="#8A8A8A","4"="#EBD4E4")
                                             ),
                                            simple_anno_size = unit(0.3, "cm"),
                                            annotation_name_gp = gpar(fontsize=8.5)),
        row_names_gp = gpar(fontface="italic",cex=0.7)
        
)

pdf("intact_genes_no_cluster.pdf",width = 6.76,height = 6.36)
draw(ht)
dev.off()


