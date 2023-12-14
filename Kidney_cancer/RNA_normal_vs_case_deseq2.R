library(readr)
library(dplyr)
library(DESeq2)
df<-read_tsv("count_raw.txt")
df<-df%>%
  select(1,7:236)
## 115N vs 115T

counts<-df%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()

colData<-tibble(sample = colnames(counts),group=c(rep("Normal",115),rep("Tumor",115)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group) 
dds <- dds[rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds)
res <- results(dds)
#resOrdered <- res[order(res$pvalue), ]
diff_gene <-subset(res, padj < 0.01 & abs(log2FoldChange) > 1) 

allresult<-as.data.frame(res)%>%
  tibble::rownames_to_column(var="Geneid")%>%
  as_tibble()

allresult<-allresult%>%
  mutate(logp=-log10(padj))

allresult<-allresult%>%
  mutate(col=case_when((log2FoldChange<(-1))&(logp>2)~"downregular",
                       (log2FoldChange>1)&(logp>2)~"upregular",
                       TRUE ~ "none"))

test<-counts(dds,normalized=T)
test<-log2(test+1)
ht_data<-as.data.frame(test[rownames(diff_gene),])

ht_data2<-t(apply(ht_data,1,function(x){scale(x,center = T,scale = T)}))
colnames(ht_data2)<-colnames(ht_data)

geneNames = bitr(allresult$Geneid,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
names(geneNames)<-c("Geneid","id")
allresult<-allresult%>%inner_join(geneNames,by="Geneid")

geneList = allresult$log2FoldChange
names(geneList) = allresult$id
geneList = sort(geneList, decreasing = TRUE)

kk2<-gseKEGG(geneList     = geneList,
             organism     = 'hsa',
             pvalueCutoff = 0.05,
             verbose      = FALSE)
kk2<-setReadable(kk2,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

paixu<-kk2@result%>%arrange(desc(abs(NES)))

topResult<-paixu[1:9,]%>%
  as_tibble()%>%
  dplyr::select(Description,core_enrichment)%>%
  tidyr::separate_rows(core_enrichment)%>%
  dplyr::rename(Geneid=core_enrichment)%>%
  distinct(Geneid,.keep_all = T)
topResult<-topResult%>%
  filter(Geneid%in%rownames(ht_data2))

tmpData<-ht_data2[topResult$Geneid,]
anno_colors<-c("#D52126", "#88CCEE", "#FEE52C","#117733", "#CC61B0","#99C945", "#2F8AC4","#332288","#E68316","#661101")
anno_colors<-anno_colors[1:9]
names(anno_colors)<-unique(topResult$Description)
ht<-Heatmap(tmpData,
        show_row_names = F,
        show_column_names = F,
        show_row_dend = F,
        cluster_columns = F,
        cluster_rows = F,
        top_annotation = HeatmapAnnotation(group = c(rep("Normal",115),rep("Tumor",115)), 
                                           col = list(group = c("Tumor"="#C5362C","Normal"="#4475A7"))),
        column_split = c(rep("Normal",115),rep("Tumor",115)),
        right_annotation = rowAnnotation(pathway=topResult$Description,col=list(pathway=anno_colors)),
        name = "Zscore"
)

pdf("diff.heatmap.pdf",width =8.15 ,height = 5.42)
draw(ht)
dev.off()

saveRDS(dds,"normal_vs_case.deseq2.RDS")
saveRDS(allresult,"normal_vs_case.diff_table.RDS")
