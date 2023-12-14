infos<-openxlsx::read.xlsx("../WGS/Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
samples<-infos$Sample2
normals<-paste0(samples,"N")
cases<-paste0(samples,"T")

genes<-read_tsv("../../Bladder/dbs/hg38.coding.bed",col_names = F)
names(genes)<-c("chr","Start","End","gene")
genes<-genes%>%
  mutate(length=End-Start)
get_genes<-function(x,dbs=genes){
  df<-read_tsv(stringr::str_glue("start_anno/cKca_{x}.startAnno.bed"),col_names = F)
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
### remove 0 value
fin<-fin%>%
  tidyr::gather(sample,value,-gene)%>%
  tidyr::replace_na(replace = list(value=0))%>%
  tidyr::pivot_wider(names_from = "sample",values_from = "value")

saveRDS(fin,"gene_drops_122.RDS")

###logFC

logfc<-fin%>%
  mutate(mean1=rowMeans(across(2:123)),
         mean2=rowMeans(across(124:245)))%>%
  select(gene,mean1,mean2)%>%
  mutate(logfc=log2(mean1)-log2(mean2))

#fin<-readRDS("gene_drops_122.RDS")
dffs<-fin%>%
  tidyr::gather(sample,TPM,-gene)%>%
  mutate(gp=if_else(sample %in% cases,"Case","Normal"))%>%
  group_by(gene)%>%
  rstatix::pairwise_wilcox_test(TPM ~ gp,p.adjust.method = "BH")

saveRDS(dffs,"gene_drops.diff.RDS")

#dffs<-readRDS("gene_drops.diff.RDS")

dffs<-dffs%>%left_join(logfc,by="gene")

plotData<-fin%>%tibble::column_to_rownames(var="gene")%>%as.matrix()
nima<-dffs%>%filter(p.adj<0.01,abs(logfc)>0.5)
#nima2<-dffs%>%filter(is.infinite(logfc))
#nima<-bind_rows(nima,nima2)
nima<-nima%>%
  arrange(logfc)

plt<-plotData[nima$gene,]
plt[is.na(plt)]<-0
plt<-log2(plt+1)
plt2<-t(apply(plt,1,function(x){scale(x)}))
colnames(plt2)<-colnames(plt)
#col_fun<-colorRamp2(c(-2,0,2),c("#2166ac","#f7f7f7","#b2182b"))
ht<-Heatmap(plt2,
        show_row_names = F,
        show_column_names = F,
        cluster_rows = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(group=c(rep("Case",122),rep("Normal",122)),
                                           col = list(group=c("Case"="#e41a1c","Normal"="#377eb8"))),
        name = "Log normalized \neccDNA count",
        row_split = c(rep("down regular",44),rep("up regular",1384)))

pdf("gene_drops.pdf",width = 5.58,height = 5.83)
draw(ht)
dev.off()

### GO and KEGG
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
topGenes<-nima%>%
  filter(logfc>0)%>%
  pull(gene)
lowGenes<-nima%>%
  filter(logfc<0)%>%
  pull(gene)

topIds<-bitr(topGenes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
lowIDs<-bitr(lowGenes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

eg_top<-enrichGO(topIds$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 readable = T)
eg_low<-enrichGO(lowIDs$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 readable = T)
ek_top<-enrichKEGG(topIds$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = 0.05)
ek_top<-setReadable(ek_top,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
ek_low<-enrichKEGG(lowIDs$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = 0.05)
ek_low<-setReadable(ek_low,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

result<-list(GO_up=eg_top@result,
             GO_DOWN=eg_low@result,
             KEGG_UP=ek_top@result,
             KEGG_DOWN=ek_low@result)

openxlsx::write.xlsx(result,"diff_genes.enrich.xlsx")
