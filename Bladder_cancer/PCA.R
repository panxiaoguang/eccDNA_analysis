eccAbundance<-readRDS("gene_drops.RDS")
dffs<-readRDS("gene_drops.diff.RDS")
logfc<-eccAbundance%>%
  mutate(mean1=rowMeans(across(2:121)),
         mean2=rowMeans(across(122:241)))%>%
  select(gene,mean1,mean2)%>%
  mutate(logfc=log2(mean1)-log2(mean2))
dffs<-dffs%>%left_join(logfc,by="gene")

diffGenes<-dffs%>%filter(p.adj<0.01,abs(logfc)>0.5)

library("FactoMineR")
library("factoextra")

data<-eccAbundance%>%
  tibble::column_to_rownames(var="gene")%>%
  as.matrix()
data2<-data[diffGenes$gene,]
data3<-t(data2)
data3<-as.data.frame(data3)

res.pca <- PCA(data3, graph = FALSE)

fviz_pca_ind(res.pca,
             geom.ind = "point" , # show points only (nbut not "text")
             col.ind = c(rep("Tumor",120),rep("NAT",120)), # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = F, # Concentration ellipses
             legend.title = "Groups",
             repel = TRUE
)
