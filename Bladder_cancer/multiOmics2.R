sourceOmics<-mat_cnv
targetOmics<-mat_ecc


#####process chrome location###


chrome_sourceOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
                        "12","13","14","15","16","17","18","19","20",
                        "21","22","X","Y")

chrome_targetOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
                        "12","13","14","15","16","17","18","19","20",
                        "21","22","X","Y")

#######Extract sub list#########
genelocate<-read_tsv("../Bladder/dbs/hg38.coding.bed",col_names = F)
genelocate<-genelocate%>%
  select(4,1,2,3)%>%
  mutate(X1=stringr::str_remove(X1,"chr"))%>%
  setNames(c("Symbol","chrom","start","end"))%>%
  as.data.frame()

genelocate_sourceOmics <- genelocate[genelocate[,2] %in%
                                       chrome_sourceOmics,]
genelocate_targetOmics <- genelocate[genelocate[,2] %in%
                                       chrome_targetOmics,]

intG <- intersect(rownames(targetOmics),genelocate_targetOmics[,1])

targetOmics <- targetOmics[intG,]

source_gene <- rownames(sourceOmics)
source_gene_locate <-intersect(unique(genelocate_sourceOmics[,1]),source_gene)
source_gene <- sourceOmics[source_gene_locate,]
genelocate_sourceOmics <- genelocate_sourceOmics[genelocate_sourceOmics[,1] %in% source_gene_locate,]
genelocate_targetOmics <- genelocate_targetOmics[genelocate_targetOmics[,1] %in% intG,]
###Calculate the correlation between cna and other omics data######
#corrArray <-calculateCorForTwoMatrices(source_gene,targetOmics,0.01)
#saveRDS(corrArray,"final.ecc2na.RDS")
#corrArray <-readRDS("final.ecc2na.RDS")
corrArray <-readRDS("RNAseq/final.cnv2ecc.RDS")
## functions#############
calculateChromLength <- function(chromLength,selectedChrom,genelocate){
  chromLength <- chromLength[chromLength[,1] %in% selectedChrom,,drop=FALSE]
  
  if(length(selectedChrom)==1){
    x <- 0
  }else{
    x <- c(0,chromLength[1:(nrow(chromLength)-1),2])
  }
  chromLength[,3] <- cumsum(as.numeric(x))
  chromLength[,4] <- cumsum(as.numeric(chromLength[,2]))
  
  genelocate <- cbind(genelocate,0,0)
  
  colnames(genelocate)[5:6] <- c("finalstart","finalend")
  
  for(i in c(1:nrow(genelocate))){
    chr <- genelocate[i,2]
    s <- genelocate[i,3]
    e <- genelocate[i,4]
    cs <- chromLength[chromLength[,1]==chr,3]
    genelocate[i,5] <- s+cs
    genelocate[i,6] <- e+cs
  }
  re <- list(chromLength=chromLength,genelocate=genelocate)
  return(re)
}
plotHeatMap <- function(corrArray,genelocate_sourceOmics,
                        chromLength_sourceOmics,genelocate_targetOmics,chromLength_targetOmics,
                        sourceOmicsName,targetOmicsName,dim=1){
  
  allChromlen_sourceOmics <-
    chromLength_sourceOmics[nrow(chromLength_sourceOmics),4]
  allChromlen_targetOmics <-
    chromLength_targetOmics[nrow(chromLength_targetOmics),4]
  
  if(dim==1){
    par(mar=c(4,4,4,0))
  }else{
    par(mar=c(0,4,4,0))
  }
  
  p <- which(corrArray!=0,arr.ind=TRUE)
  allcnagene <- rownames(corrArray)
  allovgene <- colnames(corrArray)
  
  la <- 1
  for(i in c(1:nrow(p))){
    
    cnag <- allcnagene[p[i,1]]
    ovg <- allovgene[p[i,2]]
    cnagp <- genelocate_sourceOmics[genelocate_sourceOmics[,1]==cnag,5]
    ovgp <- genelocate_targetOmics[genelocate_targetOmics[,1]==ovg,5]
    
    if(length(cnagp)==0 || length(ovgp)==0){
      next
    }
    
    cov <- corrArray[cnag,ovg]
    color <- ifelse(cov>0,"#E63126","#0932E3")
    
    if(la==1){
      if(dim==1){
        plot(cnagp,ovgp,main=paste(sourceOmicsName,"-", targetOmicsName,
                                   " correlation",sep=""), xlim=c(0,allChromlen_sourceOmics),
             ylim=c(0,allChromlen_targetOmics),xaxt="n",yaxt="n",
             frame.plot=FALSE,xlab=paste(sourceOmicsName,
                                         " chromosomal location",sep=""),ylab=paste(targetOmicsName,
                                                                                    " chromosomal location",sep=""),pch=20,col=color,cex=0.2)
        axis(side=1,at=(chromLength_sourceOmics[,4]-
                          chromLength_sourceOmics[,2]/2),
             labels=chromLength_sourceOmics[,1])
      }else{
        plot(cnagp,ovgp,main=paste(sourceOmicsName,"-", 
                                   targetOmicsName," correlation",sep=""),
             xlim=c(0,allChromlen_sourceOmics),
             ylim=c(0,allChromlen_targetOmics),xaxt="n",yaxt="n",
             frame.plot=FALSE,ylab=paste(targetOmicsName," chromosomal
location",sep=""),
             xlab="",pch=20,col=color,cex=0.2)
      }
      axis(side=2,at=(chromLength_targetOmics[,4]-
                        chromLength_targetOmics[,2]/2),
           labels=chromLength_targetOmics[,1])
      
      #abline(h=c(0,chromLength_targetOmics[,4]),v=c(0,chromLength_sourceOmics[,4]),
      #       col="gray",lty=3)
      la <- la+1
    }else{
      for(u in seq_len(length(cnagp))){
        for(v in seq_len(length(ovgp))){
          points(cnagp[u],ovgp[v],pch=20,col=color,cex=0.2)
        }
      }
    }
  }
}
plotSummaryBar <- function(corrArray,chromLength_sourceOmics,
                           genelocate_sourceOmics,sourceOmicsName){
  allChromlen <- chromLength_sourceOmics[nrow(chromLength_sourceOmics),4]
  haha<-tibble(gene=rownames(corrArray),
               numP=apply(corrArray,1,function(x){sum(x!=0)}))
  raw_gene_coord<-genelocate_sourceOmics%>%
    as_tibble()%>%
    select(1,5)%>%
    rename(gene=Symbol)
  haha<-haha%>%
    left_join(raw_gene_coord,by="gene")
  par(mar=c(4,4,0,0))
  plot(0,0,xlim=c(0,allChromlen),ylim=c(0,2000),type="n",xaxt="n", ##ylim=c(-1583,503)
       frame.plot=FALSE,xlab=paste(sourceOmicsName," chromosomal
      location",sep=""),
       ylab="Number of significant \ncorrelations")
  axis(side=1,at=(chromLength_sourceOmics[,4]-
                    chromLength_sourceOmics[,2]/2),
       labels=chromLength_sourceOmics[,1])
  axis(side=2,at=NULL,labels=TRUE)
  abline(v=c(0,chromLength_sourceOmics[,4]),col="gray",lty=3)
  points(haha$finalstart,haha$numP,cex=0.2,type="h",col="#E63126")
  #points(haha$finalstart,haha$numN,cex=0.2,type="h",col="#0932E3")
}
#####################
##Calculate the location of genes in the heatmap
chromLength <- read_tsv("../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size",col_names = F)
chromLength<-chromLength%>%
  filter(X1%in%c(paste0("chr",seq(1,22)),"chrX","chrY"))%>%
  mutate(X1=stringr::str_remove(X1,"chr"))%>%
  setNames(c("V1","V2"))
chromLength$V1<-factor(chromLength$V1,levels = c(seq(1,22),"X","Y"))
chromLength<-chromLength%>%
  arrange(V1)%>%
  as.data.frame()

re <-calculateChromLength(chromLength,chrome_sourceOmics,genelocate_sourceOmics)
genelocate_sourceOmics <- re$genelocate
chromLength_sourceOmics <- re$chromLength
re <-calculateChromLength(chromLength,chrome_targetOmics,genelocate_targetOmics)
genelocate_targetOmics <- re$genelocate
chromLength_targetOmics <- re$chromLength
###plot
png("ecc_mRNA_heatmap.png",width = 4.91,height =5.88,res=300,units = "in")
layout(matrix(c(1,2),2),heights=c(2,1))
plotHeatMap(corrArray,genelocate_sourceOmics,chromLength_sourceOmics,
            genelocate_targetOmics,chromLength_targetOmics,"eccDNA",
            "mRNA",dim=2)
plotSummaryBar(corrArray,chromLength_sourceOmics,
               genelocate_sourceOmics,"eccDNA")
dev.off()

pdf("ecc_mRNA_barplot.pdf",width = 5.38,height = 2.28)
plotSummaryBar(corrArray,chromLength_sourceOmics,
               genelocate_sourceOmics,"eccDNA")
dev.off()
## find cis
source_chrom<-genelocate_sourceOmics%>%
  select(Symbol,chrom)%>%
  as_tibble()%>%
  rename(cnvGenes=Symbol)
target_chrom<-genelocate_targetOmics%>%
  select(Symbol,chrom)%>%
  as_tibble()%>%
  rename(eccGenes=Symbol)

negs<-which(corrArray==(-1),arr.ind=T)
pos<-which(corrArray==1,arr.ind=T)
ecc_negs<-as.data.frame(negs)%>%
  as_tibble()%>%
  mutate(cnvGenes=purrr::map_chr(row,function(x) rownames(corrArray)[x]),
         eccGenes=purrr::map_chr(col,function(x) colnames(corrArray)[x]))
ecc_pos<-as.data.frame(pos)%>%
  as_tibble()%>%
  mutate(cnvGenes=purrr::map_chr(row,function(x) rownames(corrArray)[x]),
         eccGenes=purrr::map_chr(col,function(x) colnames(corrArray)[x]))
ecc_negs$cor<-(-1)
ecc_pos$cor<-1
ecc_cis<-bind_rows(ecc_negs,ecc_pos)
test2<-ecc_cis%>%
  left_join(source_chrom,by="cnvGenes")%>%
  left_join(target_chrom,by="eccGenes")%>%
  filter(chrom.x!=chrom.y)

###
neg_genes<-test2%>%
  filter(cor==(-1))%>%
  pull(eccGenes)%>%
  unique()
pos_genes<-test2%>%
  filter(cor==1)%>%
  pull(eccGenes)%>%
  unique()

allGenes<-test2%>%
  pull(eccGenes)%>%
  unique()

neg_genes_ent<-bitr(allGenes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
#pos_genes_ent<-bitr(pos_genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

eg1<-enrichGO(neg_genes_ent$ENTREZID,OrgDb =org.Hs.eg.db,ont = "BP",readable = T)
#eg2<-enrichGO(pos_genes_ent$ENTREZID,OrgDb =org.Hs.eg.db,ont = "BP",readable = T)

go_result<-eg1@result%>%
  as_tibble()%>%
  select(ID,Description,p.adjust)%>%
  mutate(lpgp=-log10(p.adjust))%>%
  arrange(desc(lpgp))%>%
  slice_head(n=20)

write_tsv(go_result,"cnv_ecc.trans.go.tsv")
