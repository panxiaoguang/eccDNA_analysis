### BAF############################################################################
get_BAF<-function(x){
  df1<-read_tsv(stringr::str_glue("BAF/circleseq/AlleleCNV/{x}T_tumourBAF.txt"))%>%
    select(1,4)%>%
    setNames(c("pos","Circle-seq"))
  df2<-read_tsv(stringr::str_glue("BAF/wgs/WGS_BAF/{x}T_tumourBAF.txt"))%>%
    select(1,4)%>%
    setNames(c("pos","WGS"))
  its<-df1%>%
    inner_join(df2,by="pos")%>%
    filter(WGS!=0)%>%
    filter(WGS!=1)%>%
    tidyr::gather(type,BAF,-pos)
  its
}

samples<-c("1","2","3","4","5","6","7",
           "8","9","10","11","12","13",
           "14","15","16","17","18","19",
           "20","21","23","24","25","26",
           "27","28","29","30","31","32",
           "33","34","35","36","37","38",
           "39","40","41","42","43","44",
           "45","46","47","48","50",
           "51","52","54","55","56","57",
           "58","59","60","62","63","64",
           "65","66","67","68","69","70",
           "71","72","74","75","76","77",
           "78","79","80","84","85","86",
           "87","88")

plotData<-do.call("bind_rows",lapply(samples, get_BAF))

plotData<-do.call("bind_rows",lapply(c("11","12"), get_BAF))

ggplot(plotData,aes(x=BAF,
               color=type,
               fill=type,
               y=..count../sum(..count..)/2))+
  geom_histogram(alpha=0.7,position = "identity")+
  scale_color_manual(values = c("Circle-seq"="#393A80","WGS"="#D1352B"))+
  scale_fill_manual(values = c("Circle-seq"="#393A80","WGS"="#D1352B"))+
  xlab("B-allele frequency")+
  ylab("Fraction of heterozygous\nSNPs")+
  theme_pubr()
ggsave("BAF.pdf",width = 4.71,height = 3.31)
##########################################oncoprint####################################################
########CCGA PART#####################################################################################
df<-openxlsx::read.xlsx("all_detective_gene_list.xlsx")
need_genes1<-openxlsx::read.xlsx("Cancer_driver genes_Bladder cancer.xlsx",sheet = 1)
need_genes1<-need_genes1%>%
  as_tibble()%>%
  filter(stringr::str_detect(Role,"oncogene"))
df1<-df%>%
  filter(gene%in%need_genes1$Gene)
mat<-df1%>%as_tibble()%>%
  select(sample_name,feature,gene)%>%
  mutate(feature=stringr::str_remove(feature,"_\\d"))%>%
  tidyr::pivot_wider(names_from = "sample_name",
                     values_from = "feature",
                     values_fn =function(x){paste(x,collapse = ";")},values_fill = "")%>%
  tibble::column_to_rownames(var="gene")%>%as.matrix()

needs<-df1%>%group_by(gene)%>%summarise(ct=n())%>%arrange(desc(ct))%>%filter(ct>2)%>%pull(gene)
needs<-needs[1:20]

mat2<-mat[needs,]
#############################################################################################
########################################TCGA PART################################################
TCGA<-openxlsx::read.xlsx("TCGA_UBC.xlsx")
TCGAbed<-TCGA%>%
  as_tibble()%>%
  mutate(info=paste(sample_barcode,
                    amplicon_index,
                    amplicon_classification,sep = "&"))%>%
  select(amplicon_intervals,info)%>%
  tidyr::separate_rows(amplicon_intervals,sep=",")%>%
  tidyr::separate(amplicon_intervals,into = c("chrom","start","end"),sep = "[:-]")%>%
  mutate(chrom=paste0("chr",chrom))

tcga_anno<-read_tsv("TCGA_data.anno.bed",col_names = F)
mat3<-tcga_anno%>%
  select(X4,X8)%>%
  tidyr::separate(X4,into = c("sample","amp","feature"),sep="&")%>%
  mutate(feature=case_when(feature == "BFB" ~ "BFB",
                           feature == "Circular" ~ "ecDNA",
                           feature == "Heavily-rearranged" ~ "Complex non-cyclic",
                           feature == "Linear" ~ "Linear amplification"))%>%
  filter(X8%in%(need_genes1$Gene))%>%
  select(sample,feature,X8)%>%
  dplyr::rename(gene=X8)%>%
  tidyr::pivot_wider(names_from = "sample",
                     values_from = "feature",
                     values_fn =function(x){paste(x,collapse = ";")},values_fill = "")%>%
  tibble::column_to_rownames(var="gene")%>%as.matrix()

mat4<-mat3[needs,]
#############################################################################################
final_data<-cbind(mat2,mat4)
nima=tibble(genes=rownames(final_data),num=as.character(apply(final_data,1,function(x){sum(x=="ecDNA")})))

col = c(BFB = "#436693", `Complex non-cyclic` = "#C58D65",
        ecDNA="#BC4137",`Linear amplification`="#A1B2C0",
        unknown="#ACABAC")
##############################################anno part#############################################
anno_table<-openxlsx::read.xlsx("TCGA-CCGA.xlsx",sheet = 2)
anno_table<-anno_table%>%
  mutate(label=stringr::str_extract(Sample,"\\d+"))
anno_table<-anno_table%>%
  filter(label%in%(colnames(mat2)))%>%
  tidyr::replace_na(replace = list("N"="No avalible","M"="No avalible"))%>%
  select(label,Survival,Gender,Age,N,M,`NMIBC/MIBC`,Grade)

anno_table2<-openxlsx::read.xlsx("TCGA-CCGA.xlsx",sheet = 1)
colnames(anno_table2)[6]<-"NMIBC/MIBC"
anno_table2<-anno_table2%>%
  filter(sample_barcode%in%(colnames(mat4)))%>%
  tidyr::replace_na(replace = list("NMIBC/MIBC"="No avalible","N"="No avalible","M"="No avalible"))%>%
  select(`sample_barcode`,Survival,Gender,Age,N,M,`NMIBC/MIBC`,Grade)%>%
  dplyr::rename(label="sample_barcode")%>%
  mutate(Grade=stringr::str_remove(Grade," "))

fin_anno<-bind_rows(anno_table,anno_table2)
fin_anno<-fin_anno%>%tibble::column_to_rownames(var="label")%>%
  as.data.frame()

fin_anno<-fin_anno[colnames(final_data),]
###############################################################################################
############################################plot oncoprint#######################################
ht<-oncoPrint(final_data,
              alter_fun = function(x, y, w, h, v) {
                n = sum(v)
                h = h*0.9
                if(n){grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h, 
                                gp = gpar(fill = col[names(which(v))],lwd=0.5,col="white"), just = "top")
                }else{grid.rect(x,y,w,h,gp = gpar(fill = "#E4E4E4",lwd=0.5,col="white"))}}, col = col,
              row_names_gp = gpar(fontface="italic",fontsize=8.5),
              column_split = c(rep("CCGA",46),rep("TCGA",68)),
              row_order = order(nima$num,decreasing = T),
              top_annotation = 
                Annotation(#cbar = anno_oncoprint_barplot(),
                Age = fin_anno$Age,
                N = fin_anno$N,
                M = fin_anno$M,
                Survival = fin_anno$Survival,
                Gender = fin_anno$Gender,
                `NMIBC/MIBC`=fin_anno$`NMIBC/MIBC`,
                Grade=fin_anno$Grade,
                simple_anno_size = unit(0.3, "cm"),
                annotation_name_gp = gpar(fontsize=8.5),
                col = list(Age=c("<=65"="#B9DFFB",">65"="#68C84D"),
                           N = c("N0"="#F3F3F4","N1"="#ABDAE4","N2"="#4B95E9","N3"="#123294","No avalible"="#747070"),
                           M = c("M0"="#F3F3F4","M1"="#F09F37","No avalible" = "#747070"),
                           Survival =c("Alive"="#F3F3F4","Death"="#010101"),
                           Gender=c("FEMALE"="#E93420","MALE"="#316DBB"),
                           `NMIBC/MIBC`=c("MIBC"="#AE2417","NMIBC"="#F3F3F4","No avalible"="#747070"),
                           Grade=c("High"="#4FADEB","Low"="#F3F3F4")
                )
              ))
pdf("oncoprint5.pdf",width = 12.55,height =4.37)
draw(ht)
dev.off()
###################################################################################################


totals_sample<-c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 21L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 30L, 31L, 32L, 33L, 34L, 35L, 36L, 37L, 38L, 39L, 40L, 41L, 42L, 43L, 44L, 45L, 46L, 47L, 48L, 50L, 51L, 52L, 54L, 55L, 56L, 57L, 58L, 59L, 60L, 62L, 63L, 64L, 65L, 66L, 67L, 68L, 69L, 70L, 71L, 72L, 74L, 75L, 76L, 77L, 78L, 79L, 80L, 84L, 85L, 86L, 87L, 88L)

####################################### dist#####################################################
df<-read_tsv("result.txt",col_names = T)
df<-df%>%
  setNames(c("chrom","start","end","ecDNA","BFB","linear","reaarange"))
df<-df%>%
  mutate(ecDNA_f=ecDNA/45,
         BFB_f=BFB/22,
         linear_f=linear/37,
         reaarange_f=reaarange/28)
chrom.size<-read_tsv("../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size",col_names = F)
names(chrom.size)<-c("chromName","chromlength")
chrom.size<-chrom.size[1:24,]
chrom.size<-chrom.size%>%
  mutate(chromName=forcats::fct_relevel(chromName,c(paste0("chr",seq(1,22)),"chrX","chrY")))%>%
  arrange(chromName)
chrom.size<-chrom.size%>%
  mutate(Chromosome=stringr::str_remove(chromName,"chr"),
         chromlengthCumsum=cumsum(chromlength),
         chromStartFrom0=c(0,cumsum(chromlength)[-24]),
         chromMiddlePosFrom0=chromStartFrom0+chromlength/2)
df<-df%>%left_join(chrom.size,by="chrom")%>%mutate(start=start+chromStartFrom0,
                                               end=end+chromStartFrom0)

plotData<-df%>%select(1:2,8:11,14)%>%tidyr::gather(type,proportion,-chrom,-start,-chromlengthCumsum)
plotData$type<-factor(plotData$type,levels = c("ecDNA_f","BFB_f","reaarange_f","linear_f"))
ggplot(plotData,aes(x=start,y=proportion))+
  geom_col(aes(color=type,fill=type))+
  theme_classic()+
  theme(axis.text.x = element_blank(),panel.spacing.x=unit(0,"cm"),
        axis.ticks.x = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.4))+
  scale_x_continuous(expand = c(0,0))+
  geom_vline(aes(xintercept=chromlengthCumsum),linetype=2)+
  facet_wrap(.~type,ncol = 1,strip.position = "right")+
  scale_fill_manual(values = c("BFB_f" = "#436693", "reaarange_f" = "#C58D65",
                               "ecDNA_f"="#BC4137","linear_f"="#A1B2C0"))+
  scale_color_manual(values = c("BFB_f" = "#436693", "reaarange_f" = "#C58D65",
                               "ecDNA_f"="#BC4137","linear_f"="#A1B2C0"))
ggsave("all_dist.pdf",width = 10.34,height = 3.56)

#####################################################################################
## gene_coverage###############################################################################
df<-read_tsv("PABPC1_geneRegion_readcoverage3.txt")
plotData<-df%>%
  mutate(coord=c(seq(100680816,100727809,by=50)[2:940],100727809))%>%
  tidyr::gather(sample,value,-coord)%>%
  mutate(group=if_else(sample %in%c("5T","16T","17T","21T","50T","84T"),"PABPC1-amplified","non-PABPC1-amplified"))

plotData$value2<-log2(plotData$value+1)

plotData3<-plotData%>%
  group_by(coord,group)%>%
  summarise(value3=list(mean_ci(value2)))%>%
  tidyr::unnest()%>%
  ungroup()

ggplot(plotData3,aes(x=coord,y=y,group=group))+
  geom_ribbon(aes(ymin = ymin, ymax = ymax),fill="#D3D2D3")+
  geom_line(aes(color=group))+
  geom_vline(xintercept = 100685816,color="grey")+
  geom_vline(xintercept = 100722809,color="grey")+
  scale_color_manual(values = c("#34327F","#CF3430"))+
  xlab("Genomic range")+
  ylab("Mean read coverage(log2)")+
  theme_pubr()

ggsave("PABPC1.plot3.pdf",width = 7.90,height = 3.70)
##############################################################################################

# transcripts -------------------------------------------------------------
### calculate TPM(https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
counts<-read_tsv("RNAseq/allSample.count.txt")
allgenes<-read_tsv("../Bladder/dbs/hg38.coding.bed",col_names = F)
counts<-counts%>%
  filter(Geneid%in%(allgenes$X4))
geneLength<-read_tsv("RNAseq/Gene_length.txt")
geneLength<-geneLength%>%
  filter(Geneid%in%(allgenes$X4))
write_tsv(counts,"Protein_coding.counts.txt")

TPMs<-counts%>%
  left_join(geneLength,by="Geneid")%>%
  mutate(across(where(is.numeric), ~((.x*1000/Length)*1000000)/sum(.x*1000/Length)))
write_tsv(TPMs,"Protein_coding.tpms.txt")
###fix samples#############################################################
TPMs<-read_tsv("Protein_coding.tpms.txt")
case<-TPMs%>%
  select(1:71)
ctrl<-TPMs%>%
  select(72:127)
case<-case%>%
  rename(Bca_19_N=Bca_19_T,
         Bca_20_N=Bca_20_T,
         Bca_39_N=Bca_39_T)
ctrl<-ctrl%>%
  rename(Bca_19_T=Bca_19_N,
         Bca_20_T=Bca_20_N,
         Bca_39_T=Bca_39_N)

TPMs2<-bind_cols(case,ctrl)
TPMs3<-TPMs2%>%
  select(1:12,88,89,15:31,104,33:87,13,14,90:103,32,105:127)
write_tsv(TPMs3,"Protein_coding.tpms.txt")
####################################################################
########fix epms#######################################################
EPMs<-readRDS("../../Bladder/gene_drops.RDS")
case<-EPMs%>%
  select(1:81)
ctrl<-EPMs%>%
  select(82:161)
case<-case%>%
  rename(cBca_19N=cBca_19T,
         cBca_20N=cBca_20T,
         cBca_39N=cBca_39T)
ctrl<-ctrl%>%
  rename(cBca_19T=cBca_19N,
         cBca_20T=cBca_20N,
         cBca_39T=cBca_39N)
EPMs2<-bind_cols(case,ctrl)
EPMs3<-EPMs2%>%
  select(1:19,100,101,22:38,119,40:99,20,21,102:118,39,120:161)
saveRDS(EPMs3,"../../Bladder/gene_drops.RDS")
##############################################################
### fix counts sample###################################################
case<-counts%>%
  dplyr::select(1:71)
ctrl<-counts%>%
  dplyr::select(72:127)
case<-case%>%
  dplyr::rename(Bca_19_N=Bca_19_T,
         Bca_20_N=Bca_20_T,
         Bca_39_N=Bca_39_T)
ctrl<-ctrl%>%
  dplyr::rename(Bca_19_T=Bca_19_N,
         Bca_20_T=Bca_20_N,
         Bca_39_T=Bca_39_N)

counts2<-bind_cols(case,ctrl)
counts3<-counts2%>%
  dplyr::select(1:12,88,89,15:31,104,33:87,13,14,90:103,32,105:127)
write_tsv(counts3,"allSample.count.txt")
# CNV ---------------------------------------------------------------------
preprocess<-function(x){
  df<-read_tsv(stringr::str_glue("CNV/annotate/{x}T.annotate.cns"))
  df%>%
    select(gene,log2)%>%
    tidyr::separate_rows(gene,sep = ",")%>%
    filter(gene%in%(allgenes$X4))%>%
    arrange(gene,desc(log2))%>%
    setNames(c("Geneid",paste0("Bca_",x,"T")))%>%
    group_by(Geneid)%>%
    slice_head(n=1)%>%
    distinct(Geneid,.keep_all = T)
}

samples<-c("1","2","3","4","5","6","7",
           "8","9","10","11","12","13",
           "14","15","16","17","18","19",
           "20","21","23","24","25","26",
           "27","28","29","30","31",
           "33","34","35","36","37","38",
           "39","40","41","42","43","44",
           "45","46","47","48","50",
           "51","52","54","55","56","57",
           "58","59","60","62","63","64",
           "65","66","67","68","69","70",
           "71","72","74","75","76","77",
           "78","79","80","84","85","86",
           "87","88")

CNVs1<-Reduce(function(x,y){full_join(x,y,by="Geneid")},lapply(samples,preprocess))
write_tsv(CNVs1,"CNV/all_CNV_bigger.tsv")

# multiOmics --------------------------------------------------------------
gongtong<-intersect(intersect(stringr::str_remove(stringr::str_remove(names(TPMs)[1:71],"Bca_"),"_T"),
                              stringr::str_remove(stringr::str_remove(names(EPMs)[1:81],"cBca_"),"T")
),stringr::str_remove(stringr::str_remove(names(CNVs),"Bca_"),"T"))

TPMs<-read_tsv("Protein_coding.tpms.txt")
mat_RNA<-TPMs%>%
  select(1:71)%>%
  setNames(stringr::str_remove(stringr::str_remove(names(TPMs)[1:71],"Bca_"),"_T"))%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()
mat_RNA<-mat_RNA[,gongtong]

EPMs<-readRDS("../../Bladder/gene_drops.RDS")
mat_ecc<-EPMs%>%
  select(1:81)%>%
  setNames(stringr::str_remove(stringr::str_remove(names(EPMs)[1:81],"cBca_"),"T"))%>%
  tibble::column_to_rownames(var="gene")%>%
  as.data.frame()
mat_ecc<-mat_ecc[,gongtong]

CNVs<-read_tsv("../CNV/annotate/all_CNV.tsv")
mat_cnv<-CNVs%>%
  setNames(stringr::str_remove(stringr::str_remove(names(CNVs),"Bca_"),"T"))%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()
mat_cnv<-mat_cnv[,gongtong]

sig <- calculateCorForTwoMatrices(matrix1=mat_ecc,matrix2=mat_RNA,fdr=0.01)
sig<-readRDS("sig2.RDS")

genelocate<-read.table("/home/panxiaoguang/R/x86_64-redhat-linux-gnu-library/4.1/multiOmicsViz/extdata/genelocate.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
chromLength <- read.table("/home/panxiaoguang/R/x86_64-redhat-linux-gnu-library/4.1/multiOmicsViz/extdata/chromLength.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
chrome_sourceOmics <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
chrome_targetOmics <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
genelocate_sourceOmics <- genelocate[genelocate[,2] %in%chrome_sourceOmics,]
genelocate_targetOmics <- genelocate[genelocate[,2] %in%chrome_targetOmics,]
re <-calculateChromLength(chromLength,chrome_sourceOmics,genelocate_sourceOmics)
genelocate_sourceOmics <- re$genelocate
chromLength_sourceOmics <- re$chromLength
re <-calculateChromLength(chromLength,chrome_targetOmics,genelocate_targetOmics)
genelocate_targetOmics <- re$genelocate
chromLength_targetOmics <- re$chromLength


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
    color <- ifelse(cov>0,"#9E3735","#48436D")
    
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
      
      abline(h=c(0,chromLength_targetOmics[,4]),v=c(0,chromLength_sourceOmics[,4]),
             col="gray",lty=3)
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
               numP=apply(corrArray,1,function(x){sum(x==1)}),
               numN=0-apply(corrArray,1,function(x){sum(x==-1)}))
  raw_gene_coord<-genelocate_sourceOmics%>%
    as_tibble()%>%
    select(1,5)%>%
    rename(gene=Symbol)
  haha<-haha%>%
    left_join(raw_gene_coord,by="gene")
  par(mar=c(4,4,0,0))
  plot(0,0,xlim=c(0,allChromlen),ylim=c(-1583,503),type="n",xaxt="n",
       frame.plot=FALSE,xlab=paste(sourceOmicsName," chromosomal
      location",sep=""),
       ylab="Number of significant \ncorrelations")
  axis(side=1,at=(chromLength_sourceOmics[,4]-
                    chromLength_sourceOmics[,2]/2),
       labels=chromLength_sourceOmics[,1])
  axis(side=2,at=NULL,labels=TRUE)
  abline(v=c(0,chromLength_sourceOmics[,4]),col="gray",lty=3)
  points(haha$finalstart,haha$numP,cex=0.2,type="h",col="#B2B1B2")
  points(haha$finalstart,haha$numN,cex=0.2,type="h",col="#060605")
}

#png("eccDNA_mRNA_heatmap.png",width = 4.71,height =7.26,units = "in",res=300)
pdf("eccDNA_mRNA_heatmap.pdf",width = 4.71,height =7.26)
layout(matrix(c(1,2),2),heights=c(2,1))
plotHeatMap(sig,genelocate_sourceOmics,chromLength_sourceOmics,
            genelocate_targetOmics,chromLength_targetOmics,"eccDNA",
            "mRNA",dim=1)
plotSummaryBar(sig,chromLength_sourceOmics,
               genelocate_sourceOmics,"eccDNA")
dev.off()
### for two matrix see "multiOmics.R"
############################################################################################
TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
EPMs<-readRDS("../Bladder/gene_drops.RDS")
CNVs<-read_tsv("CNV/annotate/all_CNV.tsv")

TPMs<-TPMs%>%
  select(1:71)%>%
  setNames(stringr::str_remove(stringr::str_remove(names(TPMs)[1:71],"Bca_"),"_T"))

EPMs<-EPMs%>%
  select(1:81)%>%
  setNames(stringr::str_remove(stringr::str_remove(names(EPMs)[1:81],"cBca_"),"T"))
gongtong<-intersect(names(TPMs),names(EPMs))

names(EPMs)[1]<-"Geneid"
TPMs<-TPMs%>%
  select("Geneid",all_of(gongtong))
EPMs<-EPMs%>%
  select("Geneid",all_of(gongtong))
TPMs<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)
EPMs<-EPMs%>%
  tidyr::gather(sample,EPM,-Geneid)
TPMs<-TPMs%>%
  tidyr::replace_na(replace = list(TPM=0))
EPMs<-EPMs%>%
  tidyr::replace_na(replace = list(EPM=0))
fin<-inner_join(TPMs,EPMs,by=c("Geneid","sample"))
fin<-na.omit(fin)
fin<-fin%>%
  mutate(logT=log2(TPM_mean+1),
         logE=log2(EPM_mean+1))
fin2<-fin%>%
  sample_frac(0.002,replace = F)
fin2<-fin2%>%
  mutate(type=if_else((logE-logT >= 5) | logT<=1,"non-enhanced","enhanced"))

ggscatterhist(
  fin, x = "logE", y = "logT", color="clu",size = 0.2, alpha = 0.6,
  palette = c("#9E3735","#48436D"),
  margin.params = list(fill="clu",color = "black", size = 0.2),
  xlab = "junction counts",ylab = "TPMs"
)
fin<-fin%>%
  mutate(type=if_else((logE-logT >= 5) | logT<=1,"non-enhanced","enhanced"))
ggplot(fin,aes(x=logE,y=logT))+
  geom_point(aes(color=type),size=0.2,alpha=0.6)+
  geom_smooth(aes(color=type),method = "lm")+
  xlab("Junction count")+
  ylab("RNA expression")

ggsave("TPM.vs.junction.pdf",width =4.19 ,height =4.39 )
write_tsv(fin2,"TPM.vs.junction.tsv")

####################if for all epms ?
epms2<-openxlsx::read.xlsx("../Bladder/Bladder_circle_numbers.xlsx")
pms<-epms2%>%
  as_tibble()%>%
  filter(group=="Tumour")%>%
  mutate(sample2=stringr::str_remove(stringr::str_remove(sample,"cBca_"),"T"))%>%
  select(sample2,EPM)%>%
  dplyr::rename(sample=sample2)
TPMs<-TPMs%>%
  select(1:71)%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  mutate(sample=stringr::str_remove(stringr::str_remove(sample,"Bca_"),"_T"))

wocai<-TPMs%>%
  left_join(pms,by="sample")

wocai<-wocai%>%
  tidyr::replace_na(replace = list(EPM=0,TPM=0))
rst<-wocai%>%
  group_by(Geneid)%>%
  rstatix::cor_test(vars = "TPM",vars2 = "EPM",method = "spearman")
rst<-rst%>%filter(!is.na(cor))
rst<-rst%>%mutate(tp=if_else(p<0.05,"yes","no"))
ggplot(rst,aes(x=cor))+geom_histogram(aes(fill=tp),color="black",bins =98)+scale_fill_manual(values = c("yes"="#08519c","no"="#feb24c"))+theme_pubr()+ylab("Gene count")+geom_vline(xintercept = 0.083,linetype="dashed",color="red")+xlab("Correlation")
ggsave("cor.pdf",width = 6.64,height =3.08 )
part1<-rst%>%filter(cor>0,p<0.05)
part2<-rst%>%filter(cor<0,p<0.05)
#################################################
###for normal
nms<-stringr::str_remove(stringr::str_remove(names(TPMs)[c(1,72:127)],"Bca_"),"_N")
TPMs<-TPMs%>%
  dplyr::select(1,72:127)%>%
  setNames(nms)
nms<-stringr::str_remove(stringr::str_remove(names(EPMs)[c(1,82:161)],"cBca_"),"N")
EPMs<-EPMs%>%
  dplyr::select(1,82:161)%>%
  setNames(nms)
gongtong<-intersect(names(TPMs),names(EPMs))

names(EPMs)[1]<-"Geneid"
TPMs<-TPMs%>%
  dplyr::select("Geneid",all_of(gongtong))
EPMs<-EPMs%>%
  dplyr::select("Geneid",all_of(gongtong))
TPMs<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)
EPMs<-EPMs%>%
  tidyr::gather(sample,EPM,-Geneid)
TPMs<-TPMs%>%
  tidyr::replace_na(replace = list(TPM=0))
EPMs<-EPMs%>%
  tidyr::replace_na(replace = list(EPM=0))
fin<-inner_join(TPMs,EPMs,by=c("Geneid","sample"))
rst<-fin%>%mutate(logT=log2(TPM+1),logE=log2(EPM+1))%>%group_by(Geneid)%>%rstatix::cor_test(vars = "logE",vars2 = "logT",method = "spearman")
rst<-rst%>%filter(!is.na(cor))
rst<-rst%>%mutate(tp=if_else(p<0.05,"yes","no"))
ggplot(rst,aes(x=cor))+geom_histogram(aes(fill=tp),color="black",bins = 120)+scale_fill_manual(values = c("yes"="#08519c","no"="#feb24c"))+theme_pubr()+ylab("Gene count")+geom_vline(xintercept = 0.068,linetype="dashed",color="red")+xlab("Correlation")

part1<-rst%>%filter(cor>0,p<0.05)
part2<-rst%>%filter(cor<0,p<0.05)
geneid<-clusterProfiler::bitr(part1$Geneid,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
ek<-clusterProfiler::enrichKEGG(gene = geneid$ENTREZID,organism = "hsa",pvalueCutoff = 0.05)
ek<-clusterProfiler::setReadable(ek,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
eg<-clusterProfiler::enrichGO(gene = geneid$ENTREZID,OrgDb = org.Hs.eg.db,ont = "BP",pvalueCutoff = 0.05)
kg<-(ek@result)%>%as_tibble()
go<-(eg@result)%>%as_tibble()
openxlsx::write.xlsx(list(kegg=kg,go=go),"Cancer_gene_cor.enrich.xlsx")
saveRDS(rst,"normal_gene_cor.rds")
####################finish################################

#### ecDNA exp .vs eccDNA exp
### use 20% coding as ecc-gene(Fig4E)
TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
## should filter
goodGenes<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  mutate(type=if_else(TPM>5,"good","bad"))%>%
  group_by(Geneid,type)%>%
  summarise(count=n())%>%
  tidyr::pivot_wider(names_from = "type",values_from = count,values_fill = 0)%>%
  filter(good==126)%>%
  pull(Geneid)
TPMs<-TPMs%>%
  filter(Geneid%in%goodGenes)
TPM2<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)
TPM3<-TPM2%>%filter(stringr::str_ends(sample,"T"))


## mean+sd
stst<-TPM3%>%
  group_by(Geneid)%>%
  summarise(mean=mean(TPM),sd=sd(TPM))
## zscore
zscores<-TPM3%>%
  left_join(stst,by="Geneid")%>%
  mutate(zscore=(TPM-mean)/sd)

zscores<-zscores%>%
  dplyr::select(Geneid,sample,zscore)%>%
  dplyr::rename(TPM=zscore)%>%
  mutate(sample=stringr::str_extract(sample,"\\d+"))

eccgenes<-read_tsv("../Bladder/allgene20.bed",col_names = F)
eccgenes<-eccgenes%>%
  mutate(sample=stringr::str_extract(X4,"cBca_\\d+[NT]"))%>%
  dplyr::select(sample,X8)%>%
  distinct(sample,X8,.keep_all = T)
names(eccgenes)[2]<-"Geneid"
eccgenes<-eccgenes%>%
  filter(stringr::str_ends(sample,"T"))%>%
  mutate(sample=stringr::str_extract(sample,"\\d+"))%>%
  left_join(zscores,by=c("sample","Geneid"))
#eccgenes<-na.omit(eccgenes)


ecgenes<-openxlsx::read.xlsx("all_detective_gene_list.xlsx")%>%
  as_tibble()%>%
  filter(stringr::str_starts(feature,"ecDNA"))%>%
  dplyr::select(sample_name,gene)%>%
  distinct(sample_name,gene,.keep_all = T)%>%
  dplyr::rename(sample=sample_name,Geneid=gene)
ecgenes$sample<-as.character(ecgenes$sample)
ecgenes<-ecgenes%>%
  left_join(zscores,by=c("sample","Geneid"))%>%
  na.omit()
ecgenes<-ecgenes%>%
  mutate(ding=paste(sample,Geneid,sep="-"))
eccgenes<-eccgenes%>%
  mutate(ding=paste(sample,Geneid,sep="-"))
### use all focal-amp genes to filter
allGenes<-openxlsx::read.xlsx("all_detective_gene_list.xlsx")%>%
  as_tibble()%>%
  dplyr::select(sample_name,gene)%>%
  distinct(sample_name,gene,.keep_all = T)%>%
  dplyr::rename(sample=sample_name,Geneid=gene)%>%
  mutate(ding=paste(sample,Geneid,sep="-"))
jiaoji<-intersect(allGenes$ding,eccgenes$ding)
eccgenes<-eccgenes%>%
  filter(!ding%in%jiaoji)
eccgenes$type<-"eccDNA"
ecgenes$type<-"ecDNA"

fin<-bind_rows(ecgenes,eccgenes)

fin<-fin%>%
  na.omit()
fin<-fin%>%
  mutate(label=if_else(ding %in% c("87-NFS1","87-CPNE1","86-CCDC91",
                                   "32-IFRD2","59-CD34","23-FOXO4"),Geneid,""))

ggplot(fin,aes(x=type,y=TPM))+
  geom_jitter(aes(color=type),
              size=1,
              alpha=0.4)+
  geom_violin(aes(fill=type))+
  geom_text_repel(aes(label=label),max.overlaps = 100000)+
  scale_fill_manual(values = c("#9E3735","#48436D"))+
  scale_color_manual(values = c("#9E3735","#48436D"))+
  geom_hline(yintercept = 0,linetype="dashed")+
  xlab("")+
  ylab("RNA expression(Z-score)")+
  theme_pubr()

ggplot(fin,aes(x=TPM,color=type))+
  geom_density()+
  scale_color_manual(values = c("#9E3735","#48436D"))+
  theme_pubr()+
  xlab("RNA expression(Z-score)")

ggsave("ecDNA_vs_eccDNA2.pdf",width = 3.55,height = 4.03)
ggsave("ecDNA_vs_eccDNA3.pdf",width = 5.52,height = 2.97)
write_tsv(fin,"ecDNA_vs_eccDNA.tsv")
####sample-65

wocao1<-TPM3%>%filter(sample=="Bca_65_T")
wocao1<-wocao1%>%mutate(rank=dense_rank(dplyr::desc(TPM)))
wocao1<-wocao1%>%mutate(lab=if_else(Geneid%in%c("FTL","RPS6","PCBP1","MYL6","AP3S2"),Geneid,""),
                lab2=if_else(Geneid%in%c("FTL","RPS6","PCBP1","MYL6","AP3S2"),"yes","no"))

p1<-ggplot(wocao1,aes(x=rank,y=log2(`TPM`+1)))+
  geom_line()+
  geom_point(aes(color=lab2))+
  ggrepel::geom_text_repel(aes(label=lab),
            max.overlaps = 100000000000)+
  scale_color_manual(values = list(yes="red",no="transparent"))+
  theme_pubr()+
  xlab("Genes ranked by expression level")+
  ylab("Expression level")+
  theme(legend.position = "none")

wocao2<-CNVs%>%select(Geneid,Bca_65T)%>%mutate(sample=65)%>%mutate(rank=dense_rank(desc(Bca_65T)))
wocao2<-wocao2%>%mutate(lab=if_else(Geneid%in%c("FTL","RPS6","PCBP1","MYL6","AP3S2"),Geneid,""),
                        lab2=if_else(Geneid%in%c("FTL","RPS6","PCBP1","MYL6","AP3S2"),"yes","no"))

p2<-ggplot(wocao2,aes(x=rank,y=`Bca_65T`))+
  geom_line()+
  geom_point(aes(color=lab2))+
  ggrepel::geom_text_repel(aes(label=lab),
                  max.overlaps = 100000000000)+
  scale_color_manual(values = list(yes="red",no="transparent"))+
  theme_pubr()+
  xlab("Genes ranked by CNA level")+
  ylab("log2 Ratio")+
  theme(
        legend.position = "none")
wocao3<-EPMs%>%select(gene,cBca_65T)%>%mutate(sample=65)
wocao3$haha<-log2(wocao3$cBca_65T)
wocao3<-wocao3%>%mutate(rank=dense_rank(desc(haha)))
wocao3<-wocao3%>%mutate(lab=if_else(gene%in%c("FTL","RPS6","PCBP1","MYL6","AP3S2"),gene,""),
                        lab2=if_else(gene%in%c("FTL","RPS6","PCBP1","MYL6","AP3S2"),"yes","no"))
p3<-ggplot(wocao3,aes(x=rank,y=haha))+
  geom_line()+
  geom_point(aes(color=lab2))+
  ggrepel::geom_text_repel(aes(label=lab),
                  max.overlaps = 100000000000)+
  scale_color_manual(values = list(yes="red",no="transparent"))+
  theme_pubr()+
  xlab("Genes ranked by circular density")+
  ylab("Junction count")+
  theme(
        legend.position = "none")
p1+p2+p3
ggsave("RNAseq/sample_65.compare.pdf",width = 8.55,height = 3.95)
####################sample16################################################################
ecgenes<-ecgenes%>%filter(sample=="16")

plotData<-TPM3%>%
  mutate(sample=stringr::str_extract(sample,"\\d+"))%>%
  filter(sample=="16")%>%
  mutate(rank=dense_rank(desc(TPM)))%>%
  mutate(ecgene=if_else(Geneid%in%c("PABPC1","BRK1","YWHAZ","TADA3","CTNNB1","OXR1","UBR5","ANKRD46"),Geneid,""))%>%
  mutate(col=if_else(ecgene=="","no","yes"))


ggplot(plotData,aes(x=rank,y=TPM))+
  geom_line()+
  geom_point(aes(color=col))+
  scale_y_log10()+
  scale_color_manual(values = list(yes="red",no="transparent"))+
  geom_text_repel(aes(label=ecgene),max.overlaps = 100000000000)+
  theme(legend.position="none")+
  xlab("Genes ranked by expression level\n(TPM>5)")+
  ylab("TPM")+
  theme_pubr()

ggsave("BAF/sample16.gene.rank.pdf",width = 3.89,height = 5.42)
#################################################################################
###only true eccDNA

data<-read_tsv("../../Bladder/bedFile/only_true_eccDNA.bed",col_names = F)
new<-data%>%
  mutate(length=X3-X2)
ecDNA<-read_tsv("../../Bladder/bedFile/all_ecDNA.bed",col_names = F)
eccDNA<-new%>%
  sample_frac(0.002)
eccDNA$type<-"eccDNA"
eccDNA<-eccDNA%>%select(1:3,5:6)
ecDNA<-ecDNA%>%
  mutate(length=X3-X2)
ecDNA$type<-"ecDNA"
plotData<-bind_rows(eccDNA,ecDNA)
plotData$logd<-log10(plotData$length)
ggviolin(plotData,x="type",y="length",
         fill="type",yscale = "log10",
         palette = c("#D1342B","#393A80"),trim = T,
         xlab="",
         ylab="Circle length (bp)")
ggsave("length_compare.ecDNA_eccDNA.pdf",width = 3.94,height = 3.51)
###different segment with sample-class
sm<-read_tsv("../class_sample.tsv")
sm<-sm%>%mutate(type=if_else(Sample_classification=="ecDNA","ecDNA","Non-ecDNA"))%>%
  mutate(sample=stringr::str_remove(Sample,"Bca_"))

segments<-read_tsv("segnumbers.tsv")
segments<-segments%>%
  mutate(sample=stringr::str_remove(sample,"N.cs.rmdup.sort.call.cns"))
sm<-sm%>%
  select(sample,type)%>%
  left_join(segments,by="sample")
sm<-sm%>%
  arrange(type)
write_tsv(sm,"type-segnumbers.tsv")

## try to plot ecc-gene with other gene
corrArray <-readRDS("final.ecc2na.RDS")
test<-tibble(gene=rownames(corrArray),num=apply(corrArray,1,function(x){sum(x!=0)}))
test<-test%>%
  arrange(desc(num))

xiangguanxing<-tibble(gene=colnames(corrArray),lnks=as.numeric(corrArray["RTKN2",]))
xiangguanxing<-xiangguanxing%>%filter(lnks!=0)
genelocate<-read_tsv("../../Bladder/dbs/hg38.coding.bed",col_names = F)
genelocate<-genelocate%>%
  select(4,1,2,3)%>%
  mutate(X1=stringr::str_remove(X1,"chr"))%>%
  setNames(c("Symbol","chrom","start","end"))

source_bed<-genelocate%>%
  filter(Symbol=="RTKN2")%>%
  select(2:4)%>%
  mutate(chrom=paste0("chr",chrom))%>%
  as.data.frame()

target_bed<-genelocate%>%
  as_tibble()%>%
  filter(Symbol%in%(xiangguanxing$gene))%>%
  select(2:4)%>%
  mutate(chrom=paste0("chr",chrom))%>%
  as.data.frame()

source_bed_2<-data.frame(chrom=rep(source_bed$chrom,nrow(target_bed)),
                         start=rep(source_bed$start,nrow(target_bed)),
                         end=rep(source_bed$end,nrow(target_bed)))

pdf("RTKN2.ecc_cor2.pdf",width = 6.53,height =4.84 )
circos.clear()
circos.initializeWithIdeogram(plotType = c("axis", "labels"),species = "hg38")
circos.genomicLink(source_bed_2, target_bed, col = if_else(xiangguanxing$lnks==1,"#E63126","#0932E3"),border = NA)
dev.off()
## gsva to mRNA
### filter gene > 0  in more than 50% samples
sm<-read_tsv("../class_sample.tsv")
sm<-sm%>%mutate(type=if_else(Sample_classification=="ecDNA","ecDNA","Non-ecDNA"),
                Sample=paste0(Sample,"_T"))
names(sm)[1]<-"sample"
TPMs<-read_tsv("Protein_coding.tpms.txt")
mat<-TPMs%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.matrix()
mat<-log2(mat+1)
flist<-kOverA(ncol(mat)*0.5,0)
newM<-genefilter(mat,flist)
haha<-mat[newM,]
geneSet<-getGmt("c2.cp.kegg.v7.5.1.symbols.gmt")
fin<-GSVA::gsva(haha,geneSet,min.sz=10,max.sz=500,method="ssgsea")

finallyData<-fin%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="geneSet")%>%
  as_tibble()%>%
  filter(geneSet=="KEGG_BLADDER_CANCER")%>%
  tidyr::gather(sample,ssGSEA,-geneSet)%>%
  filter(stringr::str_ends(sample,"_T"))%>%
  left_join(sm,by="sample")
finallyData%>%
  rstatix::t_test(ssGSEA~type,alternative = "greater")
ggboxplot(finallyData,
          x="type",
          y="ssGSEA",
          color="type",
          add = "jitter",
          palette = c("#393A80","#D1342B"))+
  stat_compare_means(method.args = list(alternative = "greater"),method = "t.test")

ggsave("bladder_cancer_diff_pathway.pdf",width = 2.84,height = 4.2)

dffs<-fin%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="geneSet")%>%
  as_tibble()%>%
  tidyr::gather(sample,ssGSEA,-geneSet)%>%
  filter(stringr::str_ends(sample,"_T"))%>%
  left_join(sm,by="sample")%>%
  group_by(geneSet)%>%
  rstatix::t_test(ssGSEA~type,alternative = "greater")
## use CNA detected arm-level case
broad<-read_tsv("../CNV/gistic2_results/broad_values_by_arm.txt")
sampleNumbers<-ncol(broad)-1
det<-broad%>%
  as_tibble()%>%
  tidyr::gather(sample,value,-`Chromosome Arm`)%>%
  mutate(cases=case_when(value>0 ~ "amp",
                         value==0 ~ "balance",
                         value<0 ~ "dep"))%>%
  group_by(`Chromosome Arm`,cases)%>%
  summarise(num=n(),
            pct=num/77*100)%>%
  arrange(desc(pct))
### tumor suppressor
database<-openxlsx::read.xlsx("../Cancer_driver genes_Bladder cancer.xlsx")
tsg_genes<-database%>%
  filter(stringr::str_detect(Role,"TSG"))%>%
  pull(Gene)%>%
  unique()
delGenes<-openxlsx::read.xlsx("../CNV/del_genes.xlsx")
detecd_del_genes<-intersect(tsg_genes,unique(delGenes$genes))
## get del genes freq
focal<-read_tsv("../CNV/gistic2_results/focal_data_by_genes.txt")
focao_dep_genes<-focal%>%
  select(-2,-3)%>%
  tidyr::gather(sample,value,-`Gene Symbol`)%>%
  mutate(cases=case_when(value>0 ~ "amp",
                         value==0 ~ "balance",
                         value<0 ~ "dep"))%>%
  group_by(`Gene Symbol`,cases)%>%
  summarise(num=n(),
            pct=num/77*100)%>%
  arrange(desc(pct))%>%
  filter(cases=="dep")%>%
  filter(`Gene Symbol`%in%detecd_del_genes)
# compare amp with no-amp ecc-count

chromSize<-read_tsv("../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size",col_names = F)
total_chrom_size<-chromSize%>%
  filter(X1%in%c(paste0("chr",seq(1,22)),"chrX","chrY"))%>%
  summarise(total=sum(X2))%>%
  pull(total)


get_stat<-function(x,total_chrom_size=3088269832){
  seed_region<-read_tsv(stringr::str_glue("../Bladder/seeds/{x}_AA_CNV_SEEDS.bed"),col_names = F)
  seed_length<-seed_region%>%
    mutate(length=X3-X2)%>%
    summarise(total=sum(length))%>%
    pull(total)
  inters<-read_tsv(stringr::str_glue("../Bladder/seeds/{x}.inters.bed"),col_names = F)
  eccs<-read_tsv(stringr::str_glue("../Bladder/bedFile/cBca_{x}T.ecc.bed"),col_names = F)
  if(nrow(inters)!=nrow(eccs)){
    cat(x,"may be have some problems!")
  }
  inters<-inters%>%
    filter(X1%in%c(paste0("chr",seq(1,22)),"chrX","chrY"))%>%
    group_by(X5)%>%
    summarise(count=n())%>%
    filter(X5<2)%>%
    arrange(X5)%>%
    mutate(length=c(total_chrom_size-seed_length,seed_length))%>%
    mutate(count2=count/length,
           count3=count2/sum(count2))%>%
    select(X5,count3)%>%
    setNames(c("type","counts"))
  inters$sample<-x
  inters
}

samples<-c("1","2","5","8","9","13","14","15","16","17","21","23","24","25","28","29","30","31","32","33","34","36","37","39","40","41","42","44","46","47","48","50","51","54","56","57","60","62","63","65","67","68","69","70","71","72","74","75","76","77","79","80","84","85","86","87","88")

fin<-do.call("bind_rows",lapply(samples, function(x) get_stat(x)))

plotData<-fin%>%
  mutate(type=if_else(type==0,"non-amp","amp"))

plotData<-plotData%>%
  tidyr::pivot_wider(names_from="sample",values_from = "counts")
write_tsv(plotData,"../ecccount_diff.tsv")
ggboxplot(plotData,x="type",
          y="counts",
          add="jitter",
          color="type",
          palette = c("#393A80","#D1342B"),
          xlab = "",
          ylab = "Reletive ")+
  stat_compare_means(method.args = list(paired = T),method = "t.test")

### CNV cor with exp
AA_derived_CNV<-openxlsx::read.xlsx("all_detective_gene_list.xlsx")
AA_derived_CNV$sample_name<-as.character(AA_derived_CNV$sample_name)

### try cnvkit cna
CNVs<-read_tsv("CNV/annotate/all_CNV.tsv")
CNVs<-CNVs%>%tidyr::gather(sample,cn,-Geneid)%>%mutate(cn=2**(cn+1))
CNVs$sample<-stringr::str_remove(stringr::str_remove(CNVs$sample,"Bca_"),"T")
names(CNVs)[1:2]<-c("gene","sample_name")
## catgory for every genes
AA_derived_CNV<-AA_derived_CNV%>%
  left_join(CNVs,by=c("gene","sample_name"))
AA_derived_CNV<-na.omit(AA_derived_CNV)

allRNA_samples<-stringr::str_remove(stringr::str_remove(names(TPMs)[2:71],"Bca_"),"_T")

TPM_td<-TPM3%>%
  dplyr::rename(sample_name=sample,gene=Geneid)
TPM_td$sample_name<-stringr::str_remove(stringr::str_remove(TPM_td$sample_name,"Bca_"),"_T")

higher_AA<-AA_derived_CNV%>%
  filter(sample_name%in%allRNA_samples)

gene_in_sample<-higher_AA%>%
  group_by(gene)%>%
  summarise(nohave=paste0(setdiff(allRNA_samples,sample_name),collapse = ","))

higher_AA<-higher_AA%>%
  left_join(TPM_td,by=c("sample_name","gene"))
higher_AA<-na.omit(higher_AA)

ecc_amp<-higher_AA%>%
  filter(stringr::str_detect(feature,"ecDNA"))
non_ecc_amp<-higher_AA%>%
  filter(!stringr::str_detect(feature,"ecDNA"))

get_fold<-function(x){
  noneed<-higher_AA%>%
    filter(gene==x)%>%
    pull(sample_name)
  tmp<-TPM_td%>%
    filter(gene==x)%>%
    filter(!(sample_name%in%noneed))%>%
    summarise(noamp=mean(TPM))
  tmp$noamp
}

ecc_amp<-ecc_amp%>%
  mutate(noampTPM=purrr::map_dbl(gene,function(x){get_fold(x)}))
non_ecc_amp<-non_ecc_amp%>%
  mutate(noampTPM=purrr::map_dbl(gene,function(x){get_fold(x)}))

ecc_amp<-ecc_amp%>%
  mutate(FC=(TPM+1)/(noampTPM+1))
non_ecc_amp<-non_ecc_amp%>%
  mutate(FC=(TPM+1)/(noampTPM+1))

ecc_amp$type<-"ecDNA"
non_ecc_amp$type<-"others"
plotData<-bind_rows(ecc_amp,non_ecc_amp)

oncogenes<-openxlsx::read.xlsx("Cancer_driver genes_Bladder cancer.xlsx")
hao<-oncogenes%>%
  #filter(Role%in%c("oncogene","oncogene, fusion"))%>%
  pull(Gene)%>%
  unique()
plotData2<-plotData%>%filter(gene%in%hao)

ggplot(plotData2,aes(x=gene_cn,y=FC,color=type))+
  geom_point()+
  scale_color_manual(values = c("#9E3735","#48436D"))+
  geom_smooth(method = "lm")+
  xlab("AA derived copy count")+
  ylab("Fold change")+
  theme_pubr()
ggplot(plotData,aes(x=gene_cn,y=FC,color=type))+
  geom_point()+
  scale_color_manual(values = c("#9E3735","#48436D"))+
  #scale_x_continuous(limits = c(0,20))+
  #scale_y_continuous(limits = c(0,10))+
  geom_smooth(method = "lm")+
  xlab("AA derived copy count")+
  ylab("Fold change")+
  theme_pubr()
ggsave("ecDNA_vs.none2.pdf",width = 4.34,height = 3.68)

###### amp_class
amp_class<-openxlsx::read.xlsx("下机信息（含barcode）/all_sum.xlsx")

newClass<-amp_class%>%
  as_tibble()%>%
  filter(AmpliconType!="No amp/Invalid")%>%
  select(1:3,7)
names(newClass)<-c("sample","amp","type","length")
breakPoints<-read_tsv("graphs/all.detective.breakpoints.tsv")

breakPoints<-breakPoints%>%
  left_join(newClass,by=c("sample","amp"))
wocao<-breakPoints%>%
  mutate(name=paste(sample,amp,sep = "_"),start=bp-1)%>%
  select(chrom,start,bp,type)
wocao<-na.omit(wocao)
#write_tsv(wocao,"all_breakpoints.bed",col_names = F)

plotData<-breakPoints%>%
  na.omit()%>%
  group_by(sample,amp,type,length)%>%
  summarise(count=n())%>%
  mutate(ratio=count/length*1000000)
plotData<-plotData%>%
  mutate(logratio=log2(ratio))

plotData<-newClass%>%
  left_join(plotData,by=c("sample","amp","type","length"))%>%
  tidyr::replace_na(replace = list(count=0,ratio=0,logratio=0))
my_comparisons <- list( c("ecDNA", "BFB"), 
                        c("Complex non-cyclic", "BFB"), 
                        c("ecDNA", "Complex non-cyclic"),
                        c("ecDNA", "Linear amplification"))
p<-ggboxplot(plotData,x="type",y="logratio",add = "jitter",color="type",ylab = "breakpoint number per million amplicon")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test")
ggpar(p,x.text.angle=30)
ggsave("break_point.firstplot.pdf",width = 3.74,height =5.88 )
write_tsv(plotData,"break_point.stat.tsv")
## circos
chrom_interval<-read_tsv("../DataBase/hg38_1kb_region.bed",col_names = F)
intersect<-read_tsv("breakpoint.intersect.bed",col_names = F)
bedlist<-intersect%>%
  group_by(X7,X1,X2,X3)%>%
  summarise(count=n())%>%
  split(.$X7)

pdf("break_point.secondplot.pdf",width =5.59 ,height =5.19)
circos.clear()
circos.par(track.margin=c(0,0),track.height=0.15)
circos.initializeWithIdeogram(species = "hg38",ideogram.height = 0.03)

circos.genomicTrack(bedlist$ecDNA[,c(2:5)],numeric.column = 4,ylim = c(0, 10),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = value, ybottom = 0,border ="#BC4137")
                    },bg.border = NA, bg.col=NA)

circos.genomicTrack(bedlist$BFB[,c(2:5)],numeric.column = 4,ylim = c(0, 10),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = value, ybottom = 0,border = "#436693")
                    },bg.border = NA, bg.col="#f0f0f0")

circos.genomicTrack(bedlist$`Complex non-cyclic`[,c(2:5)],numeric.column = 4,ylim = c(0, 10),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = value, ybottom = 0,border = "#C58D65")
                    },bg.border = NA, bg.col=NA)

circos.genomicTrack(bedlist$`Linear amplification`[,c(2:5)],numeric.column = 4,ylim = c(0, 10),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = value, ybottom = 0,border = "#A1B2C0")
                    },bg.border = NA, bg.col="#f0f0f0")
dev.off()



reaarange<-bedlist$`Complex non-cyclic`%>%group_by(count)%>%summarise(newC=n())
p3<-ggbarplot(reaarange,x="count",y="newC",position = position_identity(),fill = "#C58D65",width = 0.9,ylab = "counts",xlab = "Breakpoints per 1Kb")+
  geom_smooth(method = "glm",
              lwd=0.4,
              col="black",
              method.args = list(family = "poisson"))+
  scale_x_continuous(breaks = c(0,2,4,6,8),limits = c(0,9))

ecdist<-bedlist$ecDNA%>%group_by(count)%>%summarise(newC=n())
p1<-ggbarplot(ecdist,x="count",y="newC",position = position_identity(),fill = "#BC4137",width = 0.9,ylab = "counts",xlab = "Breakpoints per 1Kb")+
  geom_smooth(method = "glm",
              lwd=0.4,
              col="black",
              method.args = list(family = "poisson"))+
  scale_x_continuous(breaks = c(0,2,4,6,8),limits = c(0,9))

bfb<-bedlist$BFB%>%group_by(count)%>%summarise(newC=n())
p2<-ggbarplot(bfb,x="count",y="newC",position = position_identity(),fill = "#436693",width = 0.9,ylab = "counts",xlab = "Breakpoints per 1Kb")+
  geom_smooth(method = "glm",
              lwd=0.4,
              col="black",
              method.args = list(family = "poisson"))+
  scale_x_continuous(breaks = c(0,2,4,6,8),limits = c(0,9))

lin<-bedlist$`Linear amplification`%>%group_by(count)%>%summarise(newC=n())
p4<-ggbarplot(lin,x="count",y="newC",position = position_identity(),fill = "#A1B2C0",width = 0.9,ylab = "counts",xlab = "Breakpoints per 1Kb")+
  geom_smooth(method = "glm",
              lwd=0.4,
              col="black",
              method.args = list(family = "poisson"))+
  scale_x_continuous(breaks = c(0,2,4,6,8),limits = c(0,9))

ggarrange(p1,p2,p3,p4,ncol = 1,align = "v")
ggsave("break_point.thirdplot.pdf",width = 3.52,height = 5.58)

#################################################################
wocao1<-openxlsx::read.xlsx("statistic_about_amp.xlsx",sheet = 3)
test<-bind_rows(tibble(type="Linear",size=wocao1$Linear),
tibble(type="Heavily-rearranged",size=wocao1$`Heavily-rearranged`),
tibble(type="BFB",size=wocao1$BFB),
tibble(type="ecDNA",size=wocao1$ecDNA))

bk.test<-test%>%
  na.omit()%>%
  rstatix::wilcox_test(size ~ type)

openxlsx::write.xlsx(list(size=sz.test,copy=cp.test,breaks=bk.test),"statistic_about_amp2.xlsx")

####### ecc density rank
EPMs<-readRDS("../Bladder/gene_drops.RDS")
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
ggsave("eccDNA_density_rank.pdf",width = 4.12,height = 4.45)
## total=19071
## tumour-high: H4C11,GPR148,SSU72P7,H2BC3,SCYGR6
## tumor-low: USP9Y,CATSPERB,ADAM32,RPH3AL,ACACA
## Normal-High: H2AC20,ETDC,H2AC14,KRTAP19-1,LCE3A
## Normal-low: PCDH11Y,OTUD7A,USP9Y,NLGN4Y,BMERB1
#########################################
ggplot(newEPM,aes(x=logc,color=samples))+
  geom_density(size=0.3)+
  scale_color_manual(values = c(rep("#c5362c",80),rep("#4475a7",80)))+
  theme_pubr()+
  theme(legend.position = "none",axis.title.x = element_blank())+
  ylab("Density")
ggsave("ecc-abundance-density.pdf",width = 4.5,height = 2.72)
######################################################################
##calBaf

#######################WGS##################################

df<-read_tsv("MYEOV/MYEOV.WGS.allele.txt")
df<-df%>%
  filter(Good_depth>10)
wgs_alle<-df%>%
  tidyr::gather(alleleType,count,-`#CHR`,-POS,-Good_depth)%>%
  arrange(`#CHR`,POS,Good_depth,desc(count))%>%
  group_by(`#CHR`,POS,Good_depth)%>%
  slice_head(n=2)%>%
  ungroup()%>%
  mutate(tp=stringr::str_sub(alleleType,start=7,end=7))%>%
  mutate(AF=count/Good_depth)%>%
  mutate(tp2=rep(c("major","minor"),times=242))

data<-wgs_alle%>%
  select(`#CHR`,POS,AF,tp2)%>%
  tidyr::pivot_wider(names_from = "tp2",values_from = "AF")%>%
  arrange(POS)%>%
  as.data.frame()
gtrack <- GenomeAxisTrack(range=IRanges(start = 69253332,
                                        end = 70626244,
))

gr <- GRanges(seqnames = "chr8", strand = "*",
              ranges = IRanges(start = data$POS, width = 1),
              major=data[,3],minor=data[,4])

dtrack<-DataTrack(gr, name = "WGS AF",
                  groups=c("major","minor"),
                  col=c("#c5362c","#4475a7"),
                  type="p",
                  yTicksAt=c(0,0.5,1),
                  ylim=c(0,1),grid=T,lty.grid="dashed",lwd.grid=0.5,h=3)

wgsCOV<-read_tsv("PABPC1/WGS.cov.bdg",col_names = F)

cr <- GRanges(seqnames = "chr8", strand = "*",
              ranges = IRanges(start = wgsCOV$X2, end = wgsCOV$X3),
              count=wgsCOV$X4)
covtrack<-DataTrack(cr, name = "WGS cov",col="grey",type="histogram")

################################################################################
###############################RNAseq##########################################

df2<-read_tsv("PABPC1/PABPC1.RNA.allele.txt")
df2<-df2%>%
  filter(Good_depth>15)%>%
  mutate(Count_A=if_else(Count_A>2,Count_A,0),
         Count_C=if_else(Count_C>2,Count_C,0),
         Count_G=if_else(Count_G>2,Count_G,0),
         Count_T=if_else(Count_T>2,Count_T,0))%>%
  mutate(Good_depth=Count_A+Count_C+Count_G+Count_T)
RNA_alle<-df2%>%
  tidyr::gather(alleleType,count,-`#CHR`,-POS,-Good_depth)%>%
  mutate(tp=stringr::str_sub(alleleType,start=7,end=7),
         AF=count/Good_depth)%>%
  select(`#CHR`,POS,tp,AF)%>%
  dplyr::rename(RNAAF=AF)

data2<-wgs_alle%>%
  inner_join(RNA_alle,by=c("#CHR","POS","tp"))%>%
  select(`#CHR`,POS,RNAAF,tp2)%>%
  tidyr::pivot_wider(names_from = "tp2",values_from = "RNAAF")%>%
  arrange(POS)%>%
  as.data.frame()


gr2 <- GRanges(seqnames = "chr8", strand = "*",
              ranges = IRanges(start = data2$POS, width = 1),
              major=data2[,3],minor=data2[,4])

dtrack2<-DataTrack(gr2, 
                   name = "RNAcount",
                   groups=c("major","minor"),
                   col=c("#c5362c","#4475a7"),
                   yTicksAt=c(0,0.5,1),
                   ylim=c(0,1),grid=T,lty.grid="dashed",lwd.grid=0.5,h=3)

rnaCOV<-read_tsv("PABPC1/RNA.cov.bdg",col_names = F)

cr2 <- GRanges(seqnames = "chr8", strand = "*",
              ranges = IRanges(start = rnaCOV$X2, end = rnaCOV$X3),
              count=rnaCOV$X4)
covtrack2<-DataTrack(cr2, name = "RNA cov",col="grey",type="histogram",ylim=c(0,500))
##############################################################################
#####################################CIRCLE###################################
df3<-read_tsv("PABPC1/PABPC1.CIRCLE.allele.txt")
df3<-df3%>%
  filter(Good_depth>15)%>%
  mutate(Count_A=if_else(Count_A>2,Count_A,0),
         Count_C=if_else(Count_C>2,Count_C,0),
         Count_G=if_else(Count_G>2,Count_G,0),
         Count_T=if_else(Count_T>2,Count_T,0))%>%
  mutate(Good_depth=Count_A+Count_C+Count_G+Count_T)

circle_alle<-df3%>%
  tidyr::gather(alleleType,count,-`#CHR`,-POS,-Good_depth)%>%
  mutate(tp=stringr::str_sub(alleleType,start=7,end=7),
         AF=count/Good_depth)%>%
  select(`#CHR`,POS,tp,AF)%>%
  dplyr::rename(CIRCLEAF=AF)

data3<-wgs_alle%>%
  inner_join(circle_alle,by=c("#CHR","POS","tp"))%>%
  select(`#CHR`,POS,CIRCLEAF,tp2)%>%
  tidyr::pivot_wider(names_from = "tp2",values_from = "CIRCLEAF")%>%
  arrange(POS)%>%
  as.data.frame()

gr3 <- GRanges(seqnames = "chr8", strand = "*",
               ranges = IRanges(start = data3$POS, width = 1),
               major=data3[,3],minor=data3[,4])

dtrack3<-DataTrack(gr3, name = "CIRCLEcount",groups=c("major","minor"),col=c("#c5362c","#4475a7"),yTicksAt=c(0,0.5,1),
                   ylim=c(0,1),grid=T,lty.grid="dashed",lwd.grid=0.5,h=3)

circleCOV<-read_tsv("PABPC1/CIRCLE.cov.bdg",col_names = F)

cr3 <- GRanges(seqnames = "chr8", strand = "*",
               ranges = IRanges(start = circleCOV$X2, end = circleCOV$X3),
               count=circleCOV$X4)
covtrack3<-DataTrack(cr3, name = "CIRCLE cov",col="grey",type="histogram",ylim=c(0,1000))

png("test_region.png",width = 6.24,height =4.86,units = "in",res=300)
plotTracks(list(covtrack,dtrack,covtrack3,dtrack3,covtrack2,dtrack2,gtrack),
           from = 99497918,
           to=119476539,
           sizes=c(0.5,1,0.5,1,0.5,1,0.5),legend = F,
           background.title = "white",col.title="black",col.axis="black"
           )
dev.off()
###find genes
genes<-read_tsv("PABPC1/snpPos.gene.tsv")
names(genes)<-c("#CHR","POS","GENE")
genes<-genes%>%distinct(`#CHR`,POS,.keep_all=T)
test<-wgs_alle%>%
  left_join(genes,by=c("#CHR","POS"))%>%
  left_join(RNA_alle,by=c("#CHR","POS","tp"))
plotData<-data2%>%tidyr::gather(type,value,-`#CHR`,-POS)%>%left_join(genes,by=c("#CHR","POS"))
ggplot(plotData,aes(x=POS,y=value,color=type))+geom_point()+geom_text_repel(aes(label=GENE),size=1.2,max.overlaps = 10000000000)
########################################
##############Deseq2##########################################
library(DESeq2)
counts<-read_tsv("RNAseq/allSample.count.txt")
counts<-counts%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()
colData<-tibble(sample = colnames(counts),group=c(rep("Tumor",70),rep("Normal",56)))%>%
  mutate(sample2=stringr::str_remove(stringr::str_remove(sample,"_[TN]"),"Bca_"))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group+subject) 
dds <- dds[rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$pvalue), ]
diff_gene <-subset(res, padj < 0.01 & abs(log2FoldChange) > 1) 

### 火山图
allresult<-as.data.frame(res)%>%
  tibble::rownames_to_column(var="Geneid")%>%
  as_tibble()

allresult<-allresult%>%
  mutate(logp=-log10(padj))

allresult<-allresult%>%
  mutate(col=case_when((log2FoldChange<(-1))&(logp>2)~"downregular",
                       (log2FoldChange>1)&(logp>2)~"upregular",
                       TRUE ~ "none"))
ggplot(allresult,aes(x=log2FoldChange,y=logp,color=col))+
  geom_point()+
  geom_vline(xintercept = -1,linetype="dashed")+
  geom_vline(xintercept = 1,linetype="dashed")+
  geom_hline(yintercept = 2,linetype="dashed")+
  scale_color_manual(values = c("upregular"="#C5362C",
                                "downregular"="#4475A7",
                                "none"="#D8D8D8"))+
  xlab("Log2FoldChange")+
  ylab("-log P.adjust")+
  theme_pubr()

##up=1640,dowm=2236  
ggsave("vocanoplot2.pdf",width = 3.96,height =4.32)
dds<-readRDS("different_genes.deseq2.RDS")
test<-counts(dds,normalized=T)
test<-log2(test+1)
ht_data<-as.data.frame(test[rownames(diff_gene),])
pheatmap(ht_data,scale = "row")

ht_data2<-t(apply(ht_data,1,function(x){scale(x,center = T,scale = T)}))
colnames(ht_data2)<-colnames(ht_data)
##########################diffGene_with_ecDNA######################################################
counts<-read_tsv("RNAseq/allSample.count.txt")
counts<-counts%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()
counts<-counts[,1:70]
groupInfo<-read_tsv("class_sample.tsv")
groupInfo$Sample<-paste0(groupInfo$Sample,"_T")
colData<-tibble(Sample = colnames(counts))
colData<-colData%>%left_join(groupInfo,by="Sample")
colData<-colData%>%
  mutate(group=if_else(Sample_classification=="ecDNA","B","A"))%>%
  mutate(forcats::fct_relevel(group,c("A","B")))%>%
  select(Sample,group)%>%
  as.data.frame()
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group) 
dds <- dds[rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds)
res <- results(dds)
allresult<-as.data.frame(res)%>%
  tibble::rownames_to_column(var="Geneid")%>%
  as_tibble()

allresult<-allresult%>%
  mutate(logp=-log10(padj))

allresult<-allresult%>%
  mutate(col=case_when((log2FoldChange<(-1))&(logp>2)~"downregular",
                       (log2FoldChange>1)&(logp>2)~"upregular",
                       TRUE ~ "none"))
labelData<-allresult%>%filter(col!="none")%>%filter(logp>5)
ggplot(allresult,aes(x=log2FoldChange,y=logp,color=col))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label=Geneid),data = labelData,max.overlaps = 10000000000)+
  geom_vline(xintercept = -1,linetype="dashed")+
  geom_vline(xintercept = 1,linetype="dashed")+
  geom_hline(yintercept = 2,linetype="dashed")+
  scale_color_manual(values = c("upregular"="#C5362C",
                                "downregular"="#4475A7",
                                "none"="#D8D8D8"))+
  xlab("Log2FoldChange")+
  ylab("-log P.adjust")+
  theme_pubr()
ggsave("ecDNA_vs_noecDNA.diffgenes.pdf",width =5.19 ,height =5.71 )
anno_colors<-c("#D52126", "#88CCEE", "#FEE52C","#117733", "#CC61B0","#99C945", "#2F8AC4","#332288","#E68316","#661101")
names(anno_colors)<-unique(topResult$Description)
ht<-Heatmap(tmpData,
        show_row_names = F,
        show_column_names = F,
        show_row_dend = F,
        cluster_columns = T,
        cluster_rows = F,
        top_annotation = HeatmapAnnotation(group = colData$group, 
                                           col = list(group = c("A"="#C5362C","B"="#4475A7"))),
        column_split = colData$group,
        right_annotation = rowAnnotation(pathway=topResult$Description,col=list(pathway=anno_colors)),
        name = "Zscore"
)
pdf("diff.heatmap2.pdf",width = 8.73,height =4.57)
draw(ht)
dev.off()
saveRDS(dds,"different_genes.deseq2.ec_vs_noec.RDS")
####################################################################################
##################################GSEA#######################################
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

topResult<-kk2@result[1:10,]%>%
  as_tibble()%>%
  dplyr::select(Description,core_enrichment)%>%
  tidyr::separate_rows(core_enrichment)%>%
  dplyr::rename(Geneid=core_enrichment)%>%
  distinct(Geneid,.keep_all = T)
topResult<-topResult%>%
  filter(Geneid%in%rownames(ht_data2))
#####################################################################################
tmpData<-ht_data2[topResult$Geneid,]
ht<-Heatmap(tmpData,
        show_row_names = F,
        show_column_names = F,
        show_row_dend = F,
        cluster_columns = F,
        cluster_rows = F,
        top_annotation = HeatmapAnnotation(group = c(rep("Tumor",70),rep("Normal",56)), 
                                           col = list(group = c("Tumor"="#C5362C","Normal"="#4475A7"))),
        column_split = c(rep("Tumor",70),rep("Normal",56)),
        right_annotation = rowAnnotation(pathway=topResult$Description),
        name = "Zscore"
)

pdf("diff.heatmap.pdf",width =7.47 ,height = 5.12)
draw(ht)
dev.off()

##find overlap with eccDNA

ecc_diff<-dffs%>%
  filter(p.adj.signif=="****")

length(intersect(rownames(diff_gene),ecc_diff$gene))

## ceshi
allresult<-read_tsv("diffGenes.tsv")
zuobian<-allresult%>%
  filter(col!="none")%>%
  dplyr::select(Geneid,log2FoldChange)
youbian<-nima%>%
  dplyr::select(gene,logfc)
  
names(youbian)[1]<-"Geneid"
wocao<-inner_join(zuobian,youbian,by="Geneid")

wocao<-wocao%>%
  mutate(col=case_when((log2FoldChange>1)&(logfc<(-0.5)) ~ "cls1",
                       (log2FoldChange<(-1))&(logfc<(-0.5)) ~ "cls2",
                       (log2FoldChange>1)&(logfc>0.5) ~ "cls3",
                       (log2FoldChange<(-1))&(logfc>0.5) ~ "cls4"))

ggplot(wocao,aes(x=logfc,y=log2FoldChange))+
  geom_point(aes(color=col),size=1)+
  scale_x_continuous(limits = c(-2,2))+
  scale_y_continuous(limits = c(-6,6))+
  scale_color_manual(values = c("cls1"="#4daf4a","cls2"="#e41a1c",
                                "cls3"="#377eb8","cls4"="#ff7f00"))+
  theme_pubr()

ggsave("haowan.pdf",width = 4.15,height = 4.08)
### corr#################
fin<-readRDS("../Bladder/gene_drops.RDS")
fin<-fin%>%
  dplyr::select(1:81)%>%
  tidyr::gather(sample,value,-gene)%>%
  dplyr::rename(Geneid=gene,EPM=value)%>%
  tidyr::replace_na(replace = list(EPM=0))

TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
goodGenes<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  mutate(type=if_else(TPM>5,"good","bad"))%>%
  group_by(Geneid,type)%>%
  summarise(count=n())%>%
  tidyr::pivot_wider(names_from = "type",values_from = count,values_fill = 0)%>%
  filter(good==126)%>%
  pull(Geneid)

TPMs<-TPMs%>%
  filter(Geneid%in%goodGenes)%>%
  dplyr::select(1:71)%>%
  tidyr::gather(sample,TPM,-Geneid)

fin$sample<-stringr::str_remove(stringr::str_remove(fin$sample,"cBca_"),"T")
TPMs$sample<-stringr::str_remove(stringr::str_remove(TPMs$sample,"Bca_"),"_T")
fin<-inner_join(fin,TPMs,by=c("Geneid","sample"))

rst2<-fin%>%
  group_by(Geneid)%>%
  rstatix::cor_test(vars = "EPM",
                    vars2 = "TPM",
                    method = "spearman")
saveRDS(rst2,"gene_correlation_for_normal.rds")
rst<-readRDS("gene_correlation.rds")
rst2<-readRDS("gene_correlation_for_normal.rds")
rst2<-rst2%>%dplyr::filter(!is.na(cor))
rst2<-rst2%>%mutate(tp=if_else(p<0.05,"yes","no"))
p1<-ggplot(rst,aes(x=cor,y=..density..))+
  geom_histogram(fill="#BD5251",color="black",bins = 100,size=0.2)+
  scale_x_continuous(limits = c(-0.6,0.6))+
  theme(legend.position = "none")+
  theme_pubr()+ylab("Gene Density")+
  geom_vline(xintercept = 0.1,linetype="dashed",color="black")+
  xlab("Correlation")

p1<-ggplot(rst,aes(x=cor))+
  geom_density(fill="#BD5251",color="black",size=0.2,alpha=0.6)+
  scale_x_continuous(limits = c(-0.6,0.6))+
  theme(legend.position = "none")+
  theme_pubr()+ylab("Gene Density")+
  geom_vline(xintercept = 0.1,linetype="dashed",color="black")+
  xlab("Correlation metrics")
p2<-ggplot(rst2,aes(x=cor))+
  geom_density(fill="#4D4FAE",color="black",size=0.2,alpha=0.6)+
  scale_x_continuous(limits = c(-0.6,0.6))+
  theme(legend.position = "none")+
  theme_pubr()+ylab("Gene Density")+
  geom_vline(xintercept = 0.029,linetype="dashed",color="black")+
  xlab("Correlation")
ggarrange(p2,p1,nrow = 2)
ggsave("corplot.pdf",width =6.28 ,height = 3.31)

p2<-ggplot(rst2,aes(x=cor))+
  geom_histogram(fill="#4D4FAE",color="black",bins = 100,size=0.2)+
  scale_x_continuous(limits = c(-0.6,0.6))+
  theme(legend.position = "none")+
  theme_pubr()+ylab("Gene count")+
  geom_vline(xintercept = 0.029,linetype="dashed",color="black")+
  xlab("Correlation")

ggarrange(p2,p1,nrow = 2)
ggsave("corplot.pdf",width = 5.03,height = 4.53)

geneid<-clusterProfiler::bitr(rst$Geneid,
                              fromType = "SYMBOL",
                              toType = "ENTREZID",
                              OrgDb = org.Hs.eg.db)
names(geneid)<-c("Geneid","id")
rst<-rst%>%
  inner_join(geneid,by="Geneid")
###GSEA with corre as rank
## prepare geneList
geneList = rst$cor
names(geneList) = rst$id
geneList = sort(geneList, decreasing = TRUE)
## prepare database
C2_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
C2_t2n <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gs_description)
kk2<-gseKEGG(geneList     = geneList,
             organism     = 'hsa',
             pvalueCutoff = 0.5,
             verbose      = FALSE)
kk2<-setReadable(kk2,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

kk3<-GSEA(geneList,TERM2GENE = C2_t2g,TERM2NAME = C2_t2n,verbose = F)
gseaplot(kk2, geneSetID = 1, title = kk2$Description[1])
########################3
###预后
library(survival)
library(survminer)
sur_info<-openxlsx::read.xlsx("survival.xlsx")
sur_info<-as_tibble(sur_info)
names(sur_info)<-c("sample","status","days")
sur_info$sample<-stringr::str_remove(sur_info$sample,"Bca_")
allresult<-read_tsv("diffGenes.tsv")
diffs<-allresult%>%
  filter(col!="none")
TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
TPM3$sample<-stringr::str_remove(stringr::str_remove(TPM3$sample,"Bca_"),"_T")
TPM3<-TPM3%>%filter(Geneid%in%(diffs$Geneid))
TPM4<-TPM3%>%tidyr::pivot_wider(names_from = "Geneid",values_from = "TPM",values_fill = 0)
TPM5<-TPM4%>%
  mutate(across(2:7927, function(x){if_else(x<=median(x),"low","high")}))%>%
  select(-NPY2R)
sur_info<-sur_info%>%inner_join(TPM5,by="sample")
sur_info<-na.omit(sur_info)

getPvalue<-function(x){
  dff<-survdiff(Surv(days, status) ~ x,data = sur_info2)
  pValue<-1-pchisq(dff$chisq,df=1)
  pValue
}

getHR<-function(x){
  cox<-coxph(Surv(days, status) ~ x, data = sur_info2)
  HR<-(summary(cox))$coefficients[,"exp(coef)"]
  HR
}

rst<-sur_info%>%
  select(4:3878)%>%
  summarise(across(everything(),~getPvalue(.x)))
rst2<-sur_info%>%
  select(4:3878)%>%
  summarise(across(everything(),~getHR(.x)))
haha<-tibble(Geneid=names(rst),pvalue=as.numeric(rst),HR=as.numeric(rst2))
##显著相关
xiangguan<-haha%>%
  filter(pvalue<0.05)
write_tsv(xiangguan,"预后相关基因.tsv")

jiance<-survfit(Surv(days, status) ~ PABPC1,data = sur_info)

ggsurvplot(jiance, 
           data = sur_info,
           conf.int =F,
           pval = TRUE)

ggsave("MYEOV.survival.pdf",width = 4.59,height =4.16)

TPMs%>%
  filter(Geneid=="MYEOV")%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  mutate(gp=if_else(stringr::str_ends(sample,"T"),"Tumor","Normal"))%>%
  mutate(gp=forcats::fct_relevel(gp,c("Tumor","Normal")))%>%
  mutate(smlabel=if_else(sample%in%c("Bca_37_T","Bca_79_T","Bca_86_T","Bca_15_T","Bca_29_T"),sample,""))%>%
  ggplot(aes(x=gp,y=log1p(TPM),color=gp))+
  geom_boxplot()+
  geom_jitter()+
  geom_text_repel(aes(label=smlabel),max.overlaps=10000000,color="black")+
  xlab("")+
  ylab("TPM")
ggsave("MYEOV.expression.pdf",width =3.45 ,height = 3.83)


############### disscuss PABPC1 
TPM2%>%
  filter(sample%in%c("16_T","16_N","17_T","17_N","21_T","21_N","50_T","50_N"))%>%
  filter(Geneid=="PABPC1")%>%
  mutate(type=if_else(stringr::str_ends(sample,"T"),"Tumor","Normal"))%>%
  mutate(logp=log2(TPM+1))%>%
  mutate(sample2=stringr::str_remove(sample,"_[TN]"))%>%
  rstatix::t_test(logp~type,paired = T)
  ggplot(aes(x=type,y=logp,color=type))+
  geom_point()+
  geom_line(aes(group=sample2),color="#bdbdbd",size=0.3)+
  theme_pubr()
ggsave("PABPC1.paired_exp.pdf",width =2.85 ,height = 3.50)

TPM2%>%
  filter(Geneid=="PABPC1")%>%
  mutate(type=if_else(stringr::str_ends(sample,"T"),"Tumor","Normal"))%>%
  mutate(logp=log2(TPM+1))%>%
  mutate(sample2=stringr::str_remove(sample,"_[TN]"))%>%
  ggplot(aes(x=type,y=logp,color=type))+
  geom_point()+
  geom_line(aes(group=sample2),color="#bdbdbd",size=0.3)+
  theme_pubr()

chengdui<-TPM2%>%
  mutate(type=if_else(stringr::str_ends(sample,"T"),"Tumor","Normal"))%>%
  mutate(sample2=stringr::str_remove(sample,"_[TN]"))%>%
  group_by(sample2)%>%
  summarise(count=n())%>%
  filter(count==2)%>%
  pull(sample2)

ri<-TPM2%>%
  mutate(type=if_else(stringr::str_ends(sample,"T"),"Tumor","Normal"))%>%
  mutate(sample2=stringr::str_remove(sample,"_[TN]"))%>%
  select(-sample)%>%
  filter(sample2%in%chengdui)%>%
  tidyr::pivot_wider(names_from = "type",values_from = "TPM")%>%
  mutate(logfc=log2(Tumor/Normal))%>%
  arrange(desc(logfc))%>%
  mutate(fenzu=if_else(logfc>1,"high","low"))%>%
  select(Geneid,sample2,fenzu)%>%
  dplyr::rename(sample=sample2)%>%
  tidyr::pivot_wider(names_from = "Geneid",values_from = fenzu)
ri$sample<-stringr::str_remove(ri$sample,"Bca_")

sur_info2<-sur_info%>%
  filter(sample%in%(ri$sample))%>%
  left_join(ri,by="sample")
sur_info3<-sur_info%>%
  mutate(MYEOV=if_else(sample%in%c("15","29","37","79","86"),"amp","no-amp"))
jiance<-survfit(Surv(days, status) ~ MYEOV,data = sur_info3)

ggsurvplot(jiance, 
           data = sur_info3,
           conf.int =F,
           pval = TRUE)
ggsave("PABPC1.survival2.pdf",width = 5.9,height =3.61)

tmp=rst%>%mutate(rank=dense_rank(desc(cor)))%>%mutate(cosmic=if_else(Geneid%in%(need_genes1$Gene),"cosmic","non-cosmic"))%>%mutate(tp=if_else(cosmic=="cosmic","cosmic",tp))
tmp$tp<-factor(tmp$tp,levels = c("no","yes","cosmic"))

lizi<-tmp%>%
  filter(tp=="cosmic")%>%
  sample_frac(size=0.02)%>%
  arrange(rank)%>%
  pull(Geneid)

tmp<-tmp%>%
  mutate(label=if_else(Geneid%in%c("KDM5A",lizi,"PRIM2"),Geneid,""))

ggplot(tmp,aes(x=rank,y=cor))+
  geom_point(size=0.5)+
  geom_text_repel(aes(label=label),size=3.4,max.overlaps = 10000000000)+
  facet_wrap(~tp)+
  theme_pubr()
ggsave("cor.generank.pdf",width = 12.24,height = 3.73)

### 染色体碎裂
library(ShatterSeek)
library(ggplot2)

huatu<-function(x,y){
  sv = read.table(stringr::str_glue("SV/allmerges/merged.{x}.txt"),header=T)# sv_matrix
  sv$chrom1=sub("chr","",sv$chrom1)
  sv$chrom2=sub("chr","",sv$chrom2)
  a=unique(c(which(!sv$chrom1%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")),
             which(!sv$chrom2%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"))))
  if(length(a)>0){
    sv2 = sv[-a,]
    sv = sv2
  }
  cn = readr::read_tsv(stringr::str_glue("CNV/cns/{x}T.cs.rmdup.sort.call.cns"))
  cn <- as.data.frame(cn)
  cn = cn[,c(1,2,3,6)]
  colnames(cn)=c("chromosome","start","end","total_cn")
  b=unique(c(which(cn$chromosome=="chrY"),which(cn$chromosome=="chrM")))
  if(length(b)>0){
    cn2 = cn[-b,]
    cn = cn2
  }
  cn$chromosome=sub("chr","",cn$chromosome)
  CN_data = CNVsegs(chrom=as.character(cn$chromosome),start=cn$start,end=cn$end,total_cn=cn$total_cn)
  SV_data = SVs(chrom1=as.character(sv$chrom1),pos1=as.numeric(sv$start1),chrom2=as.character(sv$chrom2),pos2=as.numeric(sv$start2),SVtype=as.character(sv$svclass),strand1=as.character(sv$strand1),strand2=as.character(sv$strand2))
  chromothripsis = shatterseek(SV.sample=SV_data,seg.sample=CN_data,genome = "hg38")
  plots=plot_chromothripsis(ShatterSeek_output = chromothripsis,chr=y,genome = "hg38")
  common_ggplot2 <- theme_bw() + theme(axis.text.x=element_text(size=7,angle=0),
                                       axis.text.y=element_text(size=7),
                                       axis.title.y=element_text(size=7),
                                       axis.title.x=element_blank(),
                                       legend.position="none",
                                       legend.text = element_text(size=7),
                                       legend.key = element_blank(),
                                       plot.margin=unit(c(0.1,0.1,0,0.1),"cm"),
                                       plot.title=element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(), 
                                       legend.title=element_blank(),
                                       plot.background = element_blank(),
                                       axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                       axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))  
  
  ecDNA<-data.frame(start=c(158998165,159679867,163146300,168035674),end=c(159678918,162720504,165697799,168056404))
  SVsnow<- chromothripsis@detail$SV
  SVsnow <- unique(SVsnow[SVsnow$chrom1 == y, ])
  min_coord = min(SVsnow$pos1)
  max_coord = max(SVsnow$pos2)
  ecDNA_plot<-ggplot() +
    geom_segment(data=ecDNA, aes(x = start, y =2 , xend = end, yend = 2),size=2,color="red")+
    common_ggplot2 + ylab("ecDNA")+ xlab(NULL)+
    scale_y_continuous(limits = c(1,3),breaks = c(1,2,3))+
    scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")},limits=c(min_coord,max_coord))
  
  ecDNA_plot<-ecDNA_plot+theme(plot.margin=unit(c(0,0.5,0.5,0.15), "cm"),legend.position="none")
  
  plot = arrangeGrob(plots[[1]],plots[[2]],plots[[3]],ecDNA_plot,nrow=4,ncol=1,heights=c(0.2,.4,.4,.2))
  plot
}

plots<-huatu("80","1")


pdf("sample80.ecDNA.chromothripisis.pdf",width =5.69 ,height = 4.28)
plot(plots)
dev.off()

#################### get many samples' results
get_shatter<-function(x){
  sv = read.table(stringr::str_glue("SV/allmerges/merged.{x}.txt"),header=T)# sv_matrix
  sv$chrom1=sub("chr","",sv$chrom1)
  sv$chrom2=sub("chr","",sv$chrom2)
  a=unique(c(which(!sv$chrom1%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")),
             which(!sv$chrom2%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"))))
  if(length(a)>0){
    sv2 = sv[-a,]
    sv = sv2
  }
  cn = read_tsv(stringr::str_glue("CNV/cns/{x}T.cs.rmdup.sort.call.cns"))
  cn <-as.data.frame(cn)
  cn = cn[,c(1,2,3,6)]
  colnames(cn)=c("chromosome","start","end","total_cn")
  b=unique(c(which(cn$chromosome=="chrY"),which(cn$chromosome=="chrM")))
  if(length(b)>0){
    cn2 = cn[-b,]
    cn = cn2
  }
  cn$chromosome=sub("chr","",cn$chromosome)
  CN_data = CNVsegs(chrom=as.character(cn$chromosome),start=cn$start,end=cn$end,total_cn=cn$total_cn)
  SV_data = SVs(chrom1=as.character(sv$chrom1),pos1=as.numeric(sv$start1),chrom2=as.character(sv$chrom2),pos2=as.numeric(sv$end2),SVtype=as.character(sv$svclass),strand1=as.character(sv$strand1),strand2=as.character(sv$strand2))
  chromothripsis = shatterseek(SV.sample=SV_data,seg.sample=CN_data,genome = "hg38")
  b = chromothripsis@chromSummary
  hc_index1 = which((b$clusterSize_including_TRA-b$number_TRA)>6 & b$max_number_oscillating_CN_segments_2_states>7 & b$pval_fragment_joins<0.05 & (b$chr_breakpoint_enrichment<0.05 | b$pval_exp_cluster <0.05))
  hc_index2 = which((b$clusterSize_including_TRA-b$number_TRA)>3 & b$number_TRA>=4 & b$max_number_oscillating_CN_segments_2_states>7 & b$pval_fragment_joins<0.05)
  hc_index = unique(c(hc_index1,hc_index2))
  lc_index = which((b$clusterSize_including_TRA-b$number_TRA)>6 & b$max_number_oscillating_CN_segments_2_states_chr<7 & b$max_number_oscillating_CN_segments_2_states_chr>3 & b$pval_fragment_joins<0.05 & (b$chr_breakpoint_enrichment<0.05 | b$pval_exp_cluster <0.05))
  b_hc = b[hc_index,]
  b_lc = b[lc_index,]
  b_hc$inter_other_chroms_coords_all=gsub('\n',';',b_hc$inter_other_chroms_coords_all)
  b_lc$inter_other_chroms_coords_all=gsub('\n',';',b_lc$inter_other_chroms_coords_all)
  write.table(file=stringr::str_glue("SV/chromothripsis/chromothripis_hc_{x}.matrix"),b_hc,sep="\t",quote=FALSE,row.names=F)
  write.table(file=stringr::str_glue("SV/chromothripsis/chromothripis_lc_{x}.matrix"),b_lc,sep="\t",quote=FALSE,row.names=F)
}

sv_samples<-c("1","2","3","4","5","6","7","8","9","10",
              "11","12","13","14","15","16","17","18","19",
              "20","21","23","24","25","26","27","28","29",
              "30","31","32","33","34","35","36","37","38",
              "39","40","41","42","43","44","45","46","47",
              "48","50","51","52","54","55","56","57","58",
              "59","60","62","63","64","65","66","67","68",
              "69","70","71","72","74","75","76","77","78",
              "79","80","84","85","86","87","88")


lapply(sv_samples,get_shatter)

###################merge chromothris from all samples
mg<-function(x){
  df<-read_tsv(stringr::str_glue("SV/chromothripsis/chromothripis_hc_{x}.matrix"))
  df<-df%>%
    mutate(sample=x)%>%
    select(sample,everything())
  df
}
mg2<-function(x){
  df<-read_tsv(stringr::str_glue("SV/chromothripsis/chromothripis_lc_{x}.matrix"))
  df<-df%>%
    mutate(sample=x)%>%
    select(sample,everything())
  df
}
fin<-do.call("rbind",lapply(sv_samples,mg))
fin2<-do.call("rbind",lapply(sv_samples,mg2))
rsult<-list(high_conf=fin,low_conf=fin2)
openxlsx::write.xlsx(rsult,"SV/chromothripsis/total.chromothripsis.xlsx")
##statistics
fin<-openxlsx::read.xlsx("SV/chromothripsis/total.chromothripsis2.xlsx",sheet = 1)
fin2<-openxlsx::read.xlsx("SV/chromothripsis/total.chromothripsis2.xlsx",sheet = 2)
AMPs<-openxlsx::read.xlsx("下机信息（含barcode）/all_sum.xlsx")
AMPs%>%
  as_tibble()%>%
  select(sample,AmpliconID,AmpliconType,Intervals)%>%
  tidyr::separate_rows(Intervals,sep=",")%>%
  mutate(name=paste(sample,AmpliconID,AmpliconType,sep=":"))%>%
  tidyr::separate(Intervals,into=c("chrom","start","end"),sep="[:-]")%>%
  select(chrom,start,end,name)%>%
  arrange(chrom,start)%>%
  write_tsv("SV/statistic/AMP_intervals.bed",col_names = F)

fin%>%
  select(chrom,start,end,sample)%>%
  mutate(chrom=paste0("chr",chrom))%>%
  write_tsv("SV/statistic/chromothripisis_hc.bed",col_names = F)

fin2%>%
  select(chrom,start,end,sample)%>%
  mutate(chrom=paste0("chr",chrom))%>%
  write_tsv("SV/statistic/chromothripisis_lc.bed",col_names = F)

df1<-read_tsv("SV/statistic/AMP_with_HC.bed",col_names = F)
df1<-df1%>%
  tidyr::separate(X4,into = c("sample","amp","type"),sep = ":")
df1$X8<-as.character(df1$X8)
df1<-df1%>%
  filter(sample==X8)

df2<-read_tsv("SV/statistic/AMP_with_LC.bed",col_names = F)
df2<-df2%>%
  tidyr::separate(X4,into = c("sample","amp","type"),sep = ":")
df2$X8<-as.character(df2$X8)
df2<-df2%>%
  filter(sample==X8)
totalCounts<-AMPs%>%group_by(sample,AmpliconType)%>%summarise(count=n())%>%ungroup()
hc_counts<-df1%>%select(sample,amp,type)%>%distinct(sample,amp,type,.keep_all=T)%>%group_by(sample,type)%>%summarise(hcCount=n())%>%
  dplyr::rename(AmpliconType=type)%>%
  ungroup()
lc_counts<-df2%>%select(sample,amp,type)%>%distinct(sample,amp,type,.keep_all=T)%>%group_by(sample,type)%>%summarise(lcCount=n())%>%
  dplyr::rename(AmpliconType=type)%>%
  ungroup()

fin<-totalCounts%>%left_join(hc_counts,by=c("sample","AmpliconType"))%>%
  left_join(lc_counts,by=c("sample","AmpliconType"))

fin<-fin%>%
  tidyr::replace_na(replace = list(hcCount=0,lcCount=0))

fin<-fin%>%
  mutate(hcRatio=hcCount/count,lcRatio=lcCount/count)
write_tsv(fin,"chrom_ecDNA_statis.tsv")
######### another chromothripis
get_cnv<-function(x){
  cn<-read_tsv(stringr::str_glue("CNV/annotate/{x}T.annotate.cns"))
  cn%>%
    select(1,2,3,7,5)%>%
    mutate(sample=x)%>%
    select(sample,everything())
}

sv_samples<-c("1","2","3","4","5","6","7","8","9","10",
              "11","12","13","14","15","16","17","18","19",
              "20","21","23","24","25","26","27","28","29",
              "30","31","32","33","34","35","36","37","38",
              "39","40","41","42","43","44","45","46","47",
              "48","50","51","52","54","55","56","57","58",
              "59","60","62","63","64","65","66","67","68",
              "69","70","71","72","74","75","76","77","78",
              "79","80","84","85","86","87","88")

cnvs<-do.call("rbind",lapply(sv_samples,get_cnv))
#anno gene affacted by breakpoint####################################################
get_sv<-function(x){
  sv<-read_tsv(stringr::str_glue("SV/allmerges/merged.{x}.txt"))
  sv%>%
    dplyr::select(1,2,9,4,5,10,11,12)%>%
    mutate(sample=x)%>%
    dplyr::select(sample,everything())
}
svs<-do.call("rbind",lapply(sv_samples,get_sv))
svs<-svs%>%
  #mutate(svclass=if_else(svclass=="h2hINV"|svclass=="t2tINV","INV",svclass))%>%
  mutate(svclass=if_else(abs(start2-start1)==1,"INS",svclass))

#insertions<-svs%>%
#  filter(svclass=="INS")

#others<-svs%>%
#  filter(svclass!="INS")

part1<-svs%>%
  select(1:3,8)%>%
  mutate(start=start1-1,end=start1)%>%
  select(chrom1,start,end,sample,svclass)
part2<-svs%>%
  select(1,5,6,8)%>%
  mutate(start=start2-1,end=start2)%>%
  select(chrom2,start,end,sample,svclass)%>%
  dplyr::rename(chrom1=chrom2)

allbreaks<-rbind(part1,part2)
haha<-split(allbreaks,allbreaks$sample)

lapply(names(haha),function(x){write_tsv(haha[[x]],stringr::str_glue("SV_breakpoints/{x}.bp.bed"),col_names = F)})

allbreaks%>%
  arrange(chrom1,start)%>%
  write_tsv("Breaks.detected.sv.bed",col_names = F)
############readin anno data
annos<-read_tsv("Breaks.detected.sv.anno.bed",col_names = F)
affected_genes<-annos%>%
  select(X4,X9)%>%
  dplyr::rename(samples=X4,genes=X9)%>%
  group_by(samples)%>%
  distinct(genes,.keep_all = T)%>%
  ungroup()%>%
  group_by(genes)%>%
  summarise(sampleN=n(),farc=sampleN/80*100)
gene_pos<-read_tsv("../Bladder/dbs/hg38.coding.bed",col_names = F)
gene_pos<-gene_pos%>%
  rowwise()%>%
  mutate(BP=round(median(seq(X2,X3))))%>%
  ungroup()%>%
  select(X4,X1,BP)%>%
  dplyr::rename(genes=X4,chrm=X1)
#saveRDS(gene_pos,"gene_coordinates.rds")
gene_pos<-readRDS("gene_coordinates.rds")

test_result<-affected_genes%>%
  left_join(gene_pos,by="genes")%>%
  select(genes,chrm,BP,sampleN)%>%
  dplyr::rename(SNP=genes,CHR=chrm,p=sampleN)


highLigh_genes<-test_result%>%
  filter(p>8)%>%
  pull(SNP)
#### huatu shuju
chrom.size<-read_tsv("../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size",col_names = F)
names(chrom.size)<-c("chromName","chromlength")
chrom.size<-chrom.size[1:24,]
chrom.size<-chrom.size%>%
  mutate(chromName=forcats::fct_relevel(chromName,c(paste0("chr",seq(1,22)),"chrX","chrY")))%>%
  arrange(chromName)

chrom.size<-chrom.size%>%
  mutate(Chromosome=stringr::str_remove(chromName,"chr"),
         chromlengthCumsum=cumsum(chromlength),
         chromStartFrom0=c(0,cumsum(chromlength)[-24]),
         chromMiddlePosFrom0=chromStartFrom0+chromlength/2)
db<-chrom.size%>%select(1,5)%>%dplyr::rename(CHR=chromName)
test_result<-test_result%>%left_join(db,by="CHR")
test_result<-test_result%>%mutate(start=BP+chromStartFrom0)
test_result$CHR<-factor(test_result$CHR,levels = c(paste0("chr",seq(1,22)),"chrX","chrY"))
highLigh_genes1<-test_result%>%
  filter(p>8)

commonGenes<-intersect(highLigh_genes1$SNP,highLigh_genes2$SNP)

highLigh<-test_result%>%
  filter(SNP%in%commonGenes)

test_result<-test_result%>%
  mutate(cols = if_else(CHR%in%col_group1,"col1","col2"),
         cols=if_else(p>8,"col3",cols))

ggplot(test_result,aes(x=start,y=p))+
  geom_point(aes(color = cols))+
  geom_text_repel(aes(label=SNP),data=highLigh,size=2.5,arrow = arrow(length = unit(0.02, "npc")),min.segment.length=0)+
  scale_x_continuous(label = chrom.size$chromName, breaks = chrom.size$chromMiddlePosFrom0)+
  scale_color_manual(values = c("col1"="#252525","col2"="#737373","col3"="#e41a1c"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5),legend.position = "none")
ggsave("SV_manhadun1.pdf",width = 8.1,height = 2.44)  

#######################################################################
###############################################epm diffs########################
epms<-readRDS("../Bladder/gene_drops.RDS")
epms<-epms%>%
  dplyr::select(1:81)%>%
  tidyr::gather(sample,abundance,-gene)%>%
  tidyr::replace_na(replace = list(abundance=0))%>%
  mutate(sample=stringr::str_remove(stringr::str_remove(sample,"cBca_"),"T"))%>%
  mutate(log2ab=log2(abundance+1))

allSamples<-unique(epms$sample)
getGroups<-annos%>%
  group_by(X9)%>%
  summarise(sampleCluster=list(unique(X4)))

getPvalue<-function(x){
  tmp<-epms%>%
    filter(gene==x)
  gp1<-(getGroups%>%
          filter(X9==x)%>%pull(sampleCluster))[[1]]
  gp2<-setdiff(allSamples,gp1)
  exp1<-tmp%>%filter(sample%in%gp1)%>%pull(abundance)
  exp2<-tmp%>%filter(sample%in%gp2)%>%pull(abundance)
  wocao<-wilcox.test(exp1,exp2)
  list(p=wocao$p.value,fc=log2(mean(exp1)/mean(exp2)))
  
}

newDIFF<-test_result%>%
  filter(p>1)%>%
  filter(SNP%in%(epms$gene))%>%
  select(SNP,CHR,BP)%>%
  mutate(P=purrr::map(SNP,function(x){getPvalue(x)}))

newDIFF<-newDIFF%>%
  mutate(pvalue=purrr::map_dbl(P,function(x){x$p}),
         logfc=purrr::map_dbl(P,function(x){x$fc}))%>%
  select(-P)

highLigh_genes2<-newDIFF%>%
  filter(pvalue<0.05)
highLigh2<-newDIFF%>%
  filter(SNP%in%commonGenes)  

newDIFF<-newDIFF%>%
  left_join(db,by="CHR")

newDIFF$CHR<-factor(newDIFF$CHR,levels = c(paste0("chr",seq(1,22)),"chrX","chrY"))
newDIFF<-newDIFF%>%
  mutate(start=BP+chromStartFrom0)
newDIFF<-newDIFF%>%
  mutate(p=log10(pvalue))
col_group1<-c("chr1","chr3","chr5","chr7","chr9","chr11","chr13","chr15","chr17",
              "chr19","chr21","chrX")

newDIFF<-newDIFF%>%
  mutate(cols = if_else(CHR%in%col_group1,"col1","col2"),
         cols=if_else(pvalue<0.05,"col3",cols))


ggplot(newDIFF,aes(x=start,y=p))+
  geom_point(aes(color = cols))+
  geom_text_repel(aes(label=SNP),data=highLigh2,arrow = arrow(length = unit(0.02, "npc")),min.segment.length=0)+
  scale_x_continuous(label = chrom.size$chromName, breaks = chrom.size$chromMiddlePosFrom0)+
  scale_color_manual(values = c("col1"="#252525","col2"="#737373","col3"="#e41a1c"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5),legend.position = "none")
ggsave("SV_manhadun2.pdf",width = 8.1,height = 2.44)  

##200genes
geneid<-bitr(haha1,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
eg<-enrichGO(geneid$ENTREZID,OrgDb = "org.Hs.eg.db",ont = "BP",pvalueCutoff = 0.05,readable = T)
ek<-enrichKEGG(geneid$ENTREZID,organism = "hsa",pvalueCutoff = 0.05)
#######################################################################
###########################################GSEA for circular vs. no-circular#####################################
sampleClass<-read_tsv("class_sample.tsv")

allresult<-read_tsv("diffGenes_ecc_noecc.tsv")
geneNames = bitr(allresult$Geneid,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
names(geneNames)<-c("Geneid","id")
allresult<-allresult%>%inner_join(geneNames,by="Geneid")
geneList = allresult$log2FoldChange
names(geneList) = allresult$id
geneList = sort(geneList, decreasing = TRUE)

##download GSEA db
sub1<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:BIOCARTA")
sub2<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
sub3<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:PID")
sub4<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
C2_t2g<-rbind(sub1,sub2,sub3,sub4)%>%
  dplyr::select(gs_name, entrez_gene)
C2_t2n<-rbind(sub1,sub2,sub3,sub4)%>%
  dplyr::select(gs_name,gs_description)

em2 <- GSEA(geneList, TERM2GENE = C2_t2g,TERM2NAME=C2_t2n,minGSSize = 20,eps = 0)
saveRDS(em2,"gsea.result.ecDNA_vs_none.rds")
##################################################################################################
##########################diffGene_with_eccDNA(higher and lower)#####################################################
counts<-read_tsv("RNAseq/allSample.count.txt")
counts<-counts%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.data.frame()
counts<-counts[,1:70]
groupInfo<-read_tsv("../Bladder/EPM and clinical variables.tsv")
groupInfo<-groupInfo%>%
  select(1,4)
colData<-tibble(Sample = colnames(counts))
colData<-colData%>%left_join(groupInfo,by="Sample")
colData<-colData%>%
  mutate(group=if_else(`EPM Group`=="Low","A","B"))%>%
  select(Sample,group)%>%
  as.data.frame()
colData$group<-factor(colData$group,levels = c("A","B"))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group) 
dds <- dds[rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds)
res <- results(dds)
allresult<-as.data.frame(res)%>%
  tibble::rownames_to_column(var="Geneid")%>%
  as_tibble()

allresult<-allresult%>%
  mutate(logp=-log10(padj))

allresult<-allresult%>%
  mutate(col=case_when((log2FoldChange<(-1))&(logp>2)~"downregular",
                       (log2FoldChange>1)&(logp>2)~"upregular",
                       TRUE ~ "none"))

##############################GSVA#############################################
TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
TPMs<-TPMs%>%dplyr::select(1:71)
TPMs<-TPMs%>%
  filter(Geneid%in%(allresult$Geneid))
TPMs<-TPMs%>%
  tibble::column_to_rownames(var="Geneid")%>%
  as.matrix()
logTPM<-log2(TPMs+1)
geneSet<-getGmt("RNAseq/c2.cp.kegg.v7.5.1.symbols.gmt")
fin<-gsva(logTPM,geneSet,min.sz=20,max.sz=500)

mod <- model.matrix(~ 0+factor(colData$group))
colnames(mod) <- c("A", "B")
fit <- lmFit(fin, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.05)
tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]
plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75),
     main="", xlab="GSVA enrichment score difference", ylab=expression(-log[10]~~Raw~P-value))
abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])), col=grey(0.5), lwd=1, lty=2)
points(tt$logFC[match(DEpwys, rownames(tt))],
       -log10(tt$P.Value[match(DEpwys, rownames(tt))]), pch=".", cex=5, col="darkred")
text(max(tt$logFC)*0.85, -log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])), "1% FDR", pos=3)
fin%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var="Geneset")%>%
  as_tibble()%>%
  rowwise()%>%
  mutate(mean1=mean(c_across(2:71)),mean2=mean(c_across(72:127)))%>%
  ungroup()%>%
  dplyr::select(Geneset,mean1,mean2)

####################################GSEA#################################################

geneNames = bitr(allresult$Geneid,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
names(geneNames)<-c("Geneid","id")
allresult<-allresult%>%inner_join(geneNames,by="Geneid")
geneList = allresult$log2FoldChange
names(geneList) = allresult$id
geneList = sort(geneList, decreasing = TRUE)

##download GSEA db
sub1<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:BIOCARTA")
sub2<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
sub3<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:PID")
sub4<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
C2_t2g<-rbind(sub1,sub2,sub3,sub4)%>%
  dplyr::select(gs_name, entrez_gene)
C2_t2n<-rbind(sub1,sub2,sub3,sub4)%>%
  dplyr::select(gs_name,gs_description)

em2 <- GSEA(geneList, TERM2GENE = C2_t2g,TERM2NAME=C2_t2n,minGSSize = 20,eps = 0)
enrichplot::gseaplot2(em2, geneSetID = c(6,9,11,20,22,26,30,33,34,38,43,58),subplots = 2:3,rel_heights = c(1,1))

rankvalues<-tibble(geneid=names(rev(geneList)),foldchange=rev(geneList),order=1:18680)
needID = c(6,9,11,20,22,26,30,33,34,38,43,58)
getPlotData<-function(i){
  tibble(geneid=em2@geneSets[[em2@result[needID[i],]$ID]])%>%
    left_join(rankvalues,by="geneid")%>%
    mutate(y1=i-1,y2=i-1+0.9,col=if_else(foldchange>0,"yes","no"))%>%
    tidyr::drop_na()
}

getPlotData2<-function(i){
  tibble(geneid=stringr::str_split(em2@result[needID[i],]$core_enrichment,"/")[[1]])%>%
    left_join(rankvalues,by="geneid")%>%
    mutate(y1=i-1,y2=i-1+0.9,col=if_else(foldchange>0,"yes","no"))%>%
    tidyr::drop_na()
}

haha<-do.call("rbind",lapply(1:12, getPlotData)) 
ggplot(haha,aes(x=order,y=y1))+
  geom_segment(aes(xend=order,yend=y2,color=col))+
  scale_color_manual(values = c("no"="#5B799D","yes"="#CC726A"))+
  theme_void()

ggsave("GSEA_result.part3.png",width = 7.04,height = 2.37,dpi=300)

rankvalues<-rankvalues%>%
  mutate(col=if_else(foldchange>0,"yes","no"))

ggplot(rankvalues,aes(x=order,y=foldchange))+
  geom_area(aes(fill=col))+
  scale_fill_manual(values = c("no"="#5B799D","yes"="#CC726A"))+
  theme_pubr()
ggsave("GSEA_result.part2.pdf",width =7.04 ,height = 4.14)


ggplot(rankvalues,aes(x=foldchange))+
  geom_histogram(fill="#5B799D",color="black",bins = 100,size=0.2)+
  scale_x_continuous(limits = c(-2,2))+
  theme(legend.position = "none")+
  theme_pubr()+ylab("Gene Density")+
  xlab("Correlation")
ggsave("GSEA_result.part1.pdf",width =9.44,height = 4.14)

write_tsv(em2@result,"GSEA_result.txt")
##################################################################################################
########################################sv statistic##############################################
get_sv_numbers<-function(x){
  df<-read_tsv(stringr::str_glue("SV_breakpoints/{x}.break.bed"),col_names = F)
  df<-df%>%
    mutate(haha=paste(X1,X2,X3,sep=":"))%>%
    distinct(haha,.keep_all = T)%>%
    count(X4)
  df$sample<-x
  df
}

samples_sv<-stringr::str_remove(stringr::str_remove(Sys.glob("SV_breakpoints/*.break.bed"),".break.bed"),"SV_breakpoints/")

svnumbers<-do.call("rbind",lapply(samples_sv,get_sv_numbers))


get_ecc_breaks<-function(x){
  df<-read_tsv(stringr::str_glue("ecc_inters_sv/{x}.inter.bed"),col_names = F)
  tmp<-df%>%
    mutate(haha=paste(X1,X2,X3,sep=":"))%>%
    distinct(haha,.keep_all = T)%>%
    count(X4)
  tmp$sample<-x
  tmp
}

samples<-stringr::str_remove(stringr::str_remove(Sys.glob("ecc_inters_sv/*.inter.bed"),".inter.bed"),"ecc_inters_sv/")
fin<-do.call("rbind",lapply(samples,get_ecc_breaks))
names(fin)<-c("type","numbers","samples")

fin<-fin%>%group_by(type)%>%summarise(totaln=sum(numbers))
fin$totalsv<-sum(svnumbers$svnumbers)

get_ec_breaks<-function(x){
  df<-read_tsv(stringr::str_glue("ec_inters_sv/{x}.inter.bed"),col_names = F)
  tmp<-df%>%
    mutate(haha=paste(X1,X2,X3,sep=":"))%>%
    distinct(haha,.keep_all = T)%>%
    count(X4)
  tmp$sample<-x
  tmp
}

samples<-stringr::str_remove(stringr::str_remove(Sys.glob("ec_inters_sv/*.inter.bed"),".inter.bed"),"ec_inters_sv/")
fin2<-do.call("rbind",lapply(samples,get_ec_breaks))
names(fin2)<-c("type","numbers","samples")
fin2<-fin2%>%group_by(type)%>%summarise(totaln=sum(numbers))
fin2$totalsv<-sum(svnumbers$svnumbers)

get_ecc_breaks<-function(x){
  df<-read_tsv(stringr::str_glue("ecc_inters_sv/{x}.inter.random.bed"),col_names = F)
  tmp<-df%>%
    mutate(haha=paste(X1,X2,X3,sep=":"))%>%
    distinct(haha,.keep_all = T)%>%
    count(X5)
  tmp$sample<-x
  tmp
}

samples<-stringr::str_remove(stringr::str_remove(Sys.glob("ecc_inters_sv/*.inter.random.bed"),".inter.random.bed"),"ecc_inters_sv/")
fin3<-do.call("rbind",lapply(samples,get_ecc_breaks))
names(fin3)<-c("type","numbers","samples")
fin3<-fin3%>%group_by(type)%>%summarise(totaln=sum(numbers))
fin3$totalsv<-sum(svnumbers$svnumbers)

get_ecc_breaks<-function(x){
  df<-read_tsv(stringr::str_glue("ecc_inters_sv/{x}.inter.ecWithecc.bed"),col_names = F)
  tmp<-df%>%
    mutate(haha=paste(X1,X2,X3,sep=":"))%>%
    distinct(haha,.keep_all = T)%>%
    count(X4)
  tmp$sample<-x
  tmp
}

samples<-stringr::str_remove(stringr::str_remove(Sys.glob("ecc_inters_sv/*.inter.ecWithecc.bed"),".inter.ecWithecc.bed"),"ecc_inters_sv/")
fin4<-do.call("rbind",lapply(samples,get_ecc_breaks))
names(fin4)<-c("type","numbers","samples")
fin4<-fin4%>%group_by(type)%>%summarise(totaln=sum(numbers))
fin4$totalsv<-sum(svnumbers$svnumbers)

lst<-list(ecc_inters=fin,ec_inters=fin2,randoms=fin3)
openxlsx::write.xlsx(lst,"breaks_with_ec.xlsx")

fin<-openxlsx::read.xlsx("breaks_with_ec.xlsx",sheet = 1)
fin2<-openxlsx::read.xlsx("breaks_with_ec.xlsx",sheet = 2)
fin3<-openxlsx::read.xlsx("breaks_with_ec.xlsx",sheet = 3)
fin_sum<-fin%>%
  group_by(type)%>%
  summarise(totaln=sum(numbers),totalsv=sum(svnumbers))%>%
  mutate(pct=totaln/totalsv)
fin2_sum<-fin2%>%
  group_by(type)%>%
  summarise(totaln=sum(numbers),totalsv=sum(svnumbers))%>%
  mutate(pct=totaln/totalsv)
fin3_sum<-fin3%>%
  group_by(type)%>%
  summarise(totaln=sum(numbers),totalsv=sum(svnumbers))%>%
  mutate(pct=totaln/totalsv)
##################################################################################################
####################################some genes####################################################
genes<-c("LIG3","LIG4","POLM","POLQ","PRKDC","BRCA1","BRCA2","MSH3")
##exps
TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
TPM2<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)
TPM3<-TPM2%>%filter(stringr::str_ends(sample,"T"))
TPM3$sample<-stringr::str_remove(TPM3$sample,"_T")

##group
groupInfo<-read_tsv("class_sample.tsv")
names(groupInfo)<-c("sample","type")
groupInfo<-groupInfo%>%
  mutate(type=if_else(type=="ecDNA","ecDNA","others"))
TPM3<-TPM3%>%
  left_join(groupInfo,by="sample")

EPMS<-openxlsx::read.xlsx("../Bladder/Bladder_circle_numbers.xlsx")%>%
  as_tibble()
EPMs<-EPMS%>%
  filter(group=="Tumour")%>%
  dplyr::select(sample,EPM)
EPMs$sample<-stringr::str_remove(EPMs$sample,"c")
EPMs$sample<-stringr::str_remove(EPMs$sample,"T")
EPMs$sample<-stringr::str_remove(EPMs$sample,"N")
TPM3<-TPM3%>%
  left_join(EPMs,by="sample")

TPM3<-TPM3%>%
  mutate(logT=log2(TPM+1),logE=log2(EPM+1))

goodGenes<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)%>%
  filter(stringr::str_ends(sample,"T"))%>%
  mutate(type=if_else(TPM>1,"good","bad"))%>%
  group_by(Geneid,type)%>%
  summarise(count=n())%>%
  tidyr::pivot_wider(names_from = "type",values_from = count,values_fill = 0)%>%
  filter(good==70)%>%
  pull(Geneid)

diffGenes<-TPM3%>%
  dplyr::filter(Geneid%in%goodGenes)%>%
  dplyr::group_by(Geneid)%>%
  rstatix::wilcox_test(logT~type)

diffCors2<-TPM3%>%
  #dplyr::filter(Geneid%in%goodGenes)%>%
  dplyr::group_by(Geneid)%>%
  rstatix::cor_test(vars = "logT",vars2 = "logE",method = "spearman")

diffGene1<-diffGenes%>%
  filter(p<0.05)
diffCor1<-diffCors%>%
  filter(p<0.05)

jiaoji2<-intersect(diffGene1$Geneid,diffCor1$Geneid)

zhuanhua<-bitr(jiaoji,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

ego<-enrichGO(zhuanhua2$ENTREZID,OrgDb =org.Hs.eg.db,pvalueCutoff = 0.05,ont = 'BP',readable = T )

eko<-enrichKEGG(zhuanhua2$ENTREZID,organism = "hsa",pvalueCutoff = 0.05)

rst<-list(intersectGenes=jiaoji,GOBP=ego@result,KEGG=eko@result)
openxlsx::write.xlsx(rst,"intersect.genes2.xlsx")

### heatmap for diffgenes
plotMat<-TPM3%>%
  dplyr::select(Geneid,sample,logT)

diffCors2<-diffCors2%>%
  mutate(tp=if_else(p<0.05,"yes","no"))

ggplot(diffCors2,aes(x=cor))+
  geom_histogram(aes(fill=tp),color="black",bins =105)+
  scale_fill_manual(values = c("yes"="#08519c","no"="#feb24c"))+
  theme_pubr()+
  ylab("Gene count")+
  geom_vline(xintercept = 0.089,linetype="dashed",color="red")+
  xlab("Correlation")

ggsave("corplot2.pdf",width = 7.03,height = 3.02)

##box
exp_plot<-TPM3%>%
  filter(Geneid%in%genes)

p1<-ggplot(exp_plot%>%filter(Geneid==genes[1]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[1]} expression value"))+
  theme(legend.position = "none")

p2<-ggplot(exp_plot%>%filter(Geneid==genes[2]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[2]} expression value"))+
  theme(legend.position = "none")

p3<-ggplot(exp_plot%>%filter(Geneid==genes[3]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[3]} expression value"))+
  theme(legend.position = "none")

p4<-ggplot(exp_plot%>%filter(Geneid==genes[4]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[4]} expression value"))+
  theme(legend.position = "none")

p5<-ggplot(exp_plot%>%filter(Geneid==genes[5]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[5]} expression value"))+
  theme(legend.position = "none")

p6<-ggplot(exp_plot%>%filter(Geneid==genes[6]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[6]} expression value"))+
  theme(legend.position = "none")

p7<-ggplot(exp_plot%>%filter(Geneid==genes[7]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[7]} expression value"))+
  theme(legend.position = "none")

p8<-ggplot(exp_plot%>%filter(Geneid==genes[8]),aes(x=type,y=logT))+
  geom_boxplot(aes(fill=type),width=0.7)+
  scale_fill_manual(values = c("others"="#FBF3AA","ecDNA"="#F3B854"))+
  xlab("")+
  ylab(stringr::str_glue("{genes[8]} expression value"))+
  theme(legend.position = "none")

ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol = 4,nrow=2,align = "v")

ggsave("somegenes.plot1.pdf",width = 6.85,height = 5.3)
###xianggaun
p1<-ggplot(exp_plot%>%filter(Geneid==genes[1]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[1]} eccDNA value"))

p2<-ggplot(exp_plot%>%filter(Geneid==genes[2]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[2]} eccDNA value"))

p3<-ggplot(exp_plot%>%filter(Geneid==genes[3]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[3]} eccDNA value"))

p4<-ggplot(exp_plot%>%filter(Geneid==genes[4]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[4]} eccDNA value"))

p5<-ggplot(exp_plot%>%filter(Geneid==genes[5]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[5]} eccDNA value"))

p6<-ggplot(exp_plot%>%filter(Geneid==genes[6]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[6]} eccDNA value"))

p7<-ggplot(exp_plot%>%filter(Geneid==genes[7]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[7]} eccDNA value"))

p8<-ggplot(exp_plot%>%filter(Geneid==genes[8]),aes(x=logT,y=logE))+
  geom_point(color="#F3B854",size=2.1)+
  stat_smooth(geom = "line",method = "lm")+
  xlab("Expression value")+
  ylab(stringr::str_glue("{genes[8]} eccDNA value"))

ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol = 4,nrow=2,align = "v")
ggsave("somegenes.plot2.pdf",width = 6.85,height = 5.3)
##################################################################################################
###################################相关性通路#####################################################
eg1<-openxlsx::read.xlsx("corgenes.pathway.xlsx",sheet = 1)
ek1<-openxlsx::read.xlsx("corgenes.pathway.xlsx",sheet = 2)
eg2<-openxlsx::read.xlsx("corgenes.pathway.xlsx",sheet = 3)
ek2<-openxlsx::read.xlsx("corgenes.pathway.xlsx",sheet = 4)
my_cols <- c(
  "#fee5d9",
  "#fcbba1",
  "#fc9272",
  "#fb6a4a",
  "#de2d26",
  "#a50f15")


tmp<-rbind(eg2,ek2)%>%
  as_tibble()%>%
  arrange(p.adjust)

tmp<-tmp[1:30,]



ggballoonplot(tmp, x="Count",y="Description",fill = "pvalue",size="Count")+
  scale_fill_gradientn(colors = my_cols)
ggsave("normal_related.pathway.pdf",width =9 ,height = 9)
###################################################################
###################plot for gistic score###########################
chrom.size<-read_tsv("../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size",col_names = F)
names(chrom.size)<-c("chromName","chromlength")
chrom.size<-chrom.size[1:24,]
chrom.size<-chrom.size%>%
  mutate(chromName=forcats::fct_relevel(chromName,c(paste0("chr",seq(1,22)),"chrX","chrY")))%>%
  arrange(chromName)

chrom.size<-chrom.size%>%
  mutate(Chromosome=stringr::str_remove(chromName,"chr"),
         chromlengthCumsum=cumsum(chromlength),
         chromStartFrom0=c(0,cumsum(chromlength)[-24]),
         chromMiddlePosFrom0=chromStartFrom0+chromlength/2)
chrom.size<-chrom.size[1:22,]
chrom.size$Chromosome<-as.double(chrom.size$Chromosome)
chrom.size$ypos<-rep(c(1.5,1.75),11)
scores<-read_tsv("CNV/gistic2_results/scores.gistic")
scores<-scores%>%
  left_join(chrom.size,by="Chromosome")%>%
  mutate(Start=Start+chromStartFrom0,
         End=End+chromStartFrom0,
         `G-score`=if_else(Type=="Del",0-`G-score`,`G-score`))

ggplot(scores,aes(x=Start,y=`G-score`))+
  geom_area(aes(group=Type,fill=factor(Type,levels = c("Del","Amp"))))+
  geom_vline(aes(xintercept=chromlengthCumsum),linetype=2)+
  #geom_text(aes(x=chromMiddlePosFrom0,y=ypos,label=chromName),fonsize=1.2,fontface="italic")+
  scale_fill_manual(values = c("Del"="#325189","Amp"="#C8493A"))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank())
ggsave("gistic_results.png",width = 10.21,height = 2.40,dpi = 300)
###################################################################

###########################################Urine############################################
samples<-c("1","2","5","8","9","13",
           "31","34","36","37",
           "50","60","67","69","70",
           "72","74","75","86")

get_cal<-function(x){
  if (file.exists(stringr::str_glue("CCGA-Urine-AA/CCGA-Urine/CCGA-Urine-AC-results/CCGA-Urine-AC-results/ecDNA_inters/jiaoji.{x}.bed"))){
      df<-read_tsv(stringr::str_glue("CCGA-Urine-AA/CCGA-Urine/CCGA-Urine-AC-results/CCGA-Urine-AC-results/ecDNA_inters/jiaoji.{x}.bed"),col_names = F)
      ecNUMbs<-df%>%
      distinct(X4,X5,X6,.keep_all = T)%>%
        nrow()
      total_numbers<-read_tsv(stringr::str_glue("ecDNA_beds/{x}.merge.ecDNA.bed"),col_names = F)%>%
        nrow()
      tibble(sample=x,inters=ecNUMbs,totals=total_numbers)
  }else{
    total_numbers<-read_tsv(stringr::str_glue("ecDNA_beds/{x}.merge.ecDNA.bed"),col_names = F)%>%
      nrow()
    tibble(sample=x,inters=0,totals=total_numbers)
  }
}

fin<-do.call('rbind',lapply(samples, get_cal))
write_tsv(fin,"Urine/Urine.inters.txt")
###########################################################################################
#############################################suilie ecDNA#############################################

df<-read_tsv("SV/suilie_ec_jiaoji.bed",col_names = F)
sample_samples<-df%>%
  filter(X4==X8)
totals<-read_tsv("SV/allsamples.locus.bed",col_names = F)

#for sample

df<-read_tsv("SV/suilie_ec_sample.txt",col_names = F)
#############################################################################################
#########################################mutate vs. normal###################################
genes<-c("AKAP6","ZFHX3","WNK1","CDH23","FGFR3","CALY","TTN")
##exps
TPMs<-read_tsv("RNAseq/Protein_coding.tpms.txt")
TPM2<-TPMs%>%
  tidyr::gather(sample,TPM,-Geneid)
TPM3<-TPM2%>%filter(stringr::str_ends(sample,"T"))
TPM3$sample<-stringr::str_remove(TPM3$sample,"_T")
##group
groupInfo<-read_tsv("RNAseq/UBC_Combine_NonSyn.maf")
plt<-function(x,db=groupInfo,exp=TPM3){
  cs <- db%>%
    filter(Hugo_Symbol==x)%>%
    pull(Tumor_Sample_Barcode)%>%
    stringr::str_remove("CCGA-UBC-")%>%
    stringr::str_remove("^0+")
  cs <- paste0("Bca_",cs)
  exp_plot<-exp%>%
    filter(Geneid==x)%>%
    mutate(gp=if_else(sample%in%cs,"Mutate","wideType"))%>%
    mutate(logp=log2(TPM+1))
  exp_plot$gp<-factor(exp_plot$gp,levels = c("wideType","Mutate"))
  p<-ggboxplot(exp_plot,
               x="gp",
               y="logp",xlab = "",ylab = "Expression",
               width=0.5,
               color = "black", 
               palette =c("#325189", "#C8493A"),
               add = "jitter",
               add.params=list(color="gp",size=1.2),
               nrow=2)+stat_compare_means(label.y = 7,bracket.size=1)+
    theme(strip.background = element_rect(fill="#FAF4E3"))
  p<-ggpar(p,legend.title = "",main = x)
  p
}

p1<-plt(genes[1])
p2<-plt(genes[2])
p3<-plt(genes[3])
p4<-plt(genes[4])
p5<-plt(genes[5])
p6<-plt(genes[6])
p7<-plt(genes[7])

ggarrange(p1,p2,p3,p4,p5,p6,p7,nrow = 2,ncol = 4)
ggsave("muatate_vs_widetp.wilcox.pdf",width = 9.44,height = 7.46)

#############################################################################################

###########################groupheatmap###########################
df<-openxlsx::read.xlsx("Figure S7B.xlsx")%>%as_tibble()
df<-df%>%
  arrange(EPM)
mat<-df%>%
  dplyr::select(1,2)%>%
  tibble::column_to_rownames(var="Sample.ID")%>%
  as.matrix()%>%
  t()
ht<- Heatmap(mat,
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = F,
            bottom_annotation = HeatmapAnnotation(Age=df$Age,
                                               Gender=df$Gender,
                                               MIBC=df$`NMIBC/MIBC`,
                                               Grade=df$Grade,
                                               ecDNA=df$ecDNA,
                                               Msig=df$MSig,
                                               PR=df$`Primary/Relapsed`,
                                               MPMT=df$MPMT,
                                               Prognosis=df$Prognosis,
                                               gp=gpar(col="white",lwd=0.5),
                                               col = list(Age=c("<=65"="#B9DFFB",">65"="#68C84D"),
                                                          Gender=c("Female"="#E93420","Male"="#316DBB"),
                                                          MIBC=c("MIBC"="#AE2417","NMIBC"="#F3F3F4"),
                                                          Grade=c("High"="#4FADEB","Low"="#F3F3F4"),
                                                          ecDNA=c("Positive"="#E5CFAB","Negative"="#A5C5CE"),
                                                          Msig=c("MSig1"="#9C362F","MSig2"="#E5CFAB","MSig3"="#A5C5CE"),
                                                          PR=c("Primary"="#B0AFB0","Relapesd"="#1B1819"),
                                                          MPMT=c("Yes"="#1B1819","No"="#B0AFB0"),
                                                          Prognosis=c("Good"="#E5CFAB","Poor"="#A5C5CE")
                                               )
                                               
            ),name = "Abundance")


pdf("testPDF.pdf",width = 11.86,height = 2.23)
draw(ht)
dev.off()

