
# cal baf from ASCAT ------------------------------------------------------

readAlleleCountFiles=function(files,minCounts,keep_chr_string=F) {
  data=do.call(rbind,lapply(files,function(x) {
    tmp=data.frame(data.table::fread(x,sep='\t',showProgress=F,header=T),stringsAsFactors=F)
    tmp=tmp[tmp[,7]>=minCounts,]
    if (nrow(tmp)>0) {
      if (!keep_chr_string) tmp[,1]=gsub('^chr','',tmp[,1])
      rownames(tmp)=paste0(tmp[,1],'_',tmp[,2])
    }
    return(tmp)
  }))
  stopifnot(nrow(data)>0)
  return(data)
}


readAllelesFiles <- function(prefix, suffix, chrom_names, add_chr_string = F) {
  files <- paste0(prefix, chrom_names, suffix)
  files <- files[sapply(files, function(x) file.exists(x) && file.info(x)$size > 0)]
  stopifnot(length(files) > 0)
  data <- do.call(rbind, lapply(files, function(x) {
    tmp <- data.frame(data.table::fread(x, sep = "\t", showProgress = F, header = T))
    tmp <- tmp[!is.na(tmp[, 2] & !is.na(tmp[, 3])), ]
    tmp <- tmp[!duplicated(tmp[, 1]), ]
    tmp$chromosome <- gsub(paste0(prefix, "(", paste(chrom_names, collapse = "|"), ")", suffix), "\\1", x)
    if (add_chr_string) tmp$chromosome <- paste0("chr", tmp$chromosome)
    tmp <- tmp[, c(4, 1:3)]
    rownames(tmp) <- paste0(tmp[, 1], "_", tmp[, 2])
    return(tmp)
  }))
  stopifnot(nrow(data) > 0)
  return(data)
}

getBAFs <- function(AlleleCountsFile, alleles.prefix,chrom_names=c(1:22,'X'), minCounts=20) {
  set.seed(123)
  # Load data, only keep SNPs with enough coverage
  tumour_input_data = readAlleleCountFiles(AlleleCountsFile,minCounts,keep_chr_string = T)
  allele_data = readAllelesFiles(alleles.prefix, ".txt", chrom_names)
  # Synchronise DFs
  matched_data = Reduce(intersect, list(rownames(tumour_input_data), rownames(allele_data)))
  tumour_input_data = tumour_input_data[rownames(tumour_input_data) %in% matched_data,]
  allele_data = allele_data[rownames(allele_data) %in% matched_data,]
  rm(matched_data)
  tumour_input_data = tumour_input_data[,3:6]
  # Obtain depth for both alleles for tumour and normal
  len = nrow(allele_data)
  tumour_input_data$REF=tumour_input_data[cbind(1:len,allele_data[,3])]
  tumour_input_data$ALT=tumour_input_data[cbind(1:len,allele_data[,4])]
  # Make sure that ALT+REF fit with minimal counts
  TO_KEEP=which(tumour_input_data$REF+tumour_input_data$ALT>=minCounts)
  allele_data=allele_data[TO_KEEP,]
  tumour_input_data=tumour_input_data[TO_KEEP,]
  rm(TO_KEEP)
  # Prepare allele counts to derive BAF and logR
  len = nrow(allele_data)
  mutCount1 = tumour_input_data$REF
  mutCount2 = tumour_input_data$ALT
  totalTumour = mutCount1 + mutCount2
  rm(tumour_input_data)
  tumourBAF = vector(length=len, mode="numeric")
  # Randomise A and B alleles
  selector = round(runif(len))
  tumourBAF[which(selector==0)] = mutCount1[which(selector==0)] / totalTumour[which(selector==0)]
  tumourBAF[which(selector==1)] = mutCount2[which(selector==1)] / totalTumour[which(selector==1)]
  # Create the output data.frames
  tumor.BAF = data.frame(Chromosome=allele_data$chromosome, Position=allele_data$position, baf=tumourBAF, ID=rownames(allele_data), row.names=4, stringsAsFactors=F)
}

##################################################################
library(GenomicRanges)
library(ggbio)
library(readr)

keyRegion<-GRanges(seqnames = "chr6",IRanges(start = 38433634,end =42493758))
CCND3<-GRanges(seqnames = "chr6",IRanges(start = 41934934,end = 42050357))
####################################WGS##################################################
WGS_cov <-read_tsv("WGS.cov.bdg",col_names=F)
names(WGS_cov)<-c("seqnames","start","end","coverage")
WGS_cov_gr<-makeGRangesFromDataFrame(WGS_cov,keep.extra.columns = T)

p1<-ggplot()+
  geom_rect(keyRegion,aes(xmin=start,xmax=end,ymin=0,ymax=300),
            stat = "identity",
            fill="#525152",
            alpha=0.1,
            color="#525152")+
  stat_identity(WGS_cov_gr,aes(x=start,y=coverage),
                fill="#525152",
                geom="area")+
  scale_y_continuous(breaks = c(0,300),limits = c(0,300))+
  theme_pack_panels()

wgs_baf<-getBAFs("WGS.allele.txt","G1000_alleles_hg38_",chrom_names = "chr6")

wgs_baf<-as_tibble(wgs_baf)%>%
  mutate(end=Position)%>%
  select(Chromosome,Position,end,baf)%>%
  setNames(c("seqnames","start","end","BAF"))

##we should filter baf==0 /1

wgs_baf<-wgs_baf%>%
  filter(BAF!=0 & BAF!=1)%>%
  mutate(col=case_when(BAF>=0.75 ~ "top",
                       BAF<=0.25 ~ "bot",
                       TRUE ~ "mid"))
keep_alleles<-paste(wgs_baf$seqnames,wgs_baf$start,sep = ":")
wgs_baf_gr<-makeGRangesFromDataFrame(wgs_baf,keep.extra.columns = T)

p2<-ggplot()+
  geom_rect(keyRegion,aes(xmin=start,xmax=end,ymin=0,ymax=1),
            stat = "identity",
            fill="#525152",
            alpha=0.1,
            color="#525152")+
  stat_identity(wgs_baf_gr,aes(x=start,y=BAF,color=col),
                geom="point",size=0.5)+
  scale_color_manual(values = c("top"="#E93326","bot"="#0100F3","mid"="#3F8C37"))+
  theme_pack_panels()+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(panel.grid.major.y = element_line(linetype = "dashed",color="grey"))

####################################WGS##################################################
####################################RNA##################################################
RNA_cov <-read_tsv("RNA.cov.bdg",col_names=F)
names(RNA_cov)<-c("seqnames","start","end","coverage")
RNA_cov<-RNA_cov%>%
  mutate(coverage=if_else(coverage>300,300,coverage))
RNA_cov_gr<-makeGRangesFromDataFrame(RNA_cov,keep.extra.columns = T)

p5<-ggplot()+
  geom_rect(keyRegion,aes(xmin=start,xmax=end,ymin=0,ymax=300),
            stat = "identity",
            fill="#525152",
            alpha=0.1,
            color="#525152")+
  stat_identity(RNA_cov_gr,aes(x=start,y=coverage),
                fill="#525152",
                geom="area")+
  scale_y_continuous(breaks = c(0,300),limits = c(0,300))+
  theme_pack_panels()

rna_baf<-getBAFs("RNA.allele.txt","G1000_alleles_hg38_",chrom_names = "chr6")

rna_baf<-as_tibble(rna_baf)%>%
  mutate(end=Position)%>%
  select(Chromosome,Position,end,baf)%>%
  setNames(c("seqnames","start","end","BAF"))
rna_baf$res<-paste(rna_baf$seqnames,rna_baf$start,sep = ":")
rna_baf<-rna_baf%>%
  filter(res%in%keep_alleles)
rna_baf<-rna_baf%>%
  #filter(BAF!=0 & BAF!=1)%>%
  mutate(col=case_when(BAF>=0.75 ~ "top",
                       BAF<=0.25 ~ "bot",
                       TRUE ~ "mid"))
rna_baf_gr<-makeGRangesFromDataFrame(rna_baf,keep.extra.columns = T)
p6<-ggplot()+
  geom_rect(CCND3,aes(xmin=start,xmax=end,ymin=0,ymax=1),
               stat = "identity",
               fill="red",
               alpha=0.1,
               color="#525152")+
  geom_rect(keyRegion,aes(xmin=start,xmax=end,ymin=0,ymax=1),
            stat = "identity",
            fill="#525152",
            alpha=0.1,
            color="#525152")+
  stat_identity(rna_baf_gr,aes(x=start,y=BAF,color=col),
                geom="point",size=0.5)+
  scale_color_manual(values = c("top"="#E93326","bot"="#0100F3","mid"="#3F8C37"))+
  theme_pack_panels()+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(panel.grid.major.y = element_line(linetype = "dashed",color="grey"))
####################################WGS##################################################
####################################CIRCLE##################################################
Circle_cov <-read_tsv("CIRCLE.cov.bdg",col_names=F)
names(Circle_cov)<-c("seqnames","start","end","coverage")
Circle_cov<-Circle_cov%>%
  mutate(coverage=if_else(coverage>300,300,coverage))
Circle_cov_gr<-makeGRangesFromDataFrame(Circle_cov,keep.extra.columns = T)

p3<-ggplot()+
  geom_rect(keyRegion,aes(xmin=start,xmax=end,ymin=0,ymax=300),
            stat = "identity",
            fill="#525152",
            alpha=0.1,
            color="#525152")+
  stat_identity(Circle_cov_gr,aes(x=start,y=coverage),
                fill="#525152",
                geom="area")+
  scale_y_continuous(breaks = c(0,300),limits = c(0,300))+
  theme_pack_panels()

circle_baf<-getBAFs("CIRCLE.allele2.txt","G1000_alleles_hg38_",chrom_names = "chr6")

circle_baf<-as_tibble(circle_baf)%>%
  mutate(end=Position)%>%
  select(Chromosome,Position,end,baf)%>%
  setNames(c("seqnames","start","end","BAF"))
circle_baf$res<-paste(circle_baf$seqnames,circle_baf$start,sep = ":")
circle_baf<-circle_baf%>%
  filter(res%in%keep_alleles)
circle_baf<-circle_baf%>%
  #filter(BAF!=0 & BAF!=1)%>%
  mutate(col=case_when(BAF>=0.75 ~ "top",
                       BAF<=0.25 ~ "bot",
                       TRUE ~ "mid"))
circle_baf_gr<-makeGRangesFromDataFrame(circle_baf,keep.extra.columns = T)
p4<-ggplot()+
  geom_rect(keyRegion,aes(xmin=start,xmax=end,ymin=0,ymax=1),
            stat = "identity",
            fill="#525152",
            alpha=0.1,
            color="#525152")+
  stat_identity(circle_baf_gr,aes(x=start,y=BAF,color=col),
                geom="point",size=0.5)+
  scale_color_manual(values = c("top"="#E93326","bot"="#0100F3","mid"="#3F8C37"))+
  theme_pack_panels()+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(panel.grid.major.y = element_line(linetype = "dashed",color="grey"))
####################################CIRCLE##################################################

cairo_pdf("BAF.trakcs.pdf",width = 7.92,height = 2.91)
tracks(p1,p2,p3,p4,p5,p6,heights = c(0.4,1,0.4,1,0.4,1))
dev.off()

png("BAF.trakcs.png",width = 5.344,height = 3.34,units = "in",res = 300)
tracks(p1,p2,p3,p4,p5,p6,heights = c(0.4,1,0.4,1,0.4,1))
dev.off()
