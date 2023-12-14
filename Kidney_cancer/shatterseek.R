get_shatter<-function(x){
####################################prepare data#############################################
df<-read_tsv(stringr::str_glue("SV/merged/{x}.merged.allsv.filter.txt"))
TRA<-df%>%
  filter(svclass=="TRA")
nonTRA<-df%>%
  filter(svclass!="TRA")

sex_non_TRA<-nonTRA%>%filter(chrom1%in%c("chrX","chrY"))
non_sex_non_TRA<-nonTRA%>%filter(!(chrom1%in%c("chrX","chrY")))

sex_non_TRA<-sex_non_TRA%>%
  mutate(length=abs(start2-start1))%>%
  filter(length>0)%>%
  select(-length)

non_sex_non_TRA<-non_sex_non_TRA%>%
  mutate(length=abs(start2-start1))%>%
  filter(length>1000)%>%
  select(-length)

fin_sv<-rbind(rbind(non_sex_non_TRA,sex_non_TRA),TRA)

fin_sv<-fin_sv%>%
  mutate(svclass2=case_when(svclass=="INS/DUP"~"DUP",
                            (svclass=="INV")&(strand1=="+")&(strand2=="+") ~ "h2hINV",
                            (svclass=="INV")&(strand1=="-")&(strand2=="-") ~ "t2tINV",
                            TRUE ~ svclass))%>%
  select(chrom1,start1,end1,chrom2,start2,end2,sv_id,pe_support,strand1,strand2,svclass2,svmethod)%>%
  dplyr::rename(svclass=svclass2)
rm(df)
rm(non_sex_non_TRA)
rm(nonTRA)
rm(sex_non_TRA)
rm(TRA)
##################################shatterseek#################################
sv = fin_sv# sv_matrix
sv$chrom1=sub("chr","",sv$chrom1)
sv$chrom2=sub("chr","",sv$chrom2)
a=unique(c(which(!sv$chrom1%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")),
           which(!sv$chrom2%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"))))
if(length(a)>0){
  sv2 = sv[-a,]
  sv = sv2
}
cn = readr::read_tsv(stringr::str_glue(stringr::str_glue("CNV/cns/{x}T.cs.rmdup.sort.call.cns")))
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
###################################filter###########################################################
b = chromothripsis@chromSummary
hc_index1 = which((b$clusterSize_including_TRA-b$number_TRA)>6 & b$max_number_oscillating_CN_segments_2_states>7 & b$pval_fragment_joins<0.05 & (b$chr_breakpoint_enrichment<0.05 | b$pval_exp_cluster <0.05))
hc_index2 = which((b$clusterSize_including_TRA-b$number_TRA)>3 & b$number_TRA>=4 & b$max_number_oscillating_CN_segments_2_states>7 & b$pval_fragment_joins<0.05)
hc_index = unique(c(hc_index1,hc_index2))
lc_index = which((b$clusterSize_including_TRA-b$number_TRA)>6 & b$max_number_oscillating_CN_segments_2_states_chr<7 & b$max_number_oscillating_CN_segments_2_states_chr>3 & b$pval_fragment_joins<0.05 & (b$chr_breakpoint_enrichment<0.05 | b$pval_exp_cluster <0.05))
b_hc = b[hc_index,]
b_lc = b[lc_index,]
b_hc$inter_other_chroms_coords_all=gsub('\n',';',b_hc$inter_other_chroms_coords_all)
b_lc$inter_other_chroms_coords_all=gsub('\n',';',b_lc$inter_other_chroms_coords_all)
#plots=plot_chromothripsis(ShatterSeek_output = chromothripsis,chr="chr17",genome = "hg38")
#plot = arrangeGrob(plots[[1]],plots[[2]],plots[[3]],plots[[4]],nrow=4,ncol=1,heights=c(0.2,.4,.4,.2))
#plot(plot)
write.table(file=stringr::str_glue("SV/chromothripsis/chromothripis_hc_{x}.matrix"),b_hc,sep="\t",quote=FALSE,row.names=F)
write.table(file=stringr::str_glue("SV/chromothripsis/chromothripis_lc_{x}.matrix"),b_lc,sep="\t",quote=FALSE,row.names=F)
}

#sv_samples<-c("2","3","4","5","6","7","8",
#              "9","11","12","14","15","16",
#              "17","18","19","20","21","22",
#              "24","25","26","27","28","29",
#              "30","31","32","33","34","35",
#              "36","38","39","40","41","42",
#              "43","44","45","46","47","48",
#              "49","50","51","52","53","54",
#              "55","56","57","58","59","60",
#              "61","62","63","64","65","66",
#              "68","70","71","72","73","74",
#              "75","76","78","80","81","82",
#              "83","84","85")

sv_samples<-c("1","23","79")
lapply(sv_samples,get_shatter)
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




df<-read_tsv(stringr::str_glue("SV/merged/21.merged.allsv.filter.txt"))
TRA<-df%>%
  filter(svclass=="TRA")
nonTRA<-df%>%
  filter(svclass!="TRA")

sex_non_TRA<-nonTRA%>%filter(chrom1%in%c("chrX","chrY"))
non_sex_non_TRA<-nonTRA%>%filter(!(chrom1%in%c("chrX","chrY")))

sex_non_TRA<-sex_non_TRA%>%
  mutate(length=abs(start2-start1))%>%
  filter(length>0)%>%
  select(-length)

non_sex_non_TRA<-non_sex_non_TRA%>%
  mutate(length=abs(start2-start1))%>%
  filter(length>1000)%>%
  select(-length)

fin_sv<-rbind(rbind(non_sex_non_TRA,sex_non_TRA),TRA)

fin_sv<-fin_sv%>%
  mutate(svclass2=case_when(svclass=="INS/DUP"~"DUP",
                            (svclass=="INV")&(strand1=="+")&(strand2=="+") ~ "h2hINV",
                            (svclass=="INV")&(strand1=="-")&(strand2=="-") ~ "t2tINV",
                            TRUE ~ svclass))%>%
  select(chrom1,start1,end1,chrom2,start2,end2,sv_id,pe_support,strand1,strand2,svclass2,svmethod)%>%
  dplyr::rename(svclass=svclass2)
rm(df)
rm(non_sex_non_TRA)
rm(nonTRA)
rm(sex_non_TRA)
rm(TRA)
##################################shatterseek#################################
sv = fin_sv# sv_matrix
sv$chrom1=sub("chr","",sv$chrom1)
sv$chrom2=sub("chr","",sv$chrom2)
a=unique(c(which(!sv$chrom1%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")),
           which(!sv$chrom2%in%c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"))))
if(length(a)>0){
  sv2 = sv[-a,]
  sv = sv2
}
cn = readr::read_tsv(stringr::str_glue(stringr::str_glue("CNV/cns/21T.cs.rmdup.sort.call.cns")))
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
###################################filter###########################################################
b = chromothripsis@chromSummary
for(i in c(1,2,3,4,5,6,7,8,11,12,14,15,17,18,20,21)){
  plots=plot_chromothripsis(ShatterSeek_output = chromothripsis,chr=paste0("chr",i),genome = "hg38")
  plot = arrangeGrob(plots[[1]],plots[[2]],plots[[3]],plots[[4]],nrow=4,ncol=1,heights=c(0.2,.4,.4,.2))
  pdf(paste0("SV/chromothripsis/sample21.chr",i,".pdf"),width = 4.77,height = 4.1)
  plot(plot)
  dev.off()
}
