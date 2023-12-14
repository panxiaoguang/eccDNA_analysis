
#################################
#### CNAs
num_cnas<-function(x){
  read_tsv(stringr::str_glue("CNV/cns_call/{x}T.cs.rmdup.sort.call.cns"))%>%filter(cn!=2)%>%nrow()
}

cna_numbers<-sapply(infos$Sample2,num_cnas)

#### SVs
#num_SV<-function(x){
#  tmp<-read_tsv(stringr::str_glue("SV/merged/{x}.merged.allsv.filter.txt"))
#  stst<-tmp%>%
#    group_by(svclass)%>%
#    summarise(count=n())
#  names(stst)[2]<-x
#  stst
#}

num_SV<-function(x){
  read_tsv(stringr::str_glue("SV/merged/{x}.merged.allsv.filter.txt"))%>%
    nrow()
  
}

#hebing<-\(x,y){
#  full_join(x,y,by="svclass")
#}

SVnumbers<-sapply(infos$Sample2,num_SV)
#SVnumbers<-Reduce(hebing,lapply(infos$Sample2,num_SV))
#SVnumbers<-SVnumbers%>%
#  tibble::column_to_rownames(var="svclass")%>%
#  as.matrix()
#SVnumbers[is.na(SVnumbers)]<-0

#SVnumbers<-t(SVnumbers)
#SVnumbers_buchong<-matrix(rep(0,8),nrow = 2)
#rownames(SVnumbers_buchong)<-c("116","122")
#colnames(SVnumbers_buchong)<-colnames(SVnumbers_buchong)
#SVnumbers<-rbind(SVnumbers,SVnumbers_buchong)
##########################################################
###SNV
mutations<-df2%>%
  filter(Variant_Type%in%c("SNP","INS","DEL"))%>%
  group_by(Tumor_Sample_Barcode)%>%
  summarise(count=n())
mutations<-mutations%>%
  tibble::column_to_rownames(var = "Tumor_Sample_Barcode")%>%
  as.data.frame()
mutations<-mutations[infos$Sample2,]
##################################################3###################
df<-read_tsv("CNV/gistic2_results/broad_values_by_arm.txt")
amp_mat<-df%>%tibble::column_to_rownames(var="Chromosome Arm")%>%as.matrix()
amp_mat[amp_mat<=0.5]<-0
amp_mat[amp_mat>0.5]<-1
amp_mat<-amp_mat[c("22p","14p","13p","15p","7p","7q","5p","5q"),]
amp_names<-stringr::str_remove(colnames(amp_mat),"T")
colnames(amp_mat)<-amp_names
lessNames<-setdiff(as.character(infos$Sample2),amp_names)
amp_buchong<-matrix(rep(0,8*length(lessNames)),nrow = 8)
colnames(amp_buchong)<-lessNames
rownames(amp_buchong)<-rownames(amp_mat)
amp<-cbind(amp_mat,amp_buchong)

del_mat<-df%>%tibble::column_to_rownames(var="Chromosome Arm")%>%as.matrix()
del_mat[del_mat>=(-0.5)]<-0
del_mat[del_mat<(-0.5)]<-1
del_mat<-del_mat[c("3p","15p","14p","13p","14q","8p","9q"),]
del_names<-stringr::str_remove(colnames(del_mat),"T")
colnames(del_mat)<-del_names
lessNames<-setdiff(as.character(infos$Sample2),del_names)
del_buchong<-matrix(rep(0,7*length(lessNames)),nrow = 7)
colnames(del_buchong)<-lessNames
rownames(del_buchong)<-rownames(del_mat)
del<-cbind(del_mat,del_buchong)
del<-as.data.frame(del)
###############################################################
#############            annotation               #############
#############                                     #############
###############################################################
infos<-openxlsx::read.xlsx("Table.S1.xlsx",sheet = 5)
infos<-as_tibble(infos)
## should make age as category 
infos<-infos%>%
  mutate(Age2=case_when((Age>=17)&(Age<=39) ~ "17-39",
                       (Age>=40)&(Age<=59) ~ "40-59",
                       (Age>=60)&(Age<=84) ~ "60-84",
                       TRUE ~ "oth"))
infos$Sample2<-as.character(infos$Sample2)
infos$Grade<-as.character(infos$Grade)
infos$stage<-as.character(infos$stage)


#########################plot test
amp<-amp[,infos$Sample2]
del<-del[,infos$Sample2]
#SVnumbers<-SVnumbers[infos$Sample2,]
samplenumbers<-apply(amp,1,function(x) sum(x))
samplenumbers2<-apply(del,1,function(x) sum(x))
#load("sv-event-numbers.Rdata")
#load("cna.numbers.Rdata")
ht1<-Heatmap(amp,
        show_column_names = F,
        cluster_columns = T,
        cluster_rows = T,
        show_column_dend = F,
        show_row_dend = F,
        show_heatmap_legend = F,
        row_names_side = "left",
        row_names_gp = gpar(fontface="italic",fontsize=8),
        rect_gp = gpar(type="none"),
        #column_split = infos$Histology,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width*0.8, height = height*0.9, 
                    gp = gpar(col = "white", fill = "#F6F8FD"))
          if(amp[i,j]==1){
            grid.rect(x = x, y = y, width = width*0.8, height = height*0.9, 
                      gp = gpar(fill = "#6594CB",col="white"))
          }
        },
        right_annotation  = rowAnnotation(mutations=anno_barplot(samplenumbers,
                                                                 border=F,
                                                                 axis=T,
                                                                 bar_width = 0.45,
                                                                 axis_param = list(
                                                                                   side="top"),
                                                                 gp=gpar(fill="#6594CB",
                                                                         col="transparent")),
                                          show_annotation_name=F),
        top_annotation = HeatmapAnnotation(Histology=infos$Histology,
                                           Gender=infos$Gender,
                                           Age=infos$Age2,
                                           Grade=infos$Grade,
                                           Stage=infos$stage,
                                           Focal=infos$Focal,
                                           EPM=anno_barplot(infos$EPM,
                                                           height=unit(2, "cm"),
                                                           border = F,
                                                           axis = T,
                                                           bar_width = 1,
                                                           gp=gpar(fill="#41ab5d",
                                                                   col="transparent")),
                                           SVs=anno_barplot(svTable$sv2,
                                                             border=F,
                                                             axis=T,
                                                             height = unit(2, "cm"),
                                                             bar_width = 1,
                                                             #gp = gpar(fill = c("#a6bddb","#fdc086","#e78ac3","#a6d854"), col = "transparent")
                                                            gp = gpar(fill="#fdc086",col="transparent")
                                                             ),
                                           CNAs=anno_barplot(cna_numbers,
                                                             border=F,
                                                             axis=T,
                                                             height = unit(2, "cm"),
                                                             bar_width = 1,
                                                             gp=gpar(fill="#a6bddb",
                                                                     col="transparent")),
                                           MUTs=anno_barplot(mutationTable$muta2,
                                                             border=F,
                                                             axis=T,
                                                             height = unit(2, "cm"),
                                                             bar_width = 1,
                                                             gp=gpar(fill="#bebada",
                                                                     col="transparent")),
                                           simple_anno_size = unit(0.3, "cm"),
                                           col = list(Gender=c("Female"="#E93420","Male"="#316DBB"),
                                                      Age=c("17-39"="#D3D4E8","40-59"="#ABABCD","60-84"="#776EB1"),
                                                      Histology=c("ccRCC"="#EBABAC","CDC"="#7662A5","chRCC"="#A68575","pRCC"="#73A99B","sRCC"="#747070"),
                                                      Grade=c("0"="#747070","1"="#F5F6F7","2"="#E1E2D2","3"="#8A8A8A","4"="#EBD4E4"),
                                                      Stage=c("1"="#F5F6F7","2"="#E1E2D2","3"="#8A8A8A","4"="#EBD4E4"),
                                                      Focal=c("0"="#747070","1"="#F5F6F7"))
                                           
        )
)

ht2<-Heatmap(del,
             show_column_names = F,
             cluster_columns = T,
             cluster_rows = T,
             show_column_dend = F,
             show_row_dend = F,
             show_heatmap_legend = F,
             row_names_gp = gpar(fontface="italic",fontsize=8),
             row_names_side = "left",
             rect_gp = gpar(type="none"),
             #column_split = infos$Histology,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width*0.8, height = height*0.9, 
                         gp = gpar(col = "white", fill = "#FFFDF7"))
               if(del[i,j]==1){
                 grid.rect(x = x, y = y, width = width*0.8, height = height*0.9, 
                           gp = gpar(fill = "#D79E62",col="white"))
               }
             },
             right_annotation  = rowAnnotation(mutations=anno_barplot(samplenumbers2,
                                                                      border=F,
                                                                      axis=T,
                                                                      axis_param = list(
                                                                                        side="top"),
                                                                      bar_width = 0.45,
                                                                      gp=gpar(fill="#D79E62",
                                                                              col="transparent")),
                                               show_annotation_name=F)
)




pdf("test_oncoprint.pdf",width = 11.76,height = 4.09)
draw(ht)
dev.off()


###############################finalal##########################################
###parse_Maf
## load filter genes,
sigGenes<-openxlsx::read.xlsx("KCA_sign_gene.xlsx")

## use maftools to load maf file
Variants <- c("Missense_Mutation",
              "Nonsense_Mutation", 
              "Nonstop_Mutation",
              "Frame_Shift_Ins",
              "Frame_Shift_Del",
              "In_Frame_Del",
              "In_Frame_Ins",
              "Splice_Site",
              "Translation_Start_Site")

df2<-read_tsv("KCA_Combine.maf")

df2<-df2%>%
  filter(Variant_Classification%in%Variants)
df2$Tumor_Sample_Barcode<-stringr::str_remove(df2$Tumor_Sample_Barcode,"CCGA-RCC-")
df2$Tumor_Sample_Barcode<-stringr::str_remove(df2$Tumor_Sample_Barcode,"^0+")
df2$Tumor_Sample_Barcode<-stringr::str_remove(df2$Tumor_Sample_Barcode,"T")
#df2<-df2%>%
#  filter(!Tumor_Sample_Barcode%in%c("68","70"))
maf <- read.maf(df2)
#pdf("oncoprint.maf.pdf",width = 12,height = 6.57)
oncoplot(maf,top = 30,writeMatrix = T,genes = unique(sigGenes$Gene_hg38),removeNonMutated=T)
#dev.off()
plotData<-read.table("onco_matrix.txt",sep="\t",check.names = F)
otherSamples<-setdiff(infos$Sample2,colnames(plotData))
plotData_buchong<-matrix(rep("",20*length(otherSamples)),nrow =20 )
colnames(plotData_buchong)<-otherSamples
rownames(plotData_buchong)<-rownames(plotData)
plotData<-cbind(plotData,plotData_buchong)
plotData<-plotData[,infos$Sample2]
###################transition
trans<-df2%>%
  filter(Variant_Type=="SNP")%>%
  select(11:13)%>%
  mutate(type=paste0(Reference_Allele,">",Tumor_Seq_Allele2))%>%
  filter(type%in%c("C>T","C>G","C>A","T>A","T>C","T>G"))%>%
  group_by(Tumor_Sample_Barcode,type)%>%
  summarise(count=n())%>%
  mutate(pct=count/sum(count))%>%
  select(1,2,4)%>%
  tidyr::pivot_wider(names_from = "type",
                     values_from = pct)%>%
  tibble::column_to_rownames(var="Tumor_Sample_Barcode")%>%
  as.matrix()

trans[is.na(trans)]<-0
trans<-trans[infos$Sample2,]
col = c("Missense_Mutation" = "#a6bddb", 
        "Nonsense_Mutation" = "#fdc086", 
        "Frame_Shift_Ins"="#e78ac3",
        "Frame_Shift_Del"="#a6d854",
        "In_Frame_Del"="#96CED5",
        "In_Frame_Ins"="#CF4A31",
        "Nonstop_Mutation"="#484125",
        "Splice_Site"="#beaed4",
        "Translation_Start_Site"="#F88D51",
        "Multi_Hit"="#DC494C")

alter_fun = list(
  background = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = "#f0f0f0"),   
  Missense_Mutation = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Missense_Mutation"]),
  Nonsense_Mutation = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Nonsense_Mutation"]),
  Frame_Shift_Ins = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Frame_Shift_Ins"]),
  Frame_Shift_Del = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Frame_Shift_Del"]),
  In_Frame_Del = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["In_Frame_Del"]),
  In_Frame_Ins = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["In_Frame_Ins"]),
  Nonstop_Mutation = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Nonstop_Mutation"]),
  Splice_Site = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Splice_Site"]),
  Translation_Start_Site = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Translation_Start_Site"]),
  Multi_Hit = alter_graphic("rect", horiz_margin = unit(0.5, "pt"), vertical_margin = unit(0.5, "pt"),fill = col["Multi_Hit"])
)

ht3<-oncoPrint(plotData,
          alter_fun = alter_fun, 
          col = col,
          alter_fun_is_vectorized = FALSE,
          show_pct = F,
          row_names_side = "left",
          row_names_gp = gpar(fontface="italic",fontsize=8),
          top_annotation = NULL,
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(axis = T,axis_param = list(side="top")),
                                           show_annotation_name=F),
          bottom_annotation = HeatmapAnnotation(tbar=anno_barplot(trans,
                                                                  border=F,
                                                              axis=F,
                                                              height = unit(2, "cm"),
                                                              gp = gpar(fill = c("#5181B0","#415389","#C9534F","#64a060","#E1B057","#DB8E4E"), col = "transparent")),
                                                show_legend=T)
          )

ht<-ht1%v%ht2%v%ht3

pdf("test_muttaion_oncoplot3.pdf",width = 8.84,height = 9)
draw(ht,main_heatmap = 3)
dev.off()

save(SVnumbers,file = "sv-event-numbers.Rdata")
save(cna_numbers,file="cna.numbers.Rdata")
