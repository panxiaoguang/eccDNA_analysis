
samples<-c("100","102","103","104","105","106","107",
           "109","110","111","112","113","114","115",
           "116","117","118","119","11","120","121",
           "122","123","124","125","126","127","128",
           "129","12","130","131","132","133","134",
           "136","138","14","15","16","17","18","19",
           "1","20","21","22","23","24","25","26","27",
           "28","29","2","30","31","32","33","34","35",
           "36","38","39","3","40","41","42","43","44",
           "45","46","47","48","49","4","50","51","52",
           "53","54","55","56","57","58","59","5","60",
           "61","62","63","64","65","66","68","6","70",
           "71","72","73","74","75","76","78","79","7",
           "80","81","82","83","84","85","86","87","88",
           "89","8","90","91","92","93","94","95","96",
           "97","98","99","9")
normals<-paste0(samples,"N")
cases<-paste0(samples,"T")

getLength<-function(x){
  df<-read_tsv(stringr::str_glue("filter/cKca_{x}_circle_site.filter.tsv"))
  tibble(sample=x,length=df$length)
}

fin<-do.call("bind_rows",lapply(c(normals,cases),getLength))

fin<-fin%>%
  mutate(group=if_else(sample %in% cases,"Tumour","Normal"))
##catogray
haha2<-fin%>%
  mutate(gps=cut(length,
                 breaks = c(-1,2000,10000,Inf),
                 labels = c("<2K","2k~10K",">10K")))%>%
  group_by(sample,gps)%>%
  summarise(number=n())%>%
  mutate(ratio=number/sum(number)*100)%>%
  ungroup()

haha2<-haha2%>%
  mutate(group=if_else(sample %in% cases,"Tumour","Normal"))

haha<-haha2%>%group_by(gps)%>%rstatix::wilcox_test(ratio~group)
openxlsx::write.xlsx(list(haha,haha2),"bladder_length_category.xlsx")

save(fin,file="length_total.Rdata")
