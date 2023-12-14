#load("length_total.Rdata")
library(ggplot2)


plotData2<-fin%>%
  filter(length<=2000)

ggplot(plotData2,aes(x=length,y=after_stat(count)))+
  geom_density(aes(color=group))+
  scale_color_manual(values = c("Normal"="#415389","Tumour"="#C9534F"))+
  #scale_x_log10(breaks = c(100,500,500,1000,10000,100000),
  #             labels = c("100","300","500","1K","10K","100K"))+
  theme_pubr()

ggsave("length.plot.pdf",width =6.47 ,height =2.98 )  


testData<-plotData2%>%
  filter(sample=="100N")

### a function plot length with peak numbers
plot_length_dist<-function(iData){
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  p<-ggplot(iData,aes(x=length))+
    geom_density()
  layerData<-layer_data(p,1)%>%
    as_tibble()
  dt<-layerData$density
  lg <-layerData$x
  guess_x <-c()
  for (i in 1:length(lg)){
    if((i-5>1)&(i+5)<length(lg)){
      if((dt[i]-dt[i-5]>0)&(dt[i]-dt[i+5]>0)){
        guess_x<-c(guess_x,i)
      }
    }
  }
  md <-layerData[guess_x,]%>%
    mutate(gp=round(x/100))%>%
    group_by(gp)%>%
    summarise(x=median(x))
  labelData<-layerData[guess_x,]%>%
    filter(x%in%(md$x))
  p+geom_text_repel(aes(x=x,y=y,label=round(x)),
                    data = labelData,
                    color="red",
                    min.segment.length = 0.00000001,
                    arrow = arrow(length = unit(0.02, "npc")),
                    box.padding = 1)+
    theme_classic(base_size = 14,base_family = "arial")+
    xlab("Length distribution")+
    ylab("Density")
}

iData<-readRDS("iData_test.rds")
plot_length_dist(iData)

