library(ggplot2)
library(dplyr)

thisishere <- getwd()

setwd(thisishere)

###For Fig 5b
#regions<-c('Lombardia','Liguria','Piemonte','Valledaosta','Emilia','Marche')
#Regions<-c('Lombardia','Liguria','Piemonte',"Valle d'Aosta",'Emilia Romagna','Marche')
#missing=c(0,1,0,10,0,2);

###e.g.
regions<-Regions<-"Lombardia"
missing<-0

data<-NULL
for (i in 1:length(regions)){
  print(i)
  dat<-read.csv(paste0('excess_dead_',regions[i],'.csv'))
  dat$region<-rep(Regions[i],dim(dat)[1])
  dat<-dat[!dat$X0%%1,]
  dat<-distinct(dat, X0, .keep_all = TRUE)
  if (missing[i]>0){
    dat$date<-seq(as.Date("2020-02-27")+missing[i], as.Date("2020-02-27")+80, by="days")
  }
  else{
    dat$date<-seq(as.Date("2020-02-27")+missing[i], as.Date("2020-02-27")+80, by="days")
  }
  data<-rbind(data,dat)
}

data<-data.frame(data)
colnames(data)<-c('t','mean','sd','Regions','Date')

data[which(data$mean < 0.0001),2]<-0
data[which(data$sd < 0.0001),3]<-0

### to take as resulting from extra_dead.m
# savper<-c('10.67%','2.96%','0.2%','9.20%','3.74%','2.85%')
# pos<-sapply(1:6,function(x) {data[tail(which(data$Regions==Regions[x]),1),2]})

###e.g.
savper<-c('10.34%')
pos<-sapply(1,function(x) {data[tail(which(data$Regions==Regions[x]),1),2]})

textmat<-data.frame(Regions=Regions,percent=savper,pos=pos)

Sys.setlocale("LC_ALL", "en_GB.UTF-8")

pp<-ggplot(data,aes(x=Date,y=mean,col=Regions))+
  geom_ribbon(aes(ymin=as.vector(mean-sd),ymax=as.vector(mean+sd),fill=Regions),alpha=0.3,colour = "white")+
  geom_line(aes(y=mean))+geom_text(textmat, mapping=aes(x=as.Date("2020-05-20"), y=pos, label = percent), show.legend=F,size=4)+
  scale_y_sqrt(breaks = function(x) { sort(c(seq(0,50,10),seq(1000,4000,500)) )  } )+
  scale_x_date(date_breaks = "7 days",date_labels = "%b %d")+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),strip.background = element_rect(color="white", fill="white"), 
        strip.text=element_text(size=14), axis.text.x=element_text(size=18,angle=45,
                                                                   margin = margin(t = 5),hjust=1), legend.title=element_text(size=16),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=18), panel.background = element_rect(fill="white"), legend.text=element_text(size=20)) +
  xlab("") + ylab("Saved")
  
pdf('saved_graph.pdf',h=7,w=15)
pp
dev.off()
