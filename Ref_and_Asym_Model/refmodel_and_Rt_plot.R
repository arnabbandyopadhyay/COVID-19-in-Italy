###Fig. Reference model and time-dependent reproduction number Rt

library(data.table)
library(viridis)
library(readxl)
library(WriteXLS)
library(R.matlab)
library(ggplot2)
library(RColorBrewer)
library(knitr)
library(GGally)
library(ggpubr)
library(ggrepel)
Sys.setlocale("LC_TIME", "C")

###current folder
thisishere <- getwd()

###up to last day of fitting
timing <- seq(as.Date("2020-02-24"), as.Date("2020-07-23"), by="days")
end <- length(timing)

#nation <- "Italy"
#north <- sort(c("Lombardia","Liguria","Piemonte","Valledaosta","Bolzano","Emilia","Friuli","Veneto","Trento"))
#centre <- sort(c("Lazio","Marche","Toscana","Umbria"))
#southisl <- sort(c("Abruzzo","Basilicata","Calabria","Campania","Molise","Puglia","Sardegna","Sicilia"))

### in manin text
#regionSelection <- c("Italy","Emilia Romagna","Lombardia","Liguria","Piemonte","Valle d'Aosta","Veneto","Marche","Toscana")
### in supplementary
#supplreg <- sort(c(southisl,"Lazio","Umbria","Bolzano","Friuli","Trento"))

dataCollection <- data.frame(matrix(ncol=8,nrow=0))
colnames(dataCollection) <- c("Region","Date","quarantine","hospital","icu","recovered","dead","Type")

exgr <- "Lombardia"

#########
#simulation result folder
takesimul <- "und_results"
#data folder
takedata <- "Data"
#######

################################## fit
i <- 0
#for (region in c(nation,north,centre,southisl)) {
for (region in exgr) {
    i <- i+1
    
    print(region)
    
    setwd(takedata)
    dataReal <- read_excel(paste0(region,".xls"))
    
    daynum <- nrow(dataReal)
    missing <- length(timing) - daynum
    
    ####### collection of data
    timapp <- as.Date(timing[(missing+1):end])
    regapp <- rep(region,daynum)
    typapp <- rep("data",daynum)
    ########
    
    bindData <- as.data.frame(cbind(date=timapp,dataReal[which(colnames(dataReal)=="qua")], dataReal[which(colnames(dataReal)=="hos")], dataReal[which(colnames(dataReal)=="icu")], dataReal[which(colnames(dataReal)=="rec")], dataReal[which(colnames(dataReal)=="dead")], type=typapp))
    
    
    dataCollection <- rbind(dataCollection,cbind(Region=regapp,Date=bindData$date, quarantine=bindData$qua, hospital=bindData$hos, icu=bindData$icu, recovered=bindData$rec, dead=bindData$dead, Type=bindData$type))
    
    #colnames(dataCollection) <- c("Region","Date","quarantine","hospital","icu","recovered","dead","Type")

    setwd(thisishere)

    setwd(takesimul)
    
    simul_cap <- read.table(paste0("area_spec_",region,".csv"), sep=",")
    
    ###count the number of windows simulated; ignore the first window
    theEnd <- which(simul_cap[,7]==max(simul_cap[,7]))
    simul_cap <- simul_cap[1:theEnd,]
    total <- (ncol(simul_cap)/6)-1
    num_sim <- seq(6,total+1,6)
    
    setwd(thisishere)
    
    data_tobind <- NULL
    
    ###results dataset is arranged as [time; quar.; hosp.; icu; rec.; dead] for each window
    for (pop in 2:6) {
        tmppop <- NULL
        for (run in 1:total) {
            step <- run*6
            timeone <- which(simul_cap[,1+step]==0)
            if (run<total) {
                tmppop <- c(tmppop,simul_cap[timeone,pop+step])
            } else {
                ### for the last window take all the 7 days
                takesev <- which(simul_cap[,1+step] %in% seq(0,7))
                tmppop <- c(tmppop,simul_cap[takesev,pop+step])
            }
        }
        data_tobind <- cbind(data_tobind,tmppop)
    }
    #####
    timapp <- as.Date(timing[(missing+1):end])
    typapp <- rep("simul",daynum)
    ####
    to_append <- as.data.frame(cbind(Region=regapp,Date=timapp, quarantine=data_tobind[,1], hospital=data_tobind[,2], icu=data_tobind[,3], recovered=data_tobind[,4], dead=data_tobind[,5], Type=typapp))
    #levels(to_append$Date) <- as.Date(timing[(missing+1):end])
    dataCollection <- rbind(dataCollection,to_append)
}
levels(dataCollection$Date) <- as.Date(timing)
levels(dataCollection$Type) <- c("Data","Simul")

### adjust the regions' names
#levels(dataCollection$Region) <- c("Italy","Bolzano","Emilia Romagna","Friuli Venezia Giulia", "Liguria","Lombardia","Piemonte", "Trento", "Valle d'Aosta","Veneto",  "Lazio", "Marche", "Toscana", "Umbria","Abruzzo", "Basilicata", "Calabria", "Campania", "Molise", "Puglia", "Sardegna", "Sicilia")

####suppl (run to generate suppl. figure)
#takereg <- which(dataCollection$Region %in% supplreg)
#dataCollection <- dataCollection[takereg,]

####main (run to generate main text figure)
#takereg <- which(dataCollection$Region %in% regionSelection)
#dataCollection <- dataCollection[takereg,]



for(j in 3:7) {dataCollection[,j]<-as.numeric(as.character(dataCollection[,j]))}

sss <- dataCollection[which(dataCollection$Type=="Simul"),]
ddd <- dataCollection[which(dataCollection$Type=="Data"),]
ls <- 1.3
ps <- 1.1

fitplot <- ggplot(sss, aes(x=as.Date(Date), y=dead, color="Dead")) +
geom_line(size=ls,alpha=0.5) +
geom_line(aes(x=as.Date(Date), y=quarantine, color="Infected"),size=ls,alpha=0.5) +
geom_line(aes(x=as.Date(Date), y=hospital, color="Hosp."),size=ls,alpha=0.5) +
geom_line(aes(x=as.Date(Date), y=icu, color="ICU"),size=ls,alpha=0.5) +
geom_point(ddd, mapping=aes(x=as.Date(Date), y=dead, fill="Dead"),size=ps,col='black',shape=21) +
geom_point(ddd, mapping=aes(x=as.Date(Date), y=hospital, fill="Hosp."),size=ps,col='darkblue',shape=21) +
geom_point(ddd, mapping=aes(x=as.Date(Date), y=quarantine, fill="Infected"),size=ps,col='darkgreen',shape=21) +
geom_point(ddd, mapping=aes(x=as.Date(Date), y=icu, fill="ICU"),size=ps,col='darkred',shape=21) +
scale_y_sqrt(breaks = function(x) { sort(c(50,1000,pretty(  seq(min(x),max(x),100)  ))) } ) +
scale_colour_manual(name='', values=c('Infected'='mediumseagreen', 'Hosp.'='dodgerblue3', 'ICU'='firebrick2','Dead'='mistyrose4'), guide='legend') +
scale_fill_manual(name='', values=c('Infected'='darkgreen', 'Hosp.'='darkblue', 'ICU'='darkred','Dead'='black'), guide='legend') +
guides(fill=guide_legend(override.aes=list(fill=NA)), colour = guide_legend(override.aes = list( linetype=rep(1,4), shape=rep(NA,4), size=ls) )) +
facet_wrap(~Region,scale="free_y",ncol=3) +
scale_x_date(date_breaks = "14 days",date_labels = "%b %e") +
theme_bw() +
theme(panel.grid.minor = element_blank(),strip.background = element_rect(color="white", fill="white"), strip.text=element_text(size=30), axis.text.x=element_text(size=20,angle=45,margin = margin(t = 5),hjust=1), axis.text.y=element_text(size=20), axis.title.y=element_text(size=30), panel.background = element_rect(fill="white"), legend.text=element_text(size=30)) +
theme(panel.border = element_rect(colour = "black", size=1.5)) +
xlab("") + ylab("Numbers")

pdf('fit_Ref.pdf',w=20,h=20*0.5)
fitplot
dev.off()
#############

###########

####################### Rt comparison (Reference vs Asymptomatic model)
#folder for Rt results (Reference)
foldU <- 'r0_und'
#folder for Rt results (Asymptomatic)
foldA <- 'r0_asy'

#current folder
thisishere <- getwd()

###up to last day of fitting
timing <- seq(as.Date("2020-02-24"), as.Date("2020-07-23"), by="days")
end <- length(timing)

###
setwd(foldU)
datixls<-list.files(pattern='.csv')
reg_names <- sapply(strsplit(datixls,"[_ .]+"),function(x) x[3])
setwd(thisishere)

### regions in the main text
regions_sel <- c('Emilia','Italy','Lombardia','Liguria','Valledaosta','Marche','Piemonte','Veneto','Toscana')
### regions in the suppl.
regsuppl <- setdiff(reg_names,regions_sel)

exgr <- "Lombardia"

setwd(thisishere)
rtotU <- NULL
rtotA <- NULL
rtot<-NULL
rtline <- NULL

rr <- NULL
tt <- NULL
dayt <- NULL

###arrange the results
for (reg in exgr) {
    print(reg)
    
    ###Reference model
    setwd(foldU)
    U <- read.table(paste0('final_r0_',reg,'.csv'),sep=',')
    
    howmany <- length(timing)-7-ncol(U)
    
    if (howmany>0) {
        print(howmany)
        U <- cbind(matrix(NA,nrow=nrow(U),ncol=howmany),U)
    }
    
    for (kk in 1:ncol(U)) {
        rtot <- c(rtot,U[,kk])
        dayt <- c(dayt,rep(timing[kk],nrow(U)))
    }
    rr <- c(rr,rep(reg,nrow(U)*ncol(U)))
    tt <- c(tt,rep('Undetected',nrow(U)*ncol(U)))
    setwd(thisishere)

    ###Asymptomatic model
    setwd(foldA)
    A <- read.table(paste0('final_r0_',reg,'.csv'),sep=',')
    if (howmany>0) {
        print(howmany)
        A <- cbind(matrix(NA,nrow=nrow(A),ncol=howmany),A)
    }
    for (kk in 1:ncol(A)) {
        rtot <- c(rtot,A[,kk])
        dayt <- c(dayt,rep(timing[kk],nrow(A)))
    }
    rr <- c(rr,rep(reg,nrow(A)*ncol(A)))
    tt <- c(tt,rep('Asymptomatic',nrow(A)*ncol(A)))
    
    setwd(thisishere)
}

rtot_comp <- cbind(rtot,rr,tt,dayt)
rttot <- data.frame(rtot_comp)

### adjust regions' names (to run for main text)
#levels(rttot$rr) <- c("Emilia Romagna", "Italy", "Liguria", "Lombardia", "Marche", "Piemonte", "Toscana", "Valle d'Aosta", "Veneto")

#rttot$rr <- relevel((rttot$rr), "Italy")

###### adjust regions' names (to run for supplementary)
#levels(rttot$rr) <- c("Abruzzo", "Basilicata", "Bolzano" ,"Calabria" ,"Campania" ,"Friuli Venezia Giulia" ,"Lazio" ,"Molise" ,"Puglia" ,"Sardegna" ,"Sicilia" ,"Trento" ,"Umbria")

##############

levels(rttot$dayt) <- as.Date(timing[1:144])
rttot$rtot <- as.numeric(as.character(rttot$rtot))

###plot

Rtplot <- ggplot(rttot,aes(x=as.Date(dayt), y=rtot, colour=tt, fill=dayt)) +
geom_boxplot(outlier.shape=NA,size=0.1) +
geom_line(aes(x=as.Date("2020-03-10"), y=rtot),colour='dark green',linetype=1,size=0.7) +
geom_line(aes(x=as.Date("2020-05-18"), y=rtot),colour='magenta',linetype=1,size=0.7) +
geom_hline(aes(yintercept=1), linetype=2,colour='black',size=1) +
scale_y_sqrt(breaks = c(.01, .5, 1, 3, 5, 7, 9)) +
facet_wrap(~rr, scale="free_y", ncol=3) +
scale_colour_manual(name='Model:',values=c('Undetected'='orangered', 'Asymptomatic'='royalblue'), labels=c('Undetected'='Ref.', 'Asymptomatic'='Asympt.'),  guide='legend') +
scale_x_date(date_breaks = "14 days",date_labels = "%b %e") +
scale_fill_manual(name='',values=rep('white',length(dayt)),guide=FALSE) +
theme_bw() +
theme(panel.grid.minor = element_blank(),strip.background = element_rect(color="white", fill="white"), strip.text=element_text(size=30), axis.text.x=element_text(size=20,angle=45,margin = margin(t = 5),hjust=1), axis.text.y=element_text(size=20), axis.title.y=element_text(size=30), panel.background = element_rect(fill="white"), legend.text=element_text(size=30),legend.title=element_text(size=30)) +
theme(panel.border = element_rect(colour = "black", size=1.5)) +
xlab("") + ylab("Rt")

pdf('onlyRt.pdf',w=20,h=20*0.5)
Rtplot
dev.off()
#############

########arrange fit and Rt plot
two <- ggarrange(fitplot,Rtplot,
    labels=c("A","B"),nrow=2,font.label = list(size = 25))
    
pdf("RefFit_Rt.pdf",w=20*1.3,h=20*1.3)
two
dev.off()



