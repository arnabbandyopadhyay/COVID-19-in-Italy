library(data.table)
library("viridis")
library("readxl")
library(WriteXLS)
library(R.matlab)
library(ggplot2)


thisishere <- getwd()

timing <- seq(as.Date("2020-02-24"), as.Date("2020-05-21"), by="days")
end <- length(timing)


#north <- c("Lombardia","Liguria","Piemonte","Valledaosta","Bolzano","Emilia","Friuli","Veneto","Trento")
#centre <- sort(c("Lazio","Marche","Toscana","Umbria"))
#southisl <- sort(c("Abruzzo","Basilicata","Calabria","Campania","Molise","Puglia","Sardegna","Sicilia"))

###main text regions
#regionSelection <- c("Lombardia","Liguria","Piemonte","Valle d'Aosta","Emilia Romagna","Marche")
#regSel <- c("Lombardia","Liguria","Piemonte","Valledaosta","Emilia","Marche")

###supplementary regions
#suppl <- sort(setdiff(c(north,centre,southisl),regSel))
#supplreg <- c("Abruzzo", "Basilicata", "Bolzano", "Calabria", "Campania", "Friuli" ,"Lazio", "Molise", "Puglia" ,"Sardegna", "Sicilia","Toscana" ,"Trento" ,"Umbria" ,"Veneto")

exgr <- "Lombardia"

###set data frame
dataCollection <- data.frame(matrix(ncol=10,nrow=0))
colnames(dataCollection) <- c("Region","Date","quarantine","hospital","icu","recovered","dead","hlim","icum","Type")

######### results folder
takesimul <- "./"
takedata <- "../Data"
#######

i <- 0

# for (region in c(regSel,suppl)) {
for (region in exgr) {
    i <- i+1
    
    print(region)
    
    ### create data frame for Data
    setwd(takedata)
    dataReal <- read_excel(paste0(region,".xls"))
    setwd(thisishere)
    
    daynum <- nrow(dataReal)
    missing <- length(timing) - daynum
    
    #######
    timapp <- (timing[(missing+1):end])
    regapp <- rep(region,daynum)
    typapp <- rep("data",daynum)
    ########
    
    bindData <- as.data.frame(cbind(date=timapp,dataReal[which(colnames(dataReal)=="qua")], dataReal[which(colnames(dataReal)=="hos")], dataReal[which(colnames(dataReal)=="icu")], dataReal[which(colnames(dataReal)=="rec")], dataReal[which(colnames(dataReal)=="dead")], type=typapp))
    
    setwd(thisishere)
    
    ### collect simulations' results

    setwd(takesimul)
    
    simul_cap <- read.table(paste0("area_spec_",region,".csv"), sep=",")
    
    theEnd <- which(simul_cap[,7]==max(simul_cap[,7]))
    simul_cap <- simul_cap[1:theEnd,]
    total <- (ncol(simul_cap)/6)-1
    num_sim <- seq(6,total+1,6)
    
    
    ###############get parameters
    # a b del hlim icum
    
    dat<-readMat(paste0("parfit_",region,".mat"))
    hlim <- unlist(data.frame(dd=t(dat$savpars[4,])),use.names=F)
    icum <- unlist(data.frame(dd=t(dat$savpars[5,])),use.names=F)
    
    ###remove first window's elements, one corresponding to initial conditions, 3 corresponding to the shifting = 5
    
    hlim <- c(hlim[6:length(hlim)],rep(max(hlim),100))
    icum <- c(icum[6:length(icum)],rep(max(icum),100))
    
    ###################start appending data 
    
    dataCollection <- rbind(dataCollection,cbind(Region=regapp,Date=bindData$date, quarantine=bindData$qua, hospital=bindData$hos, icu=bindData$icu, recovered=bindData$rec, dead=bindData$dead, hlim=hlim[1:daynum],icum=icum[1:daynum], Type=bindData$type))
    
    ####
    setwd(thisishere)
    
    data_tobind <- NULL
    
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
    timapp <- (timing[(missing+1):end])
    typapp <- rep("simul",daynum)
    ####
    to_append <- as.data.frame(cbind(Region=regapp,Date=timapp, quarantine=data_tobind[,1], hospital=data_tobind[,2], icu=data_tobind[,3], recovered=data_tobind[,4], dead=data_tobind[,5], hlim=hlim[1:daynum],icum=icum[1:daynum],Type=typapp))
    
    dataCollection <- rbind(dataCollection,to_append)
}

levels(dataCollection$Date) <- as.Date(timing)
levels(dataCollection$Type) <- c("Data","Simul")

###adjust the regions' name
#levels(dataCollection$Region) <- c(regionSelection,supplreg)
###

for(j in 3:9) {dataCollection[,j]<-as.numeric(as.character(dataCollection[,j]))}

sss <- dataCollection[which(dataCollection$Type=="Simul"),]
ddd <- dataCollection[which(dataCollection$Type=="Data"),]

ls <- 2.5
ps <- 2

fitplot <- ggplot(sss, aes(x=as.Date(Date), y=dead, color="Dead")) +
geom_line(size=ls,alpha=0.5) +
geom_point(ddd, mapping=aes(x=as.Date(Date), y=dead, fill="Dead"),size=ps,col="black",shape=21) +
geom_line(aes(x=as.Date(Date), y=hospital, color="Hosp."),size=ls,alpha=0.5) +
geom_line(aes(x=as.Date(Date), y=quarantine, color="Infected"),size=ls,alpha=0.5) +
geom_line(aes(x=as.Date(Date), y=icu, color="ICU"),size=ls,alpha=0.5) +
geom_point(ddd, mapping=aes(x=as.Date(Date), y=icu, fill="ICU"),size=ps,col="darkred",shape=21) +
geom_point(ddd, mapping=aes(x=as.Date(Date), y=hospital, fill="Hosp."),size=ps,col="darkblue",shape=21) +
geom_point(ddd, mapping=aes(x=as.Date(Date), y=quarantine, fill="Infected"),size=ps,col="darkgreen",shape=21) +

geom_point(aes(x=as.Date(Date), y=hlim, fill="Hosp. Capacity"),size=2,col='darkblue',shape=23) +
geom_point(aes(x=as.Date(Date), y=icum, fill="ICU Capacity" ),size=2,col='darkred',shape=23) +
scale_y_sqrt(breaks = function(x) { sort(c(50,1000,pretty(  seq(min(x),max(x),150)  ))) } ) +
scale_colour_manual(name='', values=c('Infected'='mediumseagreen', 'Hosp.'='dodgerblue3', 'ICU'='firebrick2','Dead'='mistyrose4'), guide='legend') +
scale_fill_manual(name='', values=c('Hosp. Capacity'='dodgerblue3', 'ICU Capacity'='firebrick2', 'Infected'='darkgreen', 'Hosp.'='darkblue', 'ICU'='darkred', 'Dead'='black'),labels=c('Hosp. Capacity'='Hosp. Cap.', 'ICU Capacity'='ICU Cap.', 'Infected'='', 'Hosp.'='', 'ICU'='', 'Dead'=''), guide='legend') +
guides(fill=guide_legend(override.aes=list(fill=c(NA,NA,'dodgerblue3',NA,'firebrick3',NA), shape=c(NA,NA,23,NA,23,NA), size=4, col= c(NA,NA,'darkblue',NA,'darkred',NA) )),colour = guide_legend(override.aes = list(fill=NA, linetype=rep(1,4), shape=rep(19,4) ,size=2) )) +
facet_wrap(~Region,scale="free_y",ncol=3) +
scale_x_date(date_breaks = "14 days",date_labels = "%b %e") +
theme_bw() +
theme(panel.grid.minor = element_blank(),strip.background = element_rect(color="white", fill="white"), strip.text=element_text(size=30*1.3), axis.text.x=element_text(size=25*1.3,angle=45,margin = margin(t = 5),hjust=1), axis.text.y=element_text(size=25), axis.title.y=element_text(size=30*1.3), panel.background = element_rect(fill="white"), legend.text=element_text(size=30*1.3)) +
theme(panel.border = element_rect(colour = "black", size=1.5)) +
xlab("") + ylab("Numbers")


pdf("capfit.pdf",width=16*1.5,height=8*1.5)
# pdf("supplcap.pdf",width=16*1.7,height=15*1.7)
fitplot
dev.off()



