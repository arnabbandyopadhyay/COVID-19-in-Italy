### Fig. Impact of testing and islation on hospitalization and death (panel A and B)

library(data.table)
library(viridis)
library(readxl)
library(WriteXLS)
library(R.matlab)
library(ggplot2)
library(directlabels)
library(matrixStats)
library(ggrepel)
library(ggpubr)

thisishere <- getwd()

timing <- seq(as.Date("2020-02-24"), as.Date("2020-07-23"), by="days")
end <- length(timing)-7

north <- c("Lombardia","Liguria","Piemonte","Valledaosta","Emilia","Friuli","Veneto")
centre <- c("Lazio","Marche","Toscana","Umbria")
southisl <- c("Abruzzo","Basilicata","Calabria","Campania","Molise","Puglia","Sardegna","Sicilia")

#######
data_toplot <- data.frame(matrix(nrow=0,ncol=7))
colnames(data_toplot) <- c("Regions","Date","inf_det","hospital","dead_det","total","Model")
#######

###results' folders
##

ref_fold <- "./und_results" # results from Reference model
s22 <- "res22" #undetected decrease from 15.03
s15 <- "res15" #undetected decrease from 10.03
s8 <- "res8" #undetected decrease from 02.03


foldcoll <- c(ref_fold,s8,s15,s22)
labelll <- c('Flat Undetected Value','Decrease from 2nd of March', 'Decrease from 10th of March','Decrease from 15th of March')

labelll <- c('Constant Und.','Und. drop: 02.03', 'Und. drop: 10.03','Und. drop: 15.03')

i <- 0
f <- 0

regdone <- c("Italy","Lombardia","Veneto")

Sys.setlocale("LC_TIME", "C")

###collect results and create dataset
for (folder in foldcoll) {
    setwd(thisishere)
    
    f <- f+1
    
    setwd(folder)
    print(folder)
    for (region in regdone) {
        print(region)
        i <- i+1
        if (folder==ref_fold) {
            simul_cap <- read.table(paste0("area_spec_",region,".csv"), sep=",")
        } else {simul_cap <- read.table(paste0("area_spec_",region,"_fixed.csv"), sep=",")}
        theEnd <- which(simul_cap[,7]==max(simul_cap[,7]))
        simul_cap <- simul_cap[1:theEnd,]
        total <- (ncol(simul_cap)/6)-1
        num_sim <- seq(6,total+1,6)
       
        missing <- length(timing)-total-7
        dead_cap <- NULL
        survive <- NULL
        
        ####################### undetected
        
        total_infec <- NULL
        
        hospital <- NULL
        icu <- NULL
        dead_det <- NULL
        rowelem <- NULL
        
        for (pop in c(2,3,4,6)) {
            detected <- NULL
            count <- 0
            poptwopoints<-NULL
            for (run in 1:total) {
                step <- run*6
                timeone <- which(simul_cap[,1+step]==0)
                detected <- c(detected,simul_cap[timeone,pop+step])
            }
            rowelem <- cbind(rowelem,detected)
        }
        rowelem <- cbind(rowelem[,1],rowelem[,2]+rowelem[,3],rowelem[,4],rowelem[,1]+rowelem[,2]+rowelem[,3])
        colnames(rowelem) <- c("inf_det","hospital","dead_det","total")
        reg <- rep(region,nrow(rowelem))
        datapp <- as.Date(timing[(missing+1):end])
        typapp <- rep(folder,nrow(rowelem))
        #typapp <- rep(labelll[f],nrow(rowelem))
        data_toplot <- rbind(data_toplot,cbind(Regions=reg,Date=datapp,rowelem,Model=typapp))
        
    }
    setwd(thisishere)
}
##############

partial <- timing[1:100]
data_toplot <- data_toplot[which(data_toplot$Date %in% partial), ]

##########################quantify
uu <- data_toplot[which(data_toplot$Model=="./und_results"),]
aa <- data_toplot[which(data_toplot$Model!="./und_results"),]
#################
for (jj in 3:(ncol(data_toplot)-1)) {
    data_toplot[,jj] <- as.numeric(as.character(data_toplot[,jj]))
    uu[,jj] <- as.numeric(as.character(uu[,jj]))
    aa[,jj] <- as.numeric(as.character(aa[,jj]))
}

#################### simulations and data comparison
quant_hos<-data.frame(matrix(nrow=0,ncol=4))
quant_die<-data.frame(matrix(nrow=0,ncol=4))
colnames(quant_hos) <- c("Regions","perc","where","Model")
colnames(quant_die) <- c("Regions","perc","where","Model")
j<-0
for (rr in levels(data_toplot$Regions)) {
    print(rr)
    tmpu1 <- (uu[which(uu$Regions==rr),])
    tmpa1 <- (aa[which(aa$Regions==rr),])
    tmpu_d1 <- (uu[which(uu$Regions==rr),])
    tmpa_d1 <- (aa[which(aa$Regions==rr),])
    
    for(ff in foldcoll) {
        
        if (ff=="./und_results") {
            tmpu <- as.numeric(as.character(tmpu1[which(tmpu1$Model==ff),4]))
            tmpu_d <- as.numeric(as.character(tmpu_d1[which(tmpu_d1$Model==ff),5]))
        } else {
            j<-j+1
            tmpa <- as.numeric(as.character(tmpa1[which(tmpa1$Model==ff),4]))
            tmpa_d <- as.numeric(as.character(tmpa_d1[which(tmpa_d1$Model==ff),5]))
        
            quant_hos[j,] <- cbind(Regions=rr, perc=(max(tmpu)-max(tmpa))/max(tmpu), where=tmpa[length(tmpa)], Model=ff )
            quant_die[j,] <- cbind(Regions=rr, perc=(max(tmpu_d)-max(tmpa_d))/max(tmpu_d), where=tmpa_d[length(tmpa_d)], Model=ff)
        }
    }
    
    print(max(tmpu))
    print(max(tmpa))
     #     print(quant_hos)
}
####################

quant_hos$perc <- round(as.numeric(as.character(quant_hos$perc)),digits=2)
quant_hos$where <- as.numeric(as.character(quant_hos$where))
quant_die$perc <- round(as.numeric(as.character(quant_die$perc)),digits=2)
quant_die$where <- as.numeric(as.character(quant_die$where))

quant_hos$Model <- factor(quant_hos$Model)
quant_die$Model <- factor(quant_die$Model)

levels(data_toplot$Model) <- labelll
levels(quant_hos$Model) <- c(labelll[3],labelll[4],labelll[2])
levels(quant_die$Model) <- c(labelll[3],labelll[4],labelll[2])

save <- data_toplot
for (jj in 3:(ncol(data_toplot)-1)) {
    data_toplot[,jj] <- as.numeric(as.character(data_toplot[,jj]))
    uu[,jj] <- as.numeric(as.character(uu[,jj]))
    aa[,jj] <- as.numeric(as.character(aa[,jj]))
}

levels(data_toplot$Date) <- timing
levels(aa$Date) <- timing
levels(uu$Date) <- timing



#########

ls=1.5
ts=7

lx <- which(timing=="2020-05-28")+13


alter_x <- -1.5
alter_y <- c(-500,+100,+1000,-300,0,+500,-30,0,30)
alter_yd <- c(-1000,+700,+3000,-1000,0,+1500,-100,0,+150)

#####hosp
timeie<-0.7
comphos <- ggplot(data_toplot, aes(x=as.Date(Date), y=hospital, colour=Model)) +
geom_line(size=ls) +
geom_text(quant_hos, mapping=aes(x=timing[lx]+alter_x, y=where+alter_y, label = paste0(100*perc," %"),colour=Model), show.legend=F,size=ts) +
scale_y_sqrt(breaks = function(x) { sort(c(1500,5000,10000,pretty(seq(min(x),max(x),100))  )) } ) +
facet_wrap(~Regions, scale="free_y") +
scale_x_date(date_breaks = "14 days",date_labels = "%b %e",limits=c(timing[1],timing[lx])) +
theme_bw() +
theme(panel.grid.minor = element_blank(),strip.background = element_rect(color="white", fill="white"), strip.text=element_text(size=30*timeie), axis.text.x=element_text(size=20*timeie,angle=45,margin = margin(t = 5),hjust=1), axis.text.y=element_text(size=20*timeie), axis.title.y=element_text(size=30*timeie), panel.background = element_rect(fill="white"), legend.text=element_text(size=30*timeie),legend.title=element_text(size=30*timeie)) +
theme(panel.border = element_rect(colour = "black", size=1.5*timeie)) +
labs(colour="Testing Model") +
xlab("") + ylab("Hospitalized")


########dead
compd <- ggplot(data_toplot, aes(x=as.Date(Date), y=dead_det, colour=Model)) +
geom_line(size=ls) +
geom_text(quant_die, mapping=aes(x=timing[lx]+alter_x, y=where+alter_yd, label = paste0(100*perc," %"),colour=Model), show.legend=F,size=ts) +
scale_y_sqrt(breaks = function(x) { sort(c(1500,5000,10000,pretty(seq(min(x),max(x),100))  )) } ) +
facet_wrap(~Regions, scale="free_y") +
scale_x_date(date_breaks = "14 days",date_labels = "%b %e",limits=c(timing[1],timing[lx])) +
theme_bw() +
theme(panel.grid.minor = element_blank(),strip.background = element_rect(color="white", fill="white"), strip.text=element_text(size=30*timeie), axis.text.x=element_text(size=20*timeie,angle=45,margin = margin(t = 5),hjust=1), axis.text.y=element_text(size=20*timeie), axis.title.y=element_text(size=30*timeie), panel.background = element_rect(fill="white"), legend.text=element_text(size=30*timeie),legend.title=element_text(size=30*timeie)) +
theme(panel.border = element_rect(colour = "black", size=1.5*timeie)) +
labs(colour="Testing Model") +
xlab("") + ylab("Dead")


k<-ggarrange(comphos, compd, labels=c("A","B"), font.label=list(size=30), nrow=2)

pdf('test_model.pdf',h=10*1.5,w=16*1.5)
k
dev.off()
