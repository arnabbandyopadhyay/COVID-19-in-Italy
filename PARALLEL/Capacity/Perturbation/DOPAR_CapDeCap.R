library(foreach)
library(doParallel)
library(parallel)


library(data.table)
library(readxl)
library(WriteXLS)
library(R.matlab)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(tikzDevice)
library(latex2exp)


# Calculate the number of cores
no_cores <- detectCores() - 1
 
# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

thisishere <- getwd()


setwd('cap_res')
fordd <- list.files(pattern = "\\.mat$")
regions <- list()
for (ii in 1:length(fordd)) {
    ff <- fordd[ii]
    regions[[ii]] <- c(unlist(strsplit(unlist(strsplit(ff,'_'))[2],'[.]'))[1])
    if ((unlist(strsplit(unlist(strsplit(ff,'_'))[2],'[.]'))[1]=='Lombardia')) {minlen<-nrow(readMat(ff)$x);print(ii);print(minlen)}
}
regions <- unique(regions)
setwd(thisishere)
pertnum <- (length(fordd)/length(regions))


###combination of results list
comb <- function(x, ...) {
  lapply(seq_along(x),
#   pertnum,
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

#results, data, deaddiff, %

ResultsList <- foreach(rr = regions, .combine='comb', .multicombine=TRUE,
                      .init= list(list(),list(),list(), list(),list()) ) %dopar% {

library(data.table)
library(readxl)
library(WriteXLS)
library(R.matlab)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(tikzDevice)
library(latex2exp)

    quar <- hosp <- icu <- dead <- NULL
    quarDeC <- hospDeC <- icuDeC <- deadDeC <- NULL
    final_Dead <- SavedPerc <- NULL

    for (ii in 1:pertnum) {
        
        ff <- fordd[ii]
            tmpplo <- readMat(paste0(thisishere,'/cap_res/',ii,'_',rr,'.mat'))
            tmpploDeC <- readMat(paste0(thisishere,'/deCap_res/',ii,'_',rr,'.mat'))
        
            AA<-tmpplo$x[,5]+tmpplo$x[,6]+tmpplo$x[,15]
            BB<-tmpplo$x[,7]+tmpplo$x[,8]
            CC<-tmpplo$x[,9]+tmpplo$x[,10]
            DDcap<-tmpplo$x[,12]
            
            { quar <- cbind(quar,(AA)); hosp <- cbind(hosp,(BB)); icu <- cbind(icu,(CC)); dead <- cbind(dead,DDcap) }
            
            
            ##deCap_res
            AA<-tmpploDeC$x[,5]+tmpploDeC$x[,6]+tmpploDeC$x[,15]
            BB<-tmpploDeC$x[,7]+tmpploDeC$x[,8]
            CC<-tmpploDeC$x[,9]+tmpploDeC$x[,10]
            DDdeC<-tmpploDeC$x[,12]
            
            { quarDeC <- cbind(quarDeC,(AA)); hospDeC <- cbind(hospDeC,(BB)); icuDeC <- cbind(icuDeC,(CC)); deadDeC <- cbind(deadDeC,DDdeC) }
            
            DD <- DDcap - DDdeC
            final_Dead <- cbind(final_Dead, DD)
            SavedPerc <- cbind(SavedPerc,DD[length(DD)]*100/DDcap[length(DDcap)])
    }
    
    miss <- minlen-nrow(quar)
    print(miss)
    A <- rowMeans(quar); B <- rowMeans(hosp); C <- rowMeans(icu); D <- rowMeans(dead)
    As <- rowSds(quar); Bs <- rowSds(hosp); Cs <- rowSds(icu); Ds <- rowSds(dead)
    
    DsDeC <- rowMeans(deadDeC)
    
    
    tt <- seq(miss,(minlen-1))
    rtmp <- rep(rr,length(A))
    collect_mean <- data.frame(cbind(tt,A,B,C,D,rtmp))
    collect_sd <- data.frame(cbind(tt,As,Bs,Cs,Ds,rtmp))
    ###finaldead
    final_Dead <- cbind(tt,final_Dead,rtmp)
    PercofMean <- (mean(dead[nrow(dead),])-mean(deadDeC[nrow(deadDeC),]))*100/mean(dead[nrow(dead),])
    PercofMean <- data.frame(cbind(PercofMean,rr))
    colnames(PercofMean)<- c('Perc','Area')
    ###
    ##Data
    dataReal <- read_excel(paste0(thisishere,'/../Data/',rr,".xls"))
    dataReal <- dataReal[1:length(tt),]
    qq <- dataReal[,which(colnames(dataReal)=='qua')]; hh <- dataReal[,which(colnames(dataReal)=='hos')]; iic <- dataReal[,which(colnames(dataReal)=='icu')]; dd <- dataReal[,which(colnames(dataReal)=='dead')];
    data_coll <- data.frame(cbind(tt,qq,hh,iic,dd,rtmp))
    
    hh <- readMat(paste0(thisishere,'/rt_cap/1_',rr,'.mat'))
    hh <- hh$parsave[3:nrow(hh$parsave), 16:17]
    
    hlim <- data.frame( cbind(tt[1:nrow(hh)], hh, rep(rr,nrow(hh))) ) 
    
    MINHH <- data.frame( cbind( min(hh[,1]+hh[,2]), rr ) )
    MAXHH <- data.frame( cbind( max(hh[,1]+hh[,2]), rr ) )
    
    setwd(thisishere)
    
    
    colnames(collect_mean) <- colnames(collect_sd) <- colnames(data_coll) <- c('Time','Active Infection','Hospitalized','ICU','Dead','Area')

    colnames(hlim) <- c('Time','Hospitalized','ICU','Area')

    for (jj in 1:(ncol(hlim)-1)) {
        hlim[,jj] <- as.numeric(as.character(hlim[,jj]))
    }
    hlim <- melt(hlim, id.vars=c('Time','Area'))

    colnames(MINHH) <- colnames(MAXHH) <- c('Limit','Area')
    MINHH <- melt(MINHH, id.vars=c('Area'))
    MAXHH <- melt(MAXHH, id.vars=c('Area'))
    MINHH$value <- as.numeric(as.character(MINHH$value))
    MAXHH$value <- as.numeric(as.character(MAXHH$value))

    for (jj in 1:(ncol(data_coll)-1)) {
        data_coll[,jj] <- as.numeric(as.character(data_coll[,jj]))
        collect_mean[,jj] <- as.numeric(as.character(collect_mean[,jj]))
        collect_sd[,jj] <- as.numeric(as.character(collect_sd[,jj]))
    }

    ccm <- melt(collect_mean,id.vars=c('Time','Area'))
    ccs <- melt(collect_sd,id.vars=c('Time','Area'))
    ddc <- melt(data_coll,id.vars=c('Time','Area'))

    ccm$sd <- ccs$value

    ccm$value <- as.numeric(as.character(ccm$value))
    ccm$sd <- as.numeric(as.character(ccm$sd))

    names(ccm)[3] <- names(ddc)[3] <- names(hlim)[3] <- 'Legend'

    levels(ccm$Time) <- levels(ddc$Time) <- as.Date( as.numeric(as.character(levels(ccm$Time))), origin='2020-02-24',by='days')    
    
    ###
    
    #### 
#     final_Dead <- final_Dead[which(final_Dead[,1]==max(tt)),]
   
    final_Dead <- data.frame(final_Dead)
    colnames(final_Dead) <- c('Time',seq(1,(ncol(final_Dead)-2)),'Area')
    for (jj in 1:(ncol(final_Dead)-1)) {
        final_Dead[,jj] <- as.numeric(as.character(final_Dead[,jj]))
    }
    vf <- melt(final_Dead,id.vars=c('Time','Area'))
    ###
    ###percentage
    SavedPerc <- data.frame( cbind(paste0(round(SavedPerc,digits=2),'%'),rr) )
    colnames(SavedPerc) <- c('Perc','Area')
    #togeth <- data.frame(cbind(rbind(Ds[length(Ds)],DsDeC[length(DsDeC)]),rbind('cap','decap'),rbind(rr,rr)))
    #colnames(togeth) <- c('value','Type','Area')
    #tmp <- togeth[which(togeth$Area==rr),]
    #tmp <- tmp[which(tmp$Time==max(tmp$Time)),]
#     deca <- DsDeC[length(DsDeC)]#(tmp[which(tmp$Type=='decap'),])$value
#     cap <- Ds[length(Ds)]#tmp[which(tmp$Type=='cap'),]$value
#     
#     savepe <- (cap-deca)*100/cap
#     perce_coll <- data.frame(cbind(paste0(round(savepe,digits=2),'%'),rr) )
#     pp_coll <- data.frame(cbind(round(savepe,digits=2),rr ))
#     colnames(perce_coll) <- colnames(pp_coll) <- c('Perc','Area')
#     pp_coll$Perc <- as.numeric(as.character(pp_coll$Perc))
#     ###
    
    
    
    list(ccm,PercofMean,vf,ddc,hlim)
    #list(capacityResults, %,finaldead,data,)
}

CapModSimul<-PercResul<-DeadDiff<-DataColl<-HospLim<-NULL
for (RL in 1:length(ResultsList[[1]])) {
    CapModSimul <- rbind(CapModSimul, as.data.frame(ResultsList[[1]][[RL]]))
    PercResul <- rbind(PercResul, as.data.frame(ResultsList[[2]][[RL]]))
    DeadDiff <- rbind(DeadDiff, as.data.frame(ResultsList[[3]][[RL]]))
    DataColl <- rbind(DataColl, as.data.frame(ResultsList[[4]][[RL]]))
    HospLim <- rbind(HospLim, as.data.frame(ResultsList[[5]][[RL]]))
}

levels(CapModSimul$Area)[which(levels(CapModSimul$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')
levels(PercResul$Area)[which(levels(PercResul$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')
levels(DeadDiff$Area)[which(levels(DeadDiff$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')
levels(DataColl$Area)[which(levels(DataColl$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')
levels(HospLim$Area)[which(levels(HospLim$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')


vf <- DeadDiff[which(DeadDiff$Time==max(DeadDiff$Time)),]
selec<-where<-NULL
for (rrr in levels(vf$Area)) {#
    ttmp <- vf[which(vf$Area==rrr),]
    levels(ttmp$Area) <- droplevels(ttmp$Area)
    where <- c(where, max(ttmp$value))
#     print(summary(ttmp))
    if (median(ttmp$value)>10) {
        selec <- c(selec, rrr)
    }
}
supplselec <- setdiff(unique(vf$Area),selec)
PercResul$where <-  round(as.numeric(as.character(where)))+1

PercResul$Perc <- paste0(round(as.numeric(as.character(PercResul$Perc)),digits=2),'%')

###
# pp_coll <-
vfsel <- vf[which(vf$Area %in% selec),]
vfsuppl <- vf[which(vf$Area %in% supplselec),]

mainCapSim <- CapModSimul[which(CapModSimul$Area %in% selec),]
suppCapSim <- CapModSimul[which(CapModSimul$Area %in% supplselec),]
mainData <- DataColl[which(DataColl$Area %in% selec),]
suppData <- DataColl[which(DataColl$Area %in% supplselec),]

ppfit1 <- ggplot(mainCapSim[which(mainCapSim$Legend %in% c('Hospitalized','ICU')),] , aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Legend)) +
geom_line() +
geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill=Legend),col=NA, alpha=0.5) +
geom_point(mainData[which(mainData$Legend %in% c('Hospitalized','ICU')),] , mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Legend), size=0.5) +
geom_line(HospLim[which(HospLim$Area %in% selec),], mapping=aes(y=value, linetype=Legend) ) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
scale_y_sqrt(breaks = function(x) { c(100,500,pretty(seq(min(x),max(x),100))) } ) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Cases') + ggtitle('A')


ppfitsuppl <- ggplot(suppCapSim[which(suppCapSim$Legend %in% c('Hospitalized','ICU')),], aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Legend)) +
geom_line() +
geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill=Legend),col=NA, alpha=0.5) +
geom_point(suppData[which(suppData$Legend %in% c('Hospitalized','ICU')),], mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Legend), size=0.5) +
geom_line(HospLim[which(HospLim$Area %in% supplselec),], mapping=aes(y=value) ) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
scale_y_sqrt(breaks = function(x) { c(100,500,pretty(seq(min(x),max(x),100))) } ) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Cases') + ggtitle('A')

mainPS <- PercResul[which(PercResul$Area %in% selec),]
suppPS <- PercResul[which(PercResul$Area %in% supplselec),]

ppbp <- ggplot(vfsel, aes(x=Area, y=value)) +
geom_boxplot(aes(x=Area,y=value),outlier.shape=NA) +
scale_y_sqrt(breaks = c( 5,50, seq(100,500,100), seq(1000,7000,1000) ) ) + #5,25,40,70
geom_text(mainPS, mapping=aes(label = Perc, y=where+50 ),show.legend = FALSE) +
#facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Excess Dead')+ ggtitle('B')

supplbp <- ggplot(vfsuppl, aes(x=Area, y=value)) +
geom_boxplot(aes(x=Area,y=value),outlier.shape=NA) +
scale_y_sqrt(breaks = c( 5,50, seq(100,500,100), seq(1000,7000,1000) ) ) + #5,25,40,70
geom_text(suppPS, mapping=aes(label = Perc, y=where+50 ),show.legend = FALSE) +
#facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Excess Dead')+ ggtitle('B')


stopCluster(cl)

ggsave(
  "fit_capSuppl.pdf",
  plot = ppfitsuppl,
  width = 15,
  height = 7,
  units = "in",
  dpi = 300
)

print('Supplementary Fit plot done ...')

ggsave(
  "excedead.pdf",
  plot = grid.arrange(ppfit1,ppbp),
  width = 11,
  height = 11,
  units = "in",
  dpi = 300
)

print('Main Fit and Excess Dead plot done ...')

ggsave(
  "excesuppl.pdf",
  plot = grid.arrange(ppfitsuppl,supplbp),
  width = 11,
  height = 11,
  units = "in",
  dpi = 300
)

print('Supplementary Fit and Excess Dead plot done.')




