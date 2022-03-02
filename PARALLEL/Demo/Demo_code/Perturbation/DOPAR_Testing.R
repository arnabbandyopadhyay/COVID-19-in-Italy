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

setwd('Refer_pert_res')
RefeFolder <- getwd()

fordd <- list.files(pattern = "\\.mat$")
regions <- list()
for (ii in 1:length(fordd)) {
    ff <- fordd[ii]
    regions[[ii]] <- (unlist(strsplit(unlist(strsplit(ff,'_'))[2],'[.]'))[1])
    if ((unlist(strsplit(unlist(strsplit(ff,'_'))[2],'[.]'))[1])=='Italy') {minlen<-nrow(readMat(ff)$x)}
}
regions <- unique(regions)

pertnum <- (length(fordd)/length(regions))
setwd(thisishere)

###combination of results list
comb <- function(x, ...) {
  lapply(seq_along(x),
#   pertnum,
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


#list(simulDynamic, textwh%, Rtdynamic)
ResultsList <- foreach(rr = regions, .combine='comb', .multicombine=TRUE,
                      .init= list(list(),list(),list(),list()) ) %dopar% {
    
    library(data.table)
    library(readxl)
    library(WriteXLS)
    library(R.matlab)
    library(ggplot2)
    library(matrixStats)
    library(gridExtra)
    library(tikzDevice)
    library(latex2exp)
    
    hospbase <- hospredEA<-hospredLA <- deadredEA <-deadredLA<- deadbase <- NULL
    rt_EA <- rt_LA <- rtM <- NULL
    
    for (pp in 1:pertnum) {
        base <- readMat(paste0(thisishere,'/Refer_pert_res/',pp,'_',rr,'.mat') )
        reducearly <- readMat(paste0(thisishere,'/Testing_pert_res/red_early_',pp,'_',rr,'.mat')  )
        reduclate <- readMat(paste0(thisishere,'/Testing_pert_res/red_late_',pp,'_',rr,'.mat')  )
        tmpb <- ( base$x[,7]+base$x[,8]+base$x[,9]+base$x[,10] )
        tmprEA <-  ( reducearly$x[,7]+reducearly$x[,8]+reducearly$x[,9]+reducearly$x[,10])
        tmprLA <-  ( reduclate$x[,7]+reduclate$x[,8]+reduclate$x[,9]+reduclate$x[,10])

        hospbase <- cbind(hospbase, tmpb )
        hospredEA <- cbind(hospredEA, tmprEA)
        hospredLA <- cbind(hospredLA, tmprLA)
        
        deadbase <- cbind(deadbase, base$x[,12])
        deadredEA <- cbind(deadredEA, reducearly$x[,12])
        deadredLA <- cbind(deadredLA, reduclate$x[,12])
        
        ###Rt
        tmpplo <- readMat(paste0( 'Testing_rt_res/early_',pp,'_',rr,'.mat' )) #parsave140x15 #rtsaave140x1
        rt_EA <- cbind(rt_EA,tmpplo$rtsave)
            
        tmpplo <- readMat(paste0( 'Testing_rt_res/late_',pp,'_',rr,'.mat' )) #parsave140x15 #rtsaave140x1
        rt_LA <- cbind(rt_LA,tmpplo$rtsave)
            
        tmpplo <- readMat(paste0('Refer_rt_res/',pp,'_',rr,'.mat'  )) #parsave140x15 #rtsaave140x1
        rtM <- cbind(rtM,tmpplo$t)
        ###Rt mean sd
        missRt <- minlen-nrow(rt_EA)-5
    
        A <- rowMeans(rt_EA); B <- rowMeans(rt_LA)  #B <- rowMeans(hosp);
        As <- rowSds(rt_EA); Bs <- rowSds(rt_LA) #Bs <- rowSds(hosp);
        C <- rowMeans(rtM); Cs <- rowSds(rtM)
        ttPar <- seq(missRt,(nrow(rt_EA)+missRt-1))#seq(miss,(minlen-1))
        rtmp <- rep(rr,length(A))
        coll_rtmean <- ( cbind(ttPar,A,B,C,rtmp))
        coll_rtsd <- ( cbind(ttPar,As,Bs,Cs,rtmp))
        ###
    }
    
    hhB <- rowMeans(hospbase); hhBsd <- rowSds(hospbase)
    hhREA <- rowMeans(hospredEA); hhREAsd <- rowSds(hospredEA)
    hhRLA <- rowMeans(hospredLA); hhRLAsd <- rowSds(hospredLA)
    ddB <- rowMeans(deadbase); ddBsd <- rowSds(deadbase)
    ddREA <- rowMeans(deadredEA); ddREAsd <- rowSds(deadredEA)
    ddRLA <- rowMeans(deadredLA); ddRLAsd <- rowSds(deadredLA)
    
    
    miss <- minlen-length(hhB)
    print(miss)
    tt <- seq(miss,(minlen-1))
    rtmp <- rep(rr,length(tt))
    
    BASE <- cbind(tt,hhB,ddB,rtmp)
    sdsbase <- cbind(tt,hhBsd,ddBsd,rtmp)
    
    RED_EARLY <- cbind(tt,hhREA,ddREA,rtmp)
    
    sdsred_early <- cbind(tt,hhREAsd,ddREAsd,rtmp)
    RED_LATE <- cbind(tt,hhRLA,ddRLA,rtmp)
    sdsred_late <- cbind(tt,hhRLAsd,ddRLAsd,rtmp)
    
    HE <- round((max(hhB)-max(hhREA))/max(hhB)*100)#, digits=2)
    HL <- round((max(hhB)-max(hhRLA))/max(hhB)*100)#, digits=2)
    DE <- round((max(ddB)-max(ddREA))/max(ddB)*100)#, digits=2)
    DL <- round((max(ddB)-max(ddRLA))/max(ddB)*100)#, digits=2)
    
    keepD <- rbind( cbind( DE, rr, (max(ddREA)),'Testing -- Early' ),cbind( DL, rr, (max(ddRLA)),'Testing -- Late' )  )
    keepH <- rbind( cbind( HE, rr, (max(hhREA)) ,'Testing -- Early' ),cbind( HL, rr, (max(hhRLA)) ,'Testing -- Late' )  )
   
    collDeadCap <- ( keepD )
    collHosCap <- ( keepH )
    
    BASE <- data.frame(BASE)
    RED_EARLY <- data.frame(RED_EARLY)
    RED_LATE <- data.frame(RED_LATE)

    sdsbase <- data.frame(sdsbase)
    sdsred_early <- data.frame(sdsred_early)
    sdsred_late <- data.frame(sdsred_late)

    collDeadCap <- data.frame(collDeadCap)
    collHosCap <- data.frame(collHosCap)

    for (jj in 1:(ncol(BASE)-1)) {
        BASE[,jj] <- as.numeric(as.character(BASE[,jj]))
        RED_EARLY[,jj] <- as.numeric(as.character(RED_EARLY[,jj]))
        RED_LATE[,jj] <- as.numeric(as.character(RED_LATE[,jj]))
        sdsbase[,jj] <- as.numeric(as.character(sdsbase[,jj]))
        sdsred_early[,jj] <- as.numeric(as.character(sdsred_early[,jj]))
        sdsred_late[,jj] <- as.numeric(as.character(sdsred_late[,jj]))
    }

    for (jj in c(1,3)) {
        collDeadCap[,jj] <- as.numeric(as.character(collDeadCap[,jj]))
        collHosCap[,jj] <- as.numeric(as.character(collHosCap[,jj]))
    }


    colnames(collDeadCap) <- colnames(collHosCap) <- c('value','Area','where','Model')
    colnames(BASE) <- colnames(RED_EARLY) <- colnames(RED_LATE) <- colnames(sdsbase)<- colnames(sdsred_early) <- colnames(sdsred_late) <- c('Time','Hospitalized','Dead','Area')

    
    basemelt <- melt(BASE, id.vars=c('Time','Area') )
    redmelt_early <- melt(RED_EARLY, id.vars=c('Time','Area') )
    redmelt_late <- melt(RED_LATE, id.vars=c('Time','Area'))
    sdsbase <- melt(sdsbase,id.vars=c('Time','Area'))
    sdsred_early <- melt(sdsred_early, id.vars=c('Time','Area') )
    sdsred_late <- melt(sdsred_late, id.vars=c('Time','Area'))

    basemelt$sds <- sdsbase$value
    redmelt_early$sds <- sdsred_early$value
    redmelt_late$sds <- sdsred_late$value

    basemelt$Model <- rep('Reference',nrow(basemelt) )
    redmelt_early$Model <- rep('Testing -- Early', nrow(redmelt_early) )
    redmelt_late$Model <- rep('Testing -- Late', nrow(redmelt_late))
    collect <- rbind(basemelt,redmelt_early,redmelt_late)

    collect$Time <- as.numeric(as.character(collect$Time))
    collect$sds <- as.numeric(as.character(collect$sds))
    collect$value <- as.numeric(as.character(collect$value))

    names(collect)[3] <- 'Legend'
    
    ###Rt
        
    coll_rtmean <- data.frame(coll_rtmean)
    coll_rtsd <- data.frame(coll_rtsd)

    for (jj in 1:(ncol(coll_rtmean)-1)) {
        coll_rtmean[,jj] <- as.numeric(as.character(coll_rtmean[,jj]))
        coll_rtsd[,jj] <- as.numeric(as.character(coll_rtsd[,jj]))
    }

    colnames(coll_rtmean) <- colnames(coll_rtsd) <- c('Time','Testing -- Early','Testing -- Late','Reference','Area')

    rtmm <- melt(coll_rtmean,id.vars=c('Time','Area'))
    rtss <- melt(coll_rtsd,id.vars=c('Time','Area'))

    rtmm$sd <- rtss$value

    names(rtmm)[3] <- 'Model'

    rtmm$ci <- (rtmm$sd*1.96)/sqrt(pertnum) 


    ####

    list(collect,collHosCap,collDeadCap,rtmm)

}


collSimul<-collHosCap<-collDeadCap<-RtSimul<-NULL
for (RL in 1:length(ResultsList[[1]])) {
    collSimul <- rbind(collSimul, as.data.frame(ResultsList[[1]][[RL]]))
    collHosCap <- rbind(collHosCap, as.data.frame(ResultsList[[2]][[RL]]))
    collDeadCap <- rbind(collDeadCap, as.data.frame(ResultsList[[3]][[RL]]))
    RtSimul <- rbind(RtSimul, as.data.frame(ResultsList[[4]][[RL]]))
}


levels(collSimul$Area)[which(levels(collSimul$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')
# 
levels(collHosCap$Area)[which(levels(collHosCap$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')

levels(collDeadCap$Area)[which(levels(collDeadCap$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')

levels(RtSimul$Area)[which(levels(RtSimul$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')




sel <- c('Italy','Lombardia','Veneto')
supplsel <- setdiff(regions,sel)


collmain <- collSimul[which(collSimul$Area %in% sel),]

collsuppl <- collSimul[which(collSimul$Area %in% supplsel),]

RtMain <- RtSimul[which(RtSimul$Area %in% sel),]

ppfithos <- ggplot(collmain[which(collmain$Legend=='Hospitalized'),], aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Model)) +
geom_line() +
geom_ribbon(aes(ymin=value-sds, ymax=value+sds, fill=Model),col=NA, alpha=0.5) +
geom_text(collHosCap[which(collHosCap$Area %in% collmain$Area),], mapping=aes(x=as.Date(165, origin='2020-02-24',by='days'), y=where, label=paste0(value,'%') ) , show.legend = FALSE)+
#geom_point(ddcs, mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Legend), size=0.5) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e",expand = expansion(add = 2)) +
scale_y_sqrt(breaks = function(x) { pretty(seq(min(x),max(x),1000)) } ) +
expand_limits(x = as.Date('2020-08-21'),y=10) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Cases') + ggtitle('A')


ppfitdead <- ggplot(collmain[which(collmain$Legend=='Dead'),], aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Model)) +
geom_line() +
geom_ribbon(aes(ymin=value-sds, ymax=value+sds, fill=Model),col=NA, alpha=0.5) +
geom_text(collDeadCap[which(collDeadCap$Area %in% collmain$Area),], mapping=aes(x=as.Date(165, origin='2020-02-24',by='days'), y=where+c(-2,2), label=paste0(value,'%') ), show.legend = FALSE) +
#geom_point(ddcs, mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Legend), size=0.5) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
expand_limits(x = as.Date('2020-08-21')) +
scale_y_sqrt(breaks = function(x) { pretty(seq(min(x),max(x),1000)) } ) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Cases') + ggtitle('B')


pp_rt <- ggplot(RtMain, aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Model)) +
geom_line() +
#geom_ribbon(aes(ymin=value-ci, ymax=value+ci, fill=Model),col=NA, alpha=0.5) +
geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill=Model),col=NA, alpha=0.5) +
geom_hline(aes(yintercept=1),col='black') +
geom_vline(aes(xintercept=as.Date('2020-03-08')),col='dark red') +
geom_vline(aes(xintercept=as.Date('2020-05-17')),col='dark green') +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
#             scale_y_sqrt(breaks = function(x) { seq(100,max(x),100) } ) +
facet_wrap(~Area, scale='free_y') +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 25)) +
xlab('') + ylab( expression( paste('Reproduction Number ', italic(R)[t]) ) ) 

stopCluster(cl)




ggsave(
  'TestMod_fit_Main.pdf',
  plot = grid.arrange(ppfithos,ppfitdead,ncol=1),
  width = 11,
  height = 7,
  units = "in",
  dpi = 300
)

print('Main Test Plot done ...')

if (length(supplsel)>0) {
ggsave(
  "TestMod_Rt_Supp.pdf",
  plot = pp_rt,
  width = 22,
  height = 6,
  units = "in",
  dpi = 300
)

print('Supplementary Test Plot done.')
}