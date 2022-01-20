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
                      .init= list(list(),list(),list(),list(),list(),list()) ) %dopar% {

library(data.table)
library(readxl)
library(WriteXLS)
library(R.matlab)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(tikzDevice)
library(latex2exp)

##############################################################
################################################ 005

#################################### EARLY LATE TESTING MAXCAP 005

    ###
    hospit <- dead <- NULL
    hospLate<-deadLate<-hospEarly<-deadEarly<-deadcap<-hoscap<-deadcapM<-hoscapM<-NULL
    for (ii in 1:pertnum) {
        
        ###LATE TESTING
        tmpplo <- readMat(paste0(thisishere,'/LateTest_res/',ii,'_',rr,'_005.mat'))
        AAd<-tmpplo$x[,5]+tmpplo$x[,6]+tmpplo$x[,15]+tmpplo$x[,16]
        BBd<-tmpplo$x[,7]+tmpplo$x[,8]
        CCd<-tmpplo$x[,9]+tmpplo$x[,10]
        DDd<-tmpplo$x[,12]
        
        { hospLate <- cbind(hospLate,(BBd+CCd)) }
        
        { deadLate <- cbind(deadLate,DDd) }
        
        ###EARLY TESTING
        tmpplo <- readMat(paste0(thisishere,'/EarlyTest_res/',ii,'_',rr,'_005.mat'))
        AAc<-tmpplo$x[,5]+tmpplo$x[,6]+tmpplo$x[,15]+tmpplo$x[,16]
        BBc<-tmpplo$x[,7]+tmpplo$x[,8]
        CCc<-tmpplo$x[,9]+tmpplo$x[,10]
        DDc<-tmpplo$x[,12]
        
        ll<-min(length(DDc),length(DDd))
        DD <- DDc-DDd#[1:ll] - DDd[1:ll]
        BB <- (BBc+CCc) - (BBd+CCd)
        
        { hospEarly <- cbind(hospEarly,(BBc+CCc)) }
        
        { deadEarly <- cbind(deadEarly,DDc) }
        
        ###NORMAL CAP
        tmpplo <- readMat(paste0(thisishere,'/cap_res/',ii,'_',rr,'.mat'))
#         tmpplo <- readMat(paste0(thisishere,'/deCap_res/',ii,'_',rr,'.mat'))
        
        HHcap <- tmpplo$x[,7]+tmpplo$x[,8]+tmpplo$x[,9]+tmpplo$x[,10]
        DDcap <- tmpplo$x[,12]
        
        deadcap <- cbind(deadcap, DDcap)
        hoscap <- cbind(hoscap, HHcap)
    }
    
    ######means and sd
    miss <- minlen-nrow(deadcap)
    print(miss)
    
    tt <- seq(miss,(minlen-1))
    rtmp <- rep(rr,length(tt) )
    

    ###Early Testing + MaxCap
    hosEA <- rowMeans(hospEarly); deadCP <- rowMeans(deadEarly)
    hosSC <- rowSds(hospEarly); deadSC <- rowSds(deadEarly)

    capmean <- cbind(tt, rtmp, rep('TestCap -- Early',length(hosEA)), hosEA )
    capsd <- cbind(tt, rtmp, rep('TestCap -- Early',length(hosEA)), hosSC )
    
    deadcapm <- cbind(tt, rtmp, rep('TestCap -- Early',length(hosEA)), deadCP )
    deadcapsd <- cbind(tt, rtmp, rep('TestCap -- Early',length(hosEA)), deadSC )
    ############################
    ###Late Testing + MaxCap
    hosT <- rowMeans(hospLate); deadT <- rowMeans(deadLate)
    hosST <- rowSds(hospLate); deadST <- rowSds(deadLate)

    testmean <- cbind(tt, rtmp, rep('TestCap -- Late',length(hosEA)), hosT )
    testsd <- cbind(tt, rtmp, rep('TestCap -- Late',length(hosEA)), hosST )
    
    deadtestm <- cbind(tt, rtmp, rep('TestCap -- Late',length(hosEA)), deadT )
    deadtestsd <- cbind(tt, rtmp, rep('TestCap -- Late',length(hosEA)), deadST )
    #######################
    ###reference cap
    deadcapM <- rowMeans(deadcap); deadcapsd <- rowSds(deadcap)
    hoscapM <- rowMeans(hoscap); hoscapsds <- rowSds(hoscap)
  
    caprefe <- cbind(tt, rtmp, rep('Capacity',length(hoscap)), hoscapM )
    caprefesd <- cbind(tt, rtmp, rep('Capacity',length(hoscap)), hoscapsds )
    
    deadcaprefe <- cbind(tt, rtmp, rep('Capacity',length(hosEA)), deadcapM )
    deadcaprefesd <- cbind(tt, rtmp, rep('Capacity',length(hosEA)), deadcapsd )
    
    
    #################
    HE <- round((max(hoscapM)-max(hosEA))/max(hoscapM)*100)#, digits=2)
    HL <- round((max(hoscapM)-max(hosT))/max(hoscapM)*100)#, digits=2)
    DE <- round((max(deadcapM)-max(deadCP))/max(deadcapM)*100)#, digits=2)
    DL <- round((max(deadcapM)-max(deadT))/max(deadcapM)*100)#, digits=2)
    
    keepD <- rbind( cbind( DE, rr, (max(deadCP)),'TestCap -- Early' ),cbind( DL, rr, (max(deadT)),'TestCap -- Late' )  )
    keepH <-rbind( cbind( HE, rr, (max(hosEA)) ,'TestCap -- Early' ),cbind( HL, rr, (max(hosT)) ,'TestCap -- Late' )  )
    collDeadCap <- data.frame(( keepD ))
    collHosCap <- data.frame(( keepH ))
    
    for (jj in c(1,3)) {
        collDeadCap[,jj] <- as.numeric(as.character(collDeadCap[,jj]))
        collHosCap[,jj] <- as.numeric(as.character(collHosCap[,jj]))
    }

    colnames(collDeadCap) <- colnames(collHosCap) <- c('value','Area','where','Model')

    ###########
    
    ###
    meancollect <- data.frame( rbind( cbind(capmean, deadCP), cbind(testmean, deadT), cbind(caprefe,deadcapM) ) )
    
    sdcollect <- data.frame( rbind( cbind(capsd, deadSC), cbind(testsd, deadST), cbind(caprefesd,deadcapsd) ) )
    ############
    
    ######Data
    dataReal <- read_excel(paste0(thisishere,'/../Data/',rr,".xls"))
    dataReal <- dataReal[1:length(tt),]
    hosdata <- dataReal[,which(colnames(dataReal)=='hos')]; iicdata <- dataReal[,which(colnames(dataReal)=='icu')]; dd <- dataReal[,which(colnames(dataReal)=='dead')];
    
    data_coll <- data.frame(cbind(tt,hosdata+iicdata,dd,rtmp))
    ####
    for (jj in c(1,4,5)) {
        meancollect[,jj] <- as.numeric(as.character(meancollect[,jj]))
        sdcollect[,jj] <- as.numeric(as.character(sdcollect[,jj]))
    }

    for (jj in 1:(ncol(data_coll)-1)) {
        data_coll[,jj] <- as.numeric(as.character(data_coll[,jj]))
    }

    colnames(meancollect) <- c('Time','Area','Model','Hospitalized','Dead')
    colnames(sdcollect) <- c('Time','Area','Model','Hospitalized','Dead')
    colnames(data_coll) <- c('Time','Hospitalized','Dead','Area')
    
    datamelt <- melt(data_coll, id.vars=c('Time','Area'))

    captest_collect <- melt(meancollect, id.vars=c('Time','Area','Model'))
    sdtest_collect <- melt(sdcollect, id.vars=c('Time','Area','Model'))
    captest_collect$sd <- sdtest_collect$value
    
    ###hosp capacities
    hh <- readMat(paste0(thisishere,'/rt_cap/1_',rr,'.mat'))
    hh <- hh$parsave[3:nrow(hh$parsave), 16:17]
#     hlim <- rbind( hlim, cbind(tt, hh, rep(rr,nrow(hh))) )
    hlim <- data.frame( cbind(tt[1:nrow(hh)], hh, rep(rr,nrow(hh))) ) 
    MINHH <- data.frame( cbind( min(hh[,1]+hh[,2]), rr ) )
    MAXHH <- data.frame( cbind( max(hh[,1]+hh[,2]), rr ) )
    colnames(hlim) <- c('Time','Hospitalized','ICU','Area') 
    colnames(MINHH) <- colnames(MAXHH) <- c('Limit','Area')

    MINHH$Limit <- as.numeric(as.character(MINHH$Limit))

    MAXHH$Limit <- as.numeric(as.character(MAXHH$Limit))
    
    for (jj in 1:(ncol(hlim)-1)) {
        hlim[,jj] <- as.numeric(as.character(hlim[,jj]))
    }
    hlim <- melt(hlim, id.vars=c('Time','Area'))
    
    MINHH <- melt(MINHH, id.vars=c('Area'))
    MAXHH <- melt(MAXHH, id.vars=c('Area'))
    ###
    list(captest_collect,datamelt,collHosCap,collDeadCap,MINHH,MAXHH)
}


captest_collect<-datamelt<-collHosCap<-collDeadCap<-MINHH<-MAXHH<-NULL
for (RL in 1:length(ResultsList[[1]])) {
    captest_collect <- rbind(captest_collect, as.data.frame(ResultsList[[1]][[RL]]))
    datamelt <- rbind(datamelt, as.data.frame(ResultsList[[2]][[RL]]))
    collHosCap <- rbind(collHosCap, as.data.frame(ResultsList[[3]][[RL]]))
    collDeadCap <- rbind(collDeadCap, as.data.frame(ResultsList[[4]][[RL]]))
    MINHH <- rbind(MINHH, as.data.frame(ResultsList[[5]][[RL]]))
    MAXHH <- rbind(MAXHH, as.data.frame(ResultsList[[6]][[RL]]))
}


levels(captest_collect$Area)[which(levels(captest_collect$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')

levels(datamelt$Area)[which(levels(datamelt$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')

levels(collHosCap$Area)[which(levels(collHosCap$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')
levels(collDeadCap$Area)[which(levels(collDeadCap$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')

levels(MINHH$Area)[which(levels(MINHH$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')

levels(MAXHH$Area)[which(levels(MAXHH$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')


sel <- c('Emilia Romagna','Lombardia','Piemonte')
supplsel <- setdiff(regions,sel)

HH <-captest_collect[which(captest_collect$variable=='Hospitalized'),]
HH <- HH[which(HH$Area %in% sel),]
DD <- captest_collect[which(captest_collect$variable=='Dead'),]
DD <- DD[which(DD$Area %in% sel),]

dmsel <- datamelt[which(datamelt$Area %in% sel ),]
dmselh <- dmsel[which(dmsel$variable=='Hospitalized'),]
dmseld <- dmsel[which(dmsel$variable=='Dead'),]

MAXHHmain <- MAXHH[which(MAXHH$Area %in% sel),]
MAXHHsupp <- MAXHH[which(MAXHH$Area %in% supplsel),]
MINHHmain <- MINHH[which(MINHH$Area %in% sel),]

collselhos <- collHosCap[which(collHosCap$Area %in% HH$Area),]
collseldead <- collDeadCap[which(collDeadCap$Area %in% DD$Area),]


ppfithos <- ggplot(HH, aes(x=as.Date(Time, origin='2020-02-24',by='days',format='%Y-%m-%d'), y=value, col=Model)) +
geom_line() +
geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill=Model),col=NA, alpha=0.5) +
geom_point(dmselh, mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value), col='black', size=0.1) +
geom_hline(MINHHmain, mapping=aes(yintercept=value ) ) + #5,25,40,70
# geom_hline(MAXHHmain, mapping=aes(yintercept=value ) ) + #5,25,40,70
#geom_text(collselhos, mapping=aes(label = value, y=where, x=as.Date(100, origin='2020-02-24',by='days',format='%Y-%m-%d') ,col=Model)) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
#scale_y_sqrt(breaks = function(x) { pretty(seq(min(x),max(x),100)) } ) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Hospitalized') + ggtitle('A')

updown <- c(-1000,0, -5000,-4000, -300,-200, -150,0)
ppfitdead <- ggplot(DD, aes(x=as.Date(Time, origin='2020-02-24',by='days',format='%Y-%m-%d'), y=value, col=Model)) +
geom_line() +
geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill=Model),col=NA, alpha=0.5) +
geom_point(dmseld, mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value), col='black', size=0.1) +
geom_text(collseldead, mapping=aes(label = paste0(value,'%'), y=where,
#+updown, 
x=as.Date(95, origin='2020-02-24',by='days',format='%Y-%m-%d') ,col=Model), show.legend = FALSE) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
expand_limits(x = as.Date('2020-06-07')) +
scale_y_sqrt(breaks = function(x) { pretty(seq(min(x),max(x),100)) } ) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Dead') + ggtitle('B')

#pdf('TestMAXCap_005.pdf',h=11,w=11)
#grid.arrange(ppfithos,ppfitdead,ncol=1)
#dev.off()

ggsave(
  'TestMAXCap_005.pdf',
  plot = grid.arrange(ppfithos,ppfitdead,ncol=1),
  width = 11,
  height = 11,
  units = "in",
  dpi = 300
)


print('1 of 2 finished ...')


ResultsList <- foreach(rr = regions, .combine='comb', .multicombine=TRUE,
                      .init= list(list(),list(),list(),list(),list(),list()) ) %dopar% {

library(data.table)
library(readxl)
library(WriteXLS)
library(R.matlab)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(tikzDevice)
library(latex2exp)

##############################################################
################################################ 005

#################################### EARLY LATE TESTING increasing capacity

    ###
    hospit <- dead <- NULL
    hospLate<-deadLate<-hospEarly<-deadEarly<-deadcap<-hoscap<-deadcapM<-hoscapM<-NULL
    for (ii in 1:pertnum) {
        
        ###LATE TESTING
        tmpplo <- readMat(paste0(thisishere,'/LateTest_res/',ii,'_',rr,'_cap_005.mat'))
        AAd<-tmpplo$x[,5]+tmpplo$x[,6]+tmpplo$x[,15]+tmpplo$x[,16]
        BBd<-tmpplo$x[,7]+tmpplo$x[,8]
        CCd<-tmpplo$x[,9]+tmpplo$x[,10]
        DDd<-tmpplo$x[,12]
        
        { hospLate <- cbind(hospLate,(BBd+CCd)) }
        
        { deadLate <- cbind(deadLate,DDd) }
        
        ###EARLY TESTING
        tmpplo <- readMat(paste0(thisishere,'/EarlyTest_res/',ii,'_',rr,'_cap_005.mat'))
        AAc<-tmpplo$x[,5]+tmpplo$x[,6]+tmpplo$x[,15]+tmpplo$x[,16]
        BBc<-tmpplo$x[,7]+tmpplo$x[,8]
        CCc<-tmpplo$x[,9]+tmpplo$x[,10]
        DDc<-tmpplo$x[,12]
        
        ll<-min(length(DDc),length(DDd))
        DD <- DDc-DDd#[1:ll] - DDd[1:ll]
        BB <- (BBc+CCc) - (BBd+CCd)
        
        { hospEarly <- cbind(hospEarly,(BBc+CCc)) }
        
        { deadEarly <- cbind(deadEarly,DDc) }
        
        ###NORMAL CAP
        tmpplo <- readMat(paste0(thisishere,'/cap_res/',ii,'_',rr,'.mat'))
#         tmpplo <- readMat(paste0(thisishere,'/deCap_res/',ii,'_',rr,'.mat'))
        
        HHcap <- tmpplo$x[,7]+tmpplo$x[,8]+tmpplo$x[,9]+tmpplo$x[,10]
        DDcap <- tmpplo$x[,12]
        
        deadcap <- cbind(deadcap, DDcap)
        hoscap <- cbind(hoscap, HHcap)
    }
    
    ######means and sd
    miss <- minlen-nrow(deadcap)
    print(miss)
    
    tt <- seq(miss,(minlen-1))
    rtmp <- rep(rr,length(tt) )
    

    ###Early Testing + MaxCap
    hosEA <- rowMeans(hospEarly); deadCP <- rowMeans(deadEarly)
    hosSC <- rowSds(hospEarly); deadSC <- rowSds(deadEarly)

    capmean <- cbind(tt, rtmp, rep('TestCap -- Early',length(hosEA)), hosEA )
    capsd <- cbind(tt, rtmp, rep('TestCap -- Early',length(hosEA)), hosSC )
    
    deadcapm <- cbind(tt, rtmp, rep('TestCap -- Early',length(hosEA)), deadCP )
    deadcapsd <- cbind(tt, rtmp, rep('TestCap -- Early',length(hosEA)), deadSC )
    ############################
    ###Late Testing + MaxCap
    hosT <- rowMeans(hospLate); deadT <- rowMeans(deadLate)
    hosST <- rowSds(hospLate); deadST <- rowSds(deadLate)

    testmean <- cbind(tt, rtmp, rep('TestCap -- Late',length(hosEA)), hosT )
    testsd <- cbind(tt, rtmp, rep('TestCap -- Late',length(hosEA)), hosST )
    
    deadtestm <- cbind(tt, rtmp, rep('TestCap -- Late',length(hosEA)), deadT )
    deadtestsd <- cbind(tt, rtmp, rep('TestCap -- Late',length(hosEA)), deadST )
    #######################
    ###reference cap
    deadcapM <- rowMeans(deadcap); deadcapsd <- rowSds(deadcap)
    hoscapM <- rowMeans(hoscap); hoscapsds <- rowSds(hoscap)
  
    caprefe <- cbind(tt, rtmp, rep('Capacity',length(hoscap)), hoscapM )
    caprefesd <- cbind(tt, rtmp, rep('Capacity',length(hoscap)), hoscapsds )
    
    deadcaprefe <- cbind(tt, rtmp, rep('Capacity',length(hosEA)), deadcapM )
    deadcaprefesd <- cbind(tt, rtmp, rep('Capacity',length(hosEA)), deadcapsd )
    
    
    #################
    HE <- round((max(hoscapM)-max(hosEA))/max(hoscapM)*100)#, digits=2)
    HL <- round((max(hoscapM)-max(hosT))/max(hoscapM)*100)#, digits=2)
    DE <- round((max(deadcapM)-max(deadCP))/max(deadcapM)*100)#, digits=2)
    DL <- round((max(deadcapM)-max(deadT))/max(deadcapM)*100)#, digits=2)
    
    keepD <- rbind( cbind( DE, rr, (max(deadCP)),'TestCap -- Early' ),cbind( DL, rr, (max(deadT)),'TestCap -- Late' )  )
    keepH <-rbind( cbind( HE, rr, (max(hosEA)) ,'TestCap -- Early' ),cbind( HL, rr, (max(hosT)) ,'TestCap -- Late' )  )
    collDeadCap <- data.frame(( keepD ))
    collHosCap <- data.frame(( keepH ))
    
    for (jj in c(1,3)) {
        collDeadCap[,jj] <- as.numeric(as.character(collDeadCap[,jj]))
        collHosCap[,jj] <- as.numeric(as.character(collHosCap[,jj]))
    }

    colnames(collDeadCap) <- colnames(collHosCap) <- c('value','Area','where','Model')

    ###########
    
    ###
    meancollect <- data.frame( rbind( cbind(capmean, deadCP), cbind(testmean, deadT), cbind(caprefe,deadcapM) ) )
    
    sdcollect <- data.frame( rbind( cbind(capsd, deadSC), cbind(testsd, deadST), cbind(caprefesd,deadcapsd) ) )
    ############
    
    ######Data
    dataReal <- read_excel(paste0(thisishere,'/../Data/',rr,".xls"))
    dataReal <- dataReal[1:length(tt),]
    hosdata <- dataReal[,which(colnames(dataReal)=='hos')]; iicdata <- dataReal[,which(colnames(dataReal)=='icu')]; dd <- dataReal[,which(colnames(dataReal)=='dead')];
    
    data_coll <- data.frame(cbind(tt,hosdata+iicdata,dd,rtmp))
    ####
    for (jj in c(1,4,5)) {
        meancollect[,jj] <- as.numeric(as.character(meancollect[,jj]))
        sdcollect[,jj] <- as.numeric(as.character(sdcollect[,jj]))
    }

    for (jj in 1:(ncol(data_coll)-1)) {
        data_coll[,jj] <- as.numeric(as.character(data_coll[,jj]))
    }

    colnames(meancollect) <- c('Time','Area','Model','Hospitalized','Dead')
    colnames(sdcollect) <- c('Time','Area','Model','Hospitalized','Dead')
    colnames(data_coll) <- c('Time','Hospitalized','Dead','Area')
    
    datamelt <- melt(data_coll, id.vars=c('Time','Area'))

    captest_collect <- melt(meancollect, id.vars=c('Time','Area','Model'))
    sdtest_collect <- melt(sdcollect, id.vars=c('Time','Area','Model'))
    captest_collect$sd <- sdtest_collect$value
    
    ###hosp capacities
    hh <- readMat(paste0(thisishere,'/rt_cap/1_',rr,'.mat'))
    hh <- hh$parsave[3:nrow(hh$parsave), 16:17]
#     hlim <- rbind( hlim, cbind(tt, hh, rep(rr,nrow(hh))) )
    hlim <- data.frame( cbind(tt[1:nrow(hh)], hh, rep(rr,nrow(hh))) ) 
    MINHH <- data.frame( cbind( min(hh[,1]+hh[,2]), rr ) )
    MAXHH <- data.frame( cbind( max(hh[,1]+hh[,2]), rr ) )
    colnames(hlim) <- c('Time','Hospitalized','ICU','Area') 
    colnames(MINHH) <- colnames(MAXHH) <- c('Limit','Area')

    MINHH$Limit <- as.numeric(as.character(MINHH$Limit))

    MAXHH$Limit <- as.numeric(as.character(MAXHH$Limit))
    
    for (jj in 1:(ncol(hlim)-1)) {
        hlim[,jj] <- as.numeric(as.character(hlim[,jj]))
    }
    hlim <- melt(hlim, id.vars=c('Time','Area'))
    
    MINHH <- melt(MINHH, id.vars=c('Area'))
    MAXHH <- melt(MAXHH, id.vars=c('Area'))
    ###
    list(captest_collect,datamelt,collHosCap,collDeadCap,MINHH,MAXHH)
}


captest_collect<-datamelt<-collHosCap<-collDeadCap<-MINHH<-MAXHH<-NULL
for (RL in 1:length(ResultsList[[1]])) {
    captest_collect <- rbind(captest_collect, as.data.frame(ResultsList[[1]][[RL]]))
    datamelt <- rbind(datamelt, as.data.frame(ResultsList[[2]][[RL]]))
    collHosCap <- rbind(collHosCap, as.data.frame(ResultsList[[3]][[RL]]))
    collDeadCap <- rbind(collDeadCap, as.data.frame(ResultsList[[4]][[RL]]))
    MINHH <- rbind(MINHH, as.data.frame(ResultsList[[5]][[RL]]))
    MAXHH <- rbind(MAXHH, as.data.frame(ResultsList[[6]][[RL]]))
}


levels(captest_collect$Area)[which(levels(captest_collect$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')

levels(datamelt$Area)[which(levels(datamelt$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')

levels(collHosCap$Area)[which(levels(collHosCap$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')
levels(collDeadCap$Area)[which(levels(collDeadCap$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')

levels(MINHH$Area)[which(levels(MINHH$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')

levels(MAXHH$Area)[which(levels(MAXHH$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna', 'Friuli Venezia Giulia', 'Valle d\'Aosta')


sel <- c('Emilia Romagna','Lombardia','Piemonte')
supplsel <- setdiff(regions,sel)

HH <-captest_collect[which(captest_collect$variable=='Hospitalized'),]
HH <- HH[which(HH$Area %in% sel),]
DD <- captest_collect[which(captest_collect$variable=='Dead'),]
DD <- DD[which(DD$Area %in% sel),]

dmsel <- datamelt[which(datamelt$Area %in% sel ),]
dmselh <- dmsel[which(dmsel$variable=='Hospitalized'),]
dmseld <- dmsel[which(dmsel$variable=='Dead'),]

MAXHHmain <- MAXHH[which(MAXHH$Area %in% sel),]
MAXHHsupp <- MAXHH[which(MAXHH$Area %in% supplsel),]
MINHHmain <- MINHH[which(MINHH$Area %in% sel),]

collselhos <- collHosCap[which(collHosCap$Area %in% HH$Area),]
collseldead <- collDeadCap[which(collDeadCap$Area %in% DD$Area),]


ppfithos <- ggplot(HH, aes(x=as.Date(Time, origin='2020-02-24',by='days',format='%Y-%m-%d'), y=value, col=Model)) +
geom_line() +
geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill=Model),col=NA, alpha=0.5) +
geom_point(dmselh, mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value), col='black', size=0.1) +
geom_hline(MINHHmain, mapping=aes(yintercept=value ) ) + #5,25,40,70
#geom_hline(MAXHHmain, mapping=aes(yintercept=value ) ) + #5,25,40,70
#geom_text(collselhos, mapping=aes(label = value, y=where, x=as.Date(90, origin='2020-02-24',by='days',format='%Y-%m-%d') ,col=Model)) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
#scale_y_sqrt(breaks = function(x) { pretty(seq(min(x),max(x),100)) } ) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Hospitalized') + ggtitle('A')

updown <- c(-100,0, 0,1000, -200, 0,0,0)
ppfitdead <- ggplot(DD, aes(x=as.Date(Time, origin='2020-02-24',by='days',format='%Y-%m-%d'), y=value, col=Model)) +
geom_line() +
geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill=Model),col=NA, alpha=0.5) +
geom_point(dmseld, mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value), col='black', size=0.1) +
geom_text(collseldead, mapping=aes(label = paste0(value,'%'), y=where,
#+updown, 
x=as.Date(95, origin='2020-02-24',by='days',format='%Y-%m-%d') ,col=Model), show.legend = FALSE) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
expand_limits(x = as.Date('2020-06-07')) +
#scale_y_sqrt(breaks = function(x) { pretty(seq(min(x),max(x),100)) } ) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ylab('Dead') + ggtitle('B')

#pdf('TestCap_005.pdf',h=11,w=11)
#grid.arrange(ppfithos,ppfitdead,ncol=1)
#dev.off()

ggsave(
  'TestCap_005.pdf',
  plot = grid.arrange(ppfithos,ppfitdead,ncol=1),
  width = 11,
  height = 11,
  units = "in",
  dpi = 300
)

print('2 of 2 finished.')

stopCluster(cl)
