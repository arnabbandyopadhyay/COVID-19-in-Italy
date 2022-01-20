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

#collect_mean <- collect_sd <- data_coll <- NULL

# ResultsList <- foreach(i=1:10, .combine='comb', .multicombine=TRUE,
#                 .init=list(list(), list())) %dopar% {
#   #list(c_m,c_sd,d_c)
# }

#list(simulDynamic, data, Rtdynamic)
ResultsList <- foreach(rr = regions, .combine='comb', .multicombine=TRUE,
                      .init= list(list(),list(),list()) ) %dopar% {
                      setwd(RefeFolder)
library(data.table)
library(readxl)
library(WriteXLS)
library(R.matlab)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(tikzDevice)
library(latex2exp)

#      print(rr)
    quar <- hosp <- icu <- dead <- NULL
    rtM <- rtASY <- NULL
    for (ii in 1:length(fordd)) {
        ff <- fordd[ii]
        if (unlist(strsplit(unlist(strsplit(ff,'_'))[2],'[.]'))[1] == rr) {
            tmpplo <- readMat(ff)
            tmpploRT <- readMat(paste0('../Refer_rt_res/',ff))
            tmpploRTASY <- readMat(paste0('../Asym_rt_res/',ff))
            
            rtM <- cbind(rtM,tmpploRT$t) #Rt
            rtASY <- cbind(rtASY,tmpploRTASY$t)
            
            quar <- cbind(quar,(tmpplo$x[,5]+tmpplo$x[,6]))
            hosp <- cbind(hosp,(tmpplo$x[,7]+tmpplo$x[,8]))
            icu <- cbind(icu,(tmpplo$x[,9]+tmpplo$x[,10]))
            dead <- cbind(dead,tmpplo$x[,12])
        }
    }
    miss <- minlen-nrow(quar)
    missPAR <- minlen-nrow(rtM)-5 
    
    
    ####Rt
    meanRt <- rowMeans(rtM); meanRtAsy <- rowMeans(rtASY)
    
    sdRt <- rowSds(rtM); sdRtAsy <- rowSds(rtASY) #Bs <- rowSds(hosp);
    ttPAR <- seq(missPAR,(nrow(rtM)+missPAR-1))#seq(miss,(minlen-1))
    rtmp <- rep(rr,length(meanRt))
    coll_rtmean <- cbind(ttPAR,meanRt,meanRtAsy,rtmp)
    coll_rtsd <- cbind(ttPAR,sdRt,sdRtAsy,rtmp)
    ####
    
    ###inf.dyn
    A <- rowMeans(quar); B <- rowMeans(hosp); C <- rowMeans(icu); D <- rowMeans(dead)
    As <- rowSds(quar); Bs <- rowSds(hosp); Cs <- rowSds(icu); Ds <- rowSds(dead)
    tt <- seq(miss,(minlen-1))
    rtmp <- rep(rr,length(A))
#     collect_mean <- rbind(collect_mean,cbind(tt,A,B,C,D,rtmp))
#     collect_sd <- rbind(collect_sd,cbind(tt,As,Bs,Cs,Ds,rtmp))
    collect_mean <- cbind(tt,A,B,C,D,rtmp)
    collect_sd <- cbind(tt,As,Bs,Cs,Ds,rtmp)
    ###
    ###Data
    setwd('../../Data')
    dataReal <- read_excel(paste0(rr,".xls"))
    dataReal <- dataReal[1:length(tt),]
    qq <- dataReal[,which(colnames(dataReal)=='qua')]; hh <- dataReal[,which(colnames(dataReal)=='hos')]; iic <- dataReal[,which(colnames(dataReal)=='icu')]; dd <- dataReal[,which(colnames(dataReal)=='dead')];
    
    data_coll <- cbind(tt,qq,hh,iic,dd,rtmp)
    
    
    collect_mean <- data.frame(collect_mean)
    collect_sd <- data.frame(collect_sd)
    data_coll <- data.frame(data_coll)


    colnames(collect_mean) <- colnames(collect_sd) <- colnames(data_coll) <-c('Time','Active Infection','Hospitalized','ICU','Dead','Area')

    for (jj in 1:(ncol(data_coll)-1)) {
        data_coll[,jj] <- as.numeric(as.character(data_coll[,jj]))
        collect_mean[,jj] <- as.numeric(as.character(collect_mean[,jj]))
        collect_sd[,jj] <- as.numeric(as.character(collect_sd[,jj])) 
    }

    ccm <- melt(collect_mean,id.vars=c('Time','Area'))
    ccs <- melt(collect_sd,id.vars=c('Time','Area'))
    ddc <- melt(data_coll,id.vars=c('Time','Area'))
#eee <- melt(error,id.vars=c('Time','Area'))

#eee$value <- as.numeric(as.character(eee$value))

    ccm$sd <- ccs$value

    ccm$value <- as.numeric(as.character(ccm$value))
    ccm$sd <- as.numeric(as.character(ccm$sd))

    names(ccm)[3] <- names(ddc)[3] <- 'Legend'
    names(ccm)[4] <- names(ddc)[4] <- 'Cases'

    levels(ccm$Time) <- levels(ddc$Time) <- as.Date( as.numeric(as.character(levels(ccm$Time))), origin='2020-02-24',by='days')

    recordit <- ccm
    recdata <- ddc
    
    #####rt
    
    coll_rtmean <- data.frame(coll_rtmean)
    coll_rtsd <- data.frame(coll_rtsd)

    colnames(coll_rtmean) <- colnames(coll_rtsd) <- c('Time','Reference','Asymptomatic','Area')
    
    for (jj in 1:3 ) {
        coll_rtmean[,jj] <- as.numeric(as.character(coll_rtmean[,jj]))
        coll_rtsd[,jj] <- as.numeric(as.character(coll_rtsd[,jj]))
    }
    ##


    rtmm <- melt(coll_rtmean,id.vars=c('Time','Area'))
    rtss <- melt(coll_rtsd,id.vars=c('Time','Area'))

    rtmm$sd <- rtss$value

    rtmm$value <- as.numeric(as.character(rtmm$value))
    rtmm$sd <- as.numeric(as.character(rtmm$sd))

    names(rtmm)[3] <- 'Model'

    rtmm$ci <- (rtmm$sd*1.96)/sqrt(pertnum) 


    ###

     setwd(thisishere)
     list(recordit,recdata,rtmm)
}

simulResul<-dataResul<-RtResul<-NULL
for (RL in 1:length(ResultsList[[1]])) {
    simulResul <- rbind(simulResul, as.data.frame(ResultsList[[1]][[RL]]))
    dataResul <- rbind(dataResul, as.data.frame(ResultsList[[2]][[RL]]))
    RtResul <- rbind(RtResul, as.data.frame(ResultsList[[3]][[RL]]))
}


levels(simulResul$Area)[which(levels(simulResul$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')
# 
levels(dataResul$Area)[which(levels(dataResul$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')

levels(RtResul$Area)[which(levels(RtResul$Area)%in%c('Emilia','Friuli','Valledaosta'))] <- c('Emilia Romagna','Friuli Venezia Giulia','Valle d\'Aosta')

# 
simulResul$Area <- relevel(simulResul$Area, "Italy")
RtResul$Area <- relevel(RtResul$Area, "Italy")
# 
sel <- c('Italy', 'Lombardia', 'Marche', 'Veneto', 'Piemonte','Emilia Romagna')
supplsel <- setdiff(regions,sel)
# 
mainReg <- simulResul[which(simulResul$Area %in% sel),]
suppplo <- simulResul[which(simulResul$Area %in% supplsel),]
mainData <- dataResul[which(dataResul$Area %in% sel),]
suppData <- dataResul[which(dataResul$Area %in% supplsel),]

mainRt <- RtResul[which(RtResul$Area %in% sel),]
supplRt <- RtResul[which(RtResul$Area %in% supplsel),]

ppfit <- ggplot(mainReg, aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=Cases, col=Legend)) +
geom_line() +
geom_ribbon(aes(ymin=Cases-sd, ymax=Cases+sd, fill=Legend),col=NA, alpha=0.5) +
geom_point(mainData, mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=Cases, col=Legend), size=0.5) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
scale_y_sqrt(breaks = function(x) { pretty(seq(min(x),max(x),100)) } ) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
xlab('') + ggtitle('A')

ppfitsupp <- ggplot(suppplo, aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=Cases, col=Legend)) +
geom_line() +
geom_ribbon(aes(ymin=Cases-sd, ymax=Cases+sd, fill=Legend),col=NA, alpha=0.5) +
geom_point(suppData, mapping=aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=Cases, col=Legend), size=0.5) +
scale_x_date(date_breaks = "15 days",date_labels = "%b %e") +
scale_y_sqrt(breaks = function(x) { pretty(seq(min(x),max(x),100)) } ) +
facet_wrap(~Area, scale='free_y')+#,ncol=6) +
theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 25)) +
xlab('')

pp <- ggplot(mainRt, aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Model)) +
geom_line() +
geom_ribbon(aes(ymin=value-ci, ymax=value+ci, fill=Model),col=NA, alpha=0.5) +
#geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill=Model),col=NA, alpha=0.5) +
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
text = element_text(size = 15)
) +
xlab('') + 
# ylab((TeX('Reproduction number ($R_t$)'))) +
ylab( expression( paste('Reproduction Number ', italic(R)[t]) ) ) +
# expression(paste('Values of ', mu))
# ylab(TeX('Reproduction number $\\mathcal{R}_{t}$')) +
ggtitle('B')

ppsupp <- ggplot(supplRt, aes(x=as.Date(Time, origin='2020-02-24',by='days'), y=value, col=Model)) +
geom_line() +
geom_ribbon(aes(ymin=value-ci, ymax=value+ci, fill=Model),col=NA, alpha=0.5) +
#geom_ribbon(aes(ymin=value-sd, ymax=value+sd, fill=Model),col=NA, alpha=0.5) +
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
text = element_text(size = 25)
) +
xlab('') + 
# ylab((TeX('Reproduction number ($R_t$)'))) +
ylab( expression( paste('Reproduction Number ', italic(R)[t]) ) )


stopCluster(cl)

ggsave(
 'RefAsMod_Fit_Rt_Main.pdf',
  plot = grid.arrange(ppfit,pp,ncol=1),
  width = 11,
  height = 11,
  units = "in",
  dpi = 300
)

print('Main Fit and Rt Plot done ...')

ggsave(
  'RefAsMod_Fit_Suppl.pdf',
  plot = ppfitsupp,
  width = 22,
  height = 13,
  units = "in",
  dpi = 300
)

print('Supplementary Fit Plot done ...')

ggsave(
  'RefAsMod_Rt_Suppl.pdf',
  plot = ppsupp,
  width = 22,
  height = 8,
  units = "in",
  dpi = 300
)

print('Supplementary Rt Plot done.')
