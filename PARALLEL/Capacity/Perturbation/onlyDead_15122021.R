library(data.table)
library(readxl)
library(WriteXLS)
library(R.matlab)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(gridExtra)
library(latex2exp)
Sys.setlocale("LC_TIME", "C")

thisishere <- getwd()

alphaAA <- c('alpha02','alpha03','alpha04','alpha05')

setwd('dead_diff_alpha02')
#pertnum <- 101
fordd <- list.files(pattern = "\\.mat$")
regions <- NULL
for (ii in 1:length(fordd)) {
    ff <- fordd[ii]
    
    minlen<-nrow(readMat(ff)$x);print(ii);print(minlen)
}
regions <- 'Lombardia'#unique(regions)
setwd(thisishere)
pertnum <- length(fordd)/length(regions)
######Hlim Icum

alpha_coll <- NULL
for (aa in alphaAA) {
    setwd(thisishere)
    setwd(paste0('dead_diff_',aa))

    dead_coll <- dead_sd <- NULL
    dead_endpoint <- NULL
    for (rr in regions) {
        print(rr)
        deadTMP<-NULL
        for (ii in 1:length(fordd)) {
            ff <- fordd[ii]
            if (unlist(strsplit(unlist(strsplit(ff,'_'))[2],'[.]'))[1] == rr) 
            {
                tmpplo <- readMat(ff)
                deadTMP <- cbind(deadTMP,tmpplo$x)
            }
        }
        miss <- minlen-nrow(deadTMP)
        print(miss)
        
        A <- rowMeans(deadTMP)
        As <- rowSds(deadTMP)
        tt <- seq(miss,(minlen-1))
        rtmp <- rep(rr,length(A))
        dead_coll <- rbind(dead_coll,cbind(tt,A,rtmp))
        dead_sd <- rbind(dead_sd,cbind(tt,As,rtmp))
        
        dd <- c(deadTMP[nrow(deadTMP),],rr)
        dead_endpoint <- rbind(dead_endpoint,dd)
    }

    setwd('..')
    dead_coll <- data.frame(dead_coll)
    dead_sd <- data.frame(dead_sd)
    dead_endpoint <- data.frame(dead_endpoint)
    #error <- data.frame(error)

    colnames(dead_coll) <- colnames(dead_sd) <-
    c('Time','Dead Diff.','Area') 
    colnames(dead_endpoint) <- c(seq(1,pertnum),'Region')

    for (jj in 1:pertnum) {dead_endpoint[,jj] <- as.numeric(as.character(dead_endpoint[,jj]))}
    ddEnd <- melt(dead_endpoint,id.vars='Region')
    
    a_num <- as.numeric(as.numeric(unlist(strsplit(aa,''))[7]))

    dead_endpoint$Region <- paste('alpha=0.',a_num)
#     dead_endpoint$Region <- (a_num)
    
    alpha_coll <- rbind(alpha_coll, dead_endpoint)
    
}

alpha_coll <- data.frame(alpha_coll)
for (jj in 1:(ncol(alpha_coll)-1)) {
    alpha_coll[,jj] <- as.numeric(as.character(alpha_coll[,jj]))
}

AAc <- melt(alpha_coll,id.vars=c('Region'))

boxplExcD <- ggplot(AAc, aes(x=Region, y=value,)) + 
        geom_boxplot(#fill="slateblue", 
        alpha=0.2, outlier.shape=NA) + 
        scale_y_sqrt(breaks = c( 5,50, seq(100,500,100), seq(1000,7000,1000) ) ) + #5,25,40,70
    theme_bw() +
theme(panel.grid.minor = element_blank(),
strip.background = element_rect(color="white", fill="white"),
axis.text.x=element_text(angle=45, margin = margin(t = 5),hjust=1),
panel.background = element_rect(fill="white"), 
plot.title = element_text(face = "bold"),
text = element_text(size = 15)) +
scale_x_discrete(labels=c(expression(paste(alpha,' = ',0.2)),expression(paste(alpha,' = ',0.3)),expression(paste(alpha,' = ',0.4)),expression(paste(alpha,' = ',0.5)))) +
xlab('') + ylab('Excess Dead') + ggtitle('Lombardia')

boxplExcD <- boxplExcD #+ stat_compare_means()
#compare_means(value ~ Region,  data = AAc)


#pdf('ALPHA_Perturb.pdf',h=7,w=9)
#boxplExcD
#dev.off()



ggsave(
  'ALPHA_Perturb.pdf',
  plot = boxplExcD,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)

print('alpha perturbation Plot done.')



