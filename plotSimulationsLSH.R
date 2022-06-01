library(data.table)
library(ggplot2)
library(ggsci)


mseDatLong <- fread("/mnt/beegfs/robingrp/sojavee/LSHResults/imprecisionSimulation.txt")
h2Grid <- seq(0.15,0.65,0.1)
KGrid <- c(0.001,0.005,0.01,0.05,0.1,0.5)
ImprecisionGrid <- c(0.0,0.1)
variableGrid <- c("origH2","newH2")
#Number of BS replicates
B <- 250

#Do the bootstrapping procedure
bsDat <- NULL
for(h2Use in h2Grid){
    print(h2Use)
    for(KUse in KGrid){
        print(KUse)
        for(ImprecisionUse in ImprecisionGrid){
            for(variableUse in variableGrid){
                #Take the data set 
                h2Use <- round(h2Use,2)  #Sth very strange from R
                subDat <- mseDatLong[mseDat$h2 == h2Use & mseDat$K == KUse & mseDat$Imprecision == ImprecisionUse & mseDat$variable == variableUse,]
                for(b in 1:B){
                    reDat <- sample(subDat$value,size=2000,replace=T)
                    mseVal <- mean((reDat - h2Use)**2)
                    recVec <- c(h2Use,KUse,ImprecisionUse,variableUse,mseVal)
                    bsDat <- rbind(bsDat,recVec)
                }
            }
        }
    }
}
bsDat <- as.data.table(bsDat)
names(bsDat) <- c("h2","K","Imprecision","variable","mse")
bsDat[,mse:=as.numeric(mse)]
bsDat[,h2:=as.numeric(h2)]
bsDat[,K:=as.numeric(K)]
bsDat[,Imprecision:=as.numeric(Imprecision)]

bsDat[,q025:=quantile(mse,probs=0.025),by=list(h2,K,Imprecision,variable)]
bsDat[,q975:=quantile(mse,probs=0.975),by=list(h2,K,Imprecision,variable)]
bsDat <- unique(bsDat[,list(h2,K,Imprecision,variable,q025,q975)])

mseDatLong[,mse:=mean((value-h2)**2),by=list(h2,K,Imprecision,variable)]
mseDat <- unique(mseDatLong[,list(h2,K,Imprecision,variable,mse)])

mseDat <- merge(mseDat,bsDat,by=c("h2","K","Imprecision","variable"))

#Make a more joint plot
mseDat[variable == "origH2",variable:="Previous expression"]
mseDat[variable == "newH2", variable := "New expression"]
mseDat[,variable:=factor(variable,levels=c("Previous expression","New expression"))]
mseDat[,h_format:=paste0("h2 = ",h2)]
mseDat[,Imprecision_format:=paste0("K ~ N(K,",Imprecision," x K)")]
mseDat[,K_format:=paste0("K = ",K)]


p <- ggplot(dat = mseDat[Imprecision == 0,], aes(x = K, y = mse)) + geom_line(aes(color=variable)) + geom_ribbon(aes(ymin=q025,ymax=q975,fill=variable),alpha=0.1) +
    facet_wrap(. ~ h_format, scales = "free",ncol=3)+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, 'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12), legend.key.size = unit(0.5, 'cm'), legend.spacing.x = unit(0.2, 'cm')) +
    theme(title = element_text(color = "gray20", size = 13)) +
    theme(axis.title = element_text(color = "gray20", size = 14)) +
    theme(axis.text = element_text(color = "gray40", size = 10)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(size = 14, color = "gray10")) +
    xlab("Prevalence") + ylab("MSE") + scale_x_continuous(trans='log10') + scale_color_lancet() + scale_fill_lancet()


ggsave(p,filename="/nfs/scistore13/robingrp/sojavee/figures/LSH_integrate_N20k.pdf",width=8,height=5)



#Summarise MSE increase K~N()
mseSub <- mseDatLong[,list(h2,K,Imprecision,variable,value)]
mseSub0 <- mseSub[Imprecision == 0,list(h2,K,variable,value)]
mseSub1 <- mseSub[Imprecision == 0.1,list(h2,K,variable,value)]
mseSub0[,ind:=1:.N,by=list(h2,K,variable)]
mseSub1[,ind:=1:.N,by=list(h2,K,variable)]

mseSub_wide <- merge(mseSub0,mseSub1,by=c("h2","K","variable","ind"))
setnames(mseSub_wide,c("value.x","value.y"),c("0","0.1"))



B <- 250

#Do the bootstrapping procedure
bsDat <- NULL
for(h2Use in h2Grid){
    print(h2Use)
    for(KUse in KGrid){
        print(KUse)
            for(variableUse in variableGrid){
                #Take the data set 
                h2Use <- round(h2Use,2)  #Sth very strange from R
                subDat <- mseSub_wide[h2 == h2Use & K == KUse & variable == variableUse,]
                for(b in 1:B){
                    reDat <- sample(subDat$ind,size=2000,replace=T)
                    mseVal0 <- mean((subDat$`0`[reDat] - h2Use)**2)
                    mseVal1 <- mean((subDat$`0.1`[reDat] - h2Use)**2)
                    mseIncrease <- 100 * (mseVal1-mseVal0)/mseVal0
                    recVec <- c(h2Use,KUse,variableUse,mseIncrease)
                    bsDat <- rbind(bsDat,recVec)
                }
            }
    }
}
bsDat <- as.data.table(bsDat)
names(bsDat) <- c("h2","K","variable","mseIncrease")
bsDat[,mseIncrease:=as.numeric(mseIncrease)]
bsDat[,h2:=as.numeric(h2)]
bsDat[,K:=as.numeric(K)]

bsDat[,q025:=quantile(mseIncrease,probs=0.025),by=list(h2,K,variable)]
bsDat[,q975:=quantile(mseIncrease,probs=0.975),by=list(h2,K,variable)]
bsDat <- unique(bsDat[,list(h2,K,variable,q025,q975)])

mseSub_wide[,mse0:=mean((`0`-h2)**2),by=list(h2,K,variable)]
mseSub_wide[,mse1:=mean((`0.1`-h2)**2),by=list(h2,K,variable)]

mseSubShort <- unique(mseSub_wide[,list(h2,K,variable,mse0,mse1)])
mseSubShort[,mseDiff:=100 * (mse1-mse0)/mse0,by=list(h2,K,variable)]

mseSubShort <- merge(mseSubShort,bsDat,by=c("h2","K","variable"))
mseSubShort[,h_format:=paste0("h2 = ",h2)]
mseSubShort[,K_format:=paste0("K = ",K)]
mseSubShort[variable == "origH2",variable:="Previous expression"]
mseSubShort[variable == "newH2", variable := "New expression"]

p <- ggplot(dat = mseSubShort[h2 >= 0.15 & h2 <= 0.65 & K < 0.1], aes(x = variable, y = mseDiff, fill = variable,col=variable)) +
     geom_bar(stat = "identity",alpha=0.5) +
     geom_errorbar(aes(ymin=q025, ymax=q975), width=.3, lwd=0.3, alpha=0.5) + facet_grid(K_format ~ h_format,scales="free_y") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, 'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12), legend.key.size = unit(0.5, 'cm'), legend.spacing.x = unit(0.2, 'cm')) +
    theme(title = element_text(color = "gray20", size = 13)) +
    theme(axis.title = element_text(color = "gray20", size = 14)) +
        theme(axis.text = element_text(color = "gray40", size = 10)) +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(size = 14, color = "gray10")) +
    xlab("") + ylab("MSE increase (%)") + scale_fill_lancet() + scale_color_lancet() 

ggsave(p,filename="/nfs/scistore13/robingrp/sojavee/figures/LSH_integrate_mseIncrease.pdf",width=10,height=6.5)

