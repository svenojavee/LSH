library(data.table)
library(ggplot2)
library(ggsci)
library(ggpubr)


classicFormula <- function(h2Obs, K, P) {
    return(h2Obs * (K * (1-K))**2 /(dnorm(qnorm(K)))**2 /(P*(1-P)) )
}

newFormula <- function(h2Obs,K, precision = 0.00001, convergenceSteps = 100){
  t <- qnorm(1-K)
  z <- dnorm(t)
  xx <- h2Obs
  a1 <- 0
  a2 <- 1
  it <- 1
  while(abs(a2-a1) > precision){
    it <- it+1
    res <- (xx*z**2)/(xx*z**2 + 1-K - pbivnorm::pbivnorm(t,t,xx))  - h2Obs
    if(res < 0){ #Find the solution from between xx and a2
      a1 <- xx
      xx <- 0.5*(a2 + xx)
    }
    if(res>0){
      a2 <- xx
      xx <- 0.5*(a1 + xx)
    }
    if(it > convergenceSteps){
      print(paste0("Has not converged in ",convergenceSteps," steps"))
      return(xx)
    }
  }
  return(xx)
}

newFormulaPK <- function(h2Obs,P,K){
  tP <- qnorm(P)
  tK <- qnorm(K)
  z <- dnorm(tK)
  xx <- h2Obs
  a1 <- 0
  a2 <- 1
  it <- 1
  while(abs(a2-a1) > 0.00001){
    it <- it+1
    #res <- h2Obs * (K * (1-K))**2 / z**2 * (1/(P - pbivnorm::pbivnorm(tP,tP,xx))) - xx
    res <- h2Obs * (xx * z**2 + K - pbivnorm::pbivnorm(tK,tK,xx))**2 / z**2 * (1/(P*(1-P))) - xx
 
    if(res > 0){ #Find the solution from between xx and a2
      a1 <- xx
      xx <- 0.5*(a2 + xx)
    }
    if(res < 0){
      a2 <- xx
      xx <- 0.5*(a1 + xx)
    }
    if(it > 100){
      print("Has not converged in 100 steps")
      return(xx)
    }
  }
  return(xx)
}


wd <- "/mnt/beegfs/robingrp/sojavee/LSHResults/results/"

h2_grid <- c( 0.15, 0.25, 0.35, 0.45, 0.55, 0.65)
#P_grid <- c(0.5,0.75,0.9,0.95,1,1.25,1.5,2,4) # how many percent is P from K
P_grid <- c(0.25,0.5,0.75,0.9,1,1.25,1.5,2,4)

K_grid <- c(0.05, 0.02, 0.01, 0.005)
nSims <- 100

dat <- NULL
for(h2 in h2_grid){
    print(h2)
    for(K in K_grid){
        for(P in P_grid){
            P_use <- P * K
            for(n in 1:nSims){
              if(T | ( !( h2 %in% c(0.25) & K %in% c(0.005) & P %in% c(1.25) & n %in% c(17)) ) & 
                ( !( h2 %in% c(0.35) & K %in% c(0.005) & P %in% c(1.25) & n %in% c(16)) ) & 
                ( !( h2 %in% c(0.35) & K %in% c(0.005) & P %in% c(1.5) & n %in% c(46)) ) & 
                ( !( h2 %in% c(0.35) & K %in% c(0.005) & P %in% c(2) & n %in% c(62)) )){
                h2Obs <- fread(paste0(wd,"s",n,"_M10k_N20k_h",h2,"_K",K,"_P",P,".hsq"),fill=T)$Variance[4]
                h2Classic <- classicFormula(h2Obs, K, P_use)
                h2New <- newFormula(h2Obs,P_use)
                h2NewAdjust <- newFormulaPK(h2Obs,P_use,K)
                res_vec <- c(h2,K,P,n,h2Classic,h2New,h2NewAdjust)
                dat <- rbind(dat,res_vec)
              }
                
            }
        }
    }
}
dat <- as.data.table(dat)
names(dat) <- c("h2","K","P","n","h2Classic","h2New","h2NewAdjust")

dat[,n:=NULL]
datLong <- melt(dat,id.vars=c("h2","K","P"))

#Read the PCGC results
datPCGC <- fread("/mnt/beegfs/robingrp/sojavee/LSHResults/sumStatsPCGC/allPCGCSims.txt")
datPCGC[,n:=NULL]
datPCGC[,value:=h2Est]
datPCGC[,variable:="PCGC"]
datPCGC <- datPCGC[,list(h2,K,P,variable,value)]

datLong <- rbind(datLong,datPCGC)
#Let's take each case and bootstrap this
methodVec <- c("h2Classic","h2New","h2NewAdjust","PCGC")
B <- 250 # Number of bootstrap replicates
bsResults <- NULL
for(h2Use in h2_grid){
    print(h2Use)
    for(KUse in K_grid){
      print(KUse)
        for(PUse in P_grid){
            for(method in methodVec){
              subDat <- datLong[h2 == h2Use & KUse == K & PUse == P & method == variable,]
              for(b in 1:B){
                bsDat <- sample(subDat$value,size=100,replace=T)
                mse <- mean((bsDat - h2Use)**2)
                bias <- mean(bsDat) - h2Use
                recVec <- c(h2Use,KUse,PUse,method,mse,bias)
                bsResults <- rbind(bsResults,recVec)
              }
            }
        }
    }
    gc()
}
bsResults <- as.data.table(bsResults)
names(bsResults) <- c("h2","K","P","variable","mse","bias")
fwrite(bsResults,file="/mnt/beegfs/robingrp/sojavee/LSHResults/bsREML_PCGC.txt")

#Find the 95%CI
bsResults[,mse:=as.numeric(mse)]
bsResults[,bias:=as.numeric(bias)]
bsResults[,h2:=as.numeric(h2)]
bsResults[,P:=as.numeric(P)]
bsResults[,K:=as.numeric(K)]

bsResults[,MSEq975:=quantile(mse,probs=c(0.975)),by=list(h2,K,P,variable)]
bsResults[,MSEq025:=quantile(mse,probs=0.025),by=list(h2,K,P,variable)]
bsResults[,BIASq975:=quantile(bias,probs=c(0.975)),by=list(h2,K,P,variable)]
bsResults[,BIASq025:=quantile(bias,probs=0.025),by=list(h2,K,P,variable)]


bsResultsUnique <- unique(bsResults[,list(h2,K,P,variable,MSEq975,MSEq025,BIASq975,BIASq025)])

datLong[,meanValue:=mean(value),by=list(h2,K,P,variable)]
datLong[,sdValue:=sd(value),by=list(h2,K,P,variable)]
datLong[,mse:=mean((value-h2)**2),by=list(h2,K,P,variable)]

datLongUnique <- unique(datLong[,list(h2,K,P,variable,meanValue,sdValue,mse)])
datLongUnique <- merge(datLongUnique,bsResultsUnique,by=c("h2","K","P","variable"))


datLongUnique[,ylow:=meanValue-sdValue]
datLongUnique[,yhigh:=meanValue+sdValue]

datLongUnique[variable == "h2Classic",variable:="GREML+previous\nexpression"]
datLongUnique[variable == "h2New",variable:="GREML+new\nexpression"]
datLongUnique[variable == "h2NewAdjust",variable:="GREML+new adjusted\nexpression"]
datLongUnique[,P_perc:=P*100]
datLongUnique[,P_perc:=(P)]

datLongUnique[,K_format:=paste0("K = ",K)]
datLongUnique[,h_format:=paste0("h2 = ",h2)]


datLongUnique[,rmse:=sqrt(mse)]
datLongUnique[,bias:=meanValue-h2]

#Save the data for quick execution
fwrite(datLongUnique,file="/mnt/beegfs/robingrp/sojavee/LSHResults/summary_REML_PCGC.txt")

p_mse <- ggplot(dat = datLongUnique, aes(x = P_perc, y = mse, color = variable)) + geom_line() + 
    facet_grid(K_format ~ h_format,scales="free_y")+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, 'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 14), legend.key.size = unit(0.5, 'cm'), legend.spacing.x = unit(0.2, 'cm')) +
    theme(title = element_text(color = "gray20", size = 14)) +
    theme(axis.title = element_text(color = "gray20", size = 14)) +
    theme(axis.text = element_text(color = "gray40", size = 10)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(size = 11, color = "gray10")) +
    xlab("Sample prevalence (of population prevalence)") + ylab("MSE") + scale_x_continuous(trans='log2',breaks=c(0.25,0.5,1,2,4),labels=c("0.25","0.5","1","2","4")) +
    scale_color_lancet()

p_bias <- ggplot(dat = datLongUnique, aes(x = P_perc, y = bias, color = variable)) + geom_line() +
    facet_grid(K_format ~ h_format,scales="free_y")+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, 'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 14), legend.key.size = unit(0.5, 'cm'), legend.spacing.x = unit(0.2, 'cm')) +
    theme(title = element_text(color = "gray20", size = 14)) +
    theme(axis.title = element_text(color = "gray20", size = 14)) +
    theme(axis.text = element_text(color = "gray40", size = 10)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(size = 11, color = "gray10")) +
    xlab("Sample prevalence (of population prevalence)") + ylab("Bias") + scale_x_continuous(trans='log2',breaks=c(0.25,0.5,1,2,4),labels=c("0.25","0.5","1","2","4")) +
    scale_color_lancet()

p <- ggarrange(p_mse,p_bias,labels=c("A","B"),nrow=2,common.legend=T)


ggsave(filename="/nfs/scistore13/robingrp/sojavee/figures/remlPlot_mse.pdf",plot=p_mse,width=10,height=5)
ggsave(filename="/nfs/scistore13/robingrp/sojavee/figures/remlPlot_bias.pdf",plot=p_bias,width=10,height=5)
ggsave(filename="/nfs/scistore13/robingrp/sojavee/figures/remlPlot.pdf",plot=p,width=10,height=10)

