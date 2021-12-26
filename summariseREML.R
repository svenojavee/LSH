library(data.table)
library(ggplot2)
library(ggsci)

palette <- wes_palette("Darjeeling1", 2)


classicFormula <- function(h2Obs, K, P) {
    return(h2Obs * (K * (1-K))**2 /(dnorm(qnorm(K)))**2 /(P*(1-P)) )
}

newFormula <- function(h2Obs,K){
  t <- qnorm(1-K)
  z <- dnorm(t)
  xx <- h2Obs
  a1 <- 0
  a2 <- 1
  it <- 1
  while(abs(a2-a1) > 0.00001){
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
    if(it > 100){
      print("Has not converged in 100 steps")
      return(xx)
    }
  }
  return(xx)
}


wd <- "/mnt/beegfs/robingrp/sojavee/LSHResults/results/"

h2_grid <- c( 0.15, 0.25, 0.35, 0.45, 0.55, 0.65)
P_grid <- c(0.5,0.75,0.9,0.95,1) # how many percent is P from K
K_grid <- c(0.05, 0.02, 0.01)
nSims <- 100

dat <- NULL

for(h2 in h2_grid){
    print(h2)
    for(K in K_grid){
        for(P in P_grid){
            P_use <- P * K
            for(n in 1:nSims){
                h2Obs <- fread(paste0(wd,"s",n,"_M10k_N20k_h",h2,"_K",K,"_P",P,".hsq"),fill=T)$Variance[4]
                h2Classic <- classicFormula(h2Obs, K, P_use)
                h2New <- newFormula(h2Obs,P_use)
                res_vec <- c(h2,K,P,n,h2Classic,h2New)
                dat <- rbind(dat,res_vec)
            }
        }
    }
}
dat <- as.data.table(dat)
names(dat) <- c("h2","K","P","n","h2Classic","h2New")

dat[,n:=NULL]
datLong <- melt(dat,id.vars=c("h2","K","P"))

datLong[,meanValue:=mean(value),by=list(h2,K,P,variable)]
datLong[,sdValue:=sd(value),by=list(h2,K,P,variable)]
datLong[,mse:=mean((value-h2)**2),by=list(h2,K,P,variable)]

datLongUnique <- unique(datLong[,list(h2,K,P,variable,meanValue,sdValue,mse)])

datLongUnique[,ylow:=meanValue-sdValue]
datLongUnique[,yhigh:=meanValue+sdValue]

datLongUnique[variable == "h2Classic",variable:="Previous expression"]
datLongUnique[variable == "h2New",variable:="New expression"]
datLongUnique[,P_perc:=P*100]
datLongUnique[,K_format:=paste0("K = ",K)]
datLongUnique[,h_format:=paste0("h2 = ",h2)]


p <- ggplot(data=datLongUnique,aes(x=h2,y=meanValue,col=variable,fill=variable)) + facet_wrap(~as.character(K)+as.character(P)) + 
    geom_line() + geom_ribbon(aes(ymin = ylow, ymax = yhigh), alpha = 0.1) + geom_abline(intercept=0,slope=1,lty=2)


p_mse <- ggplot(dat = datLongUnique, aes(x = P_perc, y = mse, color = variable)) + geom_line() + facet_grid(K_format ~ h_format)+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, 'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 14), legend.key.size = unit(0.5, 'cm'), legend.spacing.x = unit(0.2, 'cm')) +
    theme(title = element_text(color = "gray20", size = 14)) +
    theme(axis.title = element_text(color = "gray20", size = 14)) +
    theme(axis.text = element_text(color = "gray40", size = 10)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(size = 14, color = "gray10")) +
    xlab("Sample prevalence (% of population prevalence)") + ylab("MSE") + scale_x_continuous(trans='log10') +
    scale_color_lancet()


ggsave(filename="/nfs/scistore13/robingrp/sojavee/LSH/remlPlot_mse.pdf",plot=p_mse,width=10,height=5)

