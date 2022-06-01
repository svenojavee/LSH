library(data.table)

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
  while(abs(a2-a1) > 0.000001){
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


LSH_ascertain <- function(h2Obs,P,K){
  tP <- qnorm(P)
  tK <- qnorm(K)
  z <- dnorm(tK)
  xx <- h2Obs
  a1 <- 0
  a2 <- 1
  it <- 1
  while(abs(a2-a1) > 0.00000001){
    it <- it+1
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



h2Dat <- fread("/nfs/scistore13/robingrp/sojavee/LSH/h2_results.csv")
prevDat <- fread("/nfs/scistore13/robingrp/sojavee/LSH/Prevalences.csv")
prevDat[,PopulationPrevalence:=PopulationPrevalence/100]
prevDat[,SamplePrevalence:=SamplePrevalence/100]

dat <- merge(prevDat,h2Dat,by.x="Code",by.y="disease")

datLong <- melt(dat,id.vars=c("Code","Disease","SamplePrevalence","PopulationPrevalence","Notes"))
datLong[,LSH_classicFull:=classicFormula(value,PopulationPrevalence,SamplePrevalence)]
datLong[,LSH_classicCompare:=classicFormula(value,SamplePrevalence,SamplePrevalence)]

datLong[,LSH_new:=0]
for(i in 1:nrow(datLong)){
    if(datLong$PopulationPrevalence[i] == datLong$SamplePrevalence[i]){
      print(datLong$Disease[i])
      datLong$LSH_new[i] <- newFormula(datLong$value[i],datLong$PopulationPrevalence[i])
    }else{
      datLong$LSH_new[i] <- LSH_ascertain(datLong$value[i],datLong$SamplePrevalence[i],datLong$PopulationPrevalence[i])
    }
}


newDatLong <- datLong[,list(Code,Disease,variable,LSH_new)]
newDatLong[,LSH_new:=round(LSH_new,3)]
newDat <- dcast(newDatLong,Code+Disease~variable)
newDat[,valPresent:=paste0(h2," (",lCI,", ",uCI,")")]

classicDatLong <- datLong[,list(Code,Disease,variable,LSH_classicCompare)]
classicDatLong[,LSH_classicCompare:=round(LSH_classicCompare,3)]
classicDat <- dcast(classicDatLong,Code+Disease~variable)
classicDat[,valPresent:=paste0(h2," (",lCI,", ",uCI,")")]


classic2DatLong <- datLong[,list(Code,Disease,variable,LSH_classicFull)]
classic2DatLong[,LSH_classicFull:=round(LSH_classicFull,3)]
classic2Dat <- dcast(classic2DatLong,Code+Disease~variable)
classic2Dat[,valPresent:=paste0(h2," (",lCI,", ",uCI,")")]


plotData <- merge(classic2Dat[,list(Code,Disease,h2,lCI,uCI)],newDat[,list(Code,Disease,h2,lCI,uCI)],by=c("Code","Disease"))
plotData <- plotData[!(Code %in% c("F20","F42","F33")),]
library(ggplot2)
require("ggrepel")

plotData[,list(Code,h2.x,lCI.x,uCI.x,h2.y,lCI.y,uCI.y)]
#plotData[,prMode:=paste0(h2.y," (",lCI.y,", ",uCI.y,")")]
plotData[,ciClassic:=uCI.x-lCI.x]
plotData[,ciNew:=uCI.y-lCI.y]
plotData[,CIreduction:=round(100*(ciClassic-ciNew)/ciClassic)]

p <- ggplot(data=plotData,aes(x=h2.x,y=h2.y,label=Disease))  + geom_text_repel()  +
    geom_errorbarh(aes(xmax = uCI.x, xmin = lCI.x),alpha=0.3) + geom_errorbar(aes(ymax = uCI.y, ymin = lCI.y),alpha=0.3)  + geom_point(col=2)+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, 'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12), legend.key.size = unit(0.5, 'cm'), legend.spacing.x = unit(0.2, 'cm')) +
    theme(title = element_text(color = "gray20", size = 10)) +
    theme(axis.title = element_text(color = "gray20", size = 12)) +
    theme(axis.text = element_text(color = "gray40", size = 10)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(size = 8, color = "gray10")) +
    xlab("Previous estimate") + ylab("New estimate")  +  geom_abline(intercept = 0, slope = 1,lty=3)

ggsave(p,filename="/nfs/scistore13/robingrp/sojavee/figures/RealDataCompare.pdf",width=6,height=5)
