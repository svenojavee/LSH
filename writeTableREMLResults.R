library(data.table)
library(WriteXLS)


dat <- fread("/mnt/beegfs/robingrp/sojavee/LSHResults/summary_REML_PCGC.txt")


#Keep only the important columns
dat <- dat[,list(h2,K,P,variable,mse,MSEq025,MSEq975,bias,BIASq025,BIASq975)]

dat[,mseWrite:=paste0(round(mse,4)," (",round(MSEq025,4),", ",round(MSEq975,4),")")]
dat[,biasWrite:=paste0(round(bias,3)," (",round(BIASq025,3),", ",round(BIASq975,3),")")]

datWideBias <- dcast(dat[,list(h2,K,P,variable,mseWrite,biasWrite)],h2+K+P~variable,value.var="biasWrite")
datWideMSE <- dcast(dat[,list(h2,K,P,variable,mseWrite,biasWrite)],h2+K+P~variable,value.var="mseWrite")

names(datWideBias) <- c("h2","K","P (times K)","GREML+new expression","GREML+new adjusted expression","GREML+previous expression","PCGC")
names(datWideMSE) <- c("h2","K","P (times K)","GREML+new expression","GREML+new adjusted expression","GREML+previous expression","PCGC")


WriteXLS(c("datWideMSE","datWideBias"),ExcelFileName="/mnt/beegfs/robingrp/sojavee/LSHResults/supplementaryTables/Simulation_GREML_PCGC.xlsx",SheetNames=c("MSE (95% CI)","Bias (95% CI)"),BoldHeaderRow=T,AdjWidth=T)