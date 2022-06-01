library(data.table)


wd <- "/mnt/beegfs/robingrp/sojavee/LSHResults/sumStatsPCGC"

P_vec <- c(0.25,0.5,0.75,0.9,1,1.25,1.5,2,4)
K_vec <- c(0.005,0.01,0.02,0.05)
h2_vec <- c(0.15,0.25,0.35,0.45,0.55,0.65)

sims <- 1:100

datAll <- NULL
for(h2 in h2_vec){
    print(h2)
    for(K in K_vec){
        for(P in P_vec){
            print(P)
            for(n in sims){
                tmp <- fread(paste0(wd,"/s",n,"_M10k_N20k_h",h2,"_K",K,"_P",P,".s",n,"_M10k_N20k_h",h2,"_K",K,"_P",P,".output"))
                h2Est <- tmp$Value[1]
                tmpVec <- c(h2,K,P,n,h2Est)
                datAll <- rbind(datAll,tmpVec)
            }
        }
    }
}

datAll <- as.data.table(datAll)
names(datAll) <- c("h2","K","P","n","h2Est")
fwrite(datAll,file="/mnt/beegfs/robingrp/sojavee/LSHResults/sumStatsPCGC/allPCGCSims.txt")

