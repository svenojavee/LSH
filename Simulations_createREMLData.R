library(data.table)
library(genio)



M <- 10000
NTilde <- 20700 #Generate (20540 should be enough for K=0.05 and P=0.025 - worst case sim scenario)

mafs <- runif(M, 0.05, 0.5)
X_vec <- rbinom(NTilde * M, 2, prob = rep(mafs, NTilde))
X_unscaled <- (matrix(X_vec,nrow=M))

write_plink(file="/mnt/beegfs/robingrp/sojavee/LSHResults/simM10k_N20k",X_unscaled)


