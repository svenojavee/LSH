library(data.table)
library(genio)



M <- 10000
NTilde <- 26000 #Generate 22500 individuals

    mafs <- runif(M, 0.05, 0.5)
    X_vec <- rbinom(NTilde * M, 2, prob = rep(mafs, NTilde))
    X_unscaled <- (matrix(X_vec,nrow=M))

write_plink(file="/home/svenerik.ojavee/LSH_sim/simM10k_N75k",X_unscaled)


