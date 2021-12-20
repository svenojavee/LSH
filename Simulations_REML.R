library(data.table)
library(genio)


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


#M <- 10000
N <- 25000

tmp1 <- function(sigG,sigE){
  return(sigG/(sigG+sigE))
}


simulateFunction_SNP <- function(dataSet, h2, N, K, P) {
  NTilde <- nrow(dataSet)
  M <- ncol(dataSet)
    
    betEff <- rnorm(M,mean=0,sd=sqrt(h2/M))
    g <- X_mat %*% betEff

    threshold <- qnorm(1-K)

    varG <- var(g)
    e <- rnorm(NTilde,mean=0,sd=sqrt(1-varG))
    l <- g + e
    y <- rep(0, NTilde)
    y[l > threshold] <- 1

  # Select new sample with the defined probabilities
  p1 <- N * (1 - P) / NTilde / (1 - K) 
  p2 <- N * P / NTilde / K

  N1 <- N * (1 - P) 
  N2 <- N * P 

    Indices2 <- (1:length(y))[y == 1]
    Indices1 <- (1:length(y))[y == 0]

    final2 <- sample(Indices2,N2)
    final1 <- sample(Indices1,N1)

  y_final <- c(y[final1], y[final2])
  id_final <- c(final1, final2)
    l_final <- c(l[final1], l[final2])

  retDat <- data.table(id_final, id_final, y_final)
  names(retDat) <- c("IID", "FID", "phen")
  retDat <- retDat[order(IID),]
  #g_final <- c(g[final1], g[final2])

  return(retDat)
}


h2_grid <- c( 0.15, 0.25, 0.35, 0.45, 0.55, 0.65)
P_grid <- c(0.5,0.75) # how many percent is P from K
K_grid <- c(0.05, 0.02, 0.01)


#Read the genetic data

#d1 <- read_plink("/Users/admin/Documents/plotLSH/simREML")
d1 <- read_plink("/home/svenerik.ojavee/LSH_sim/simM10k_N25k")
#scale
X_mat <- scale(t(d1$X))


for (h2_val in h2_grid) {
    print(h2_val)
    for (K_val in K_grid) {
      for (P_val in P_grid) {
        P_use <- K_val * P_val

          #Generate the data
          samp_val <- simulateFunction_SNP(X_mat, h2_val, N, K_val, P_use)
          fwrite(samp_val,file="/home/svenerik.ojavee/LSH_sim/tmp.phen",sep=" ",col.names=F)
         system(paste0("/data/software/gcta_1.93.3beta2/gcta64 --grm simM10k_N75k --keep tmp.phen --pheno tmp.phen --reml --threads 20 --out /home/svenerik.ojavee/LSH_sim/resDir/res_",h2_val,"_",K_val,"_",P_val))
      }
    }
}

