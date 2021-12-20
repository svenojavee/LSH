library(data.table)
library(MASS)
library(doParallel)
library(doSNOW)


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




#1. Simulate weakly related individuals such that P < K
N <- 10000
blockSize <- 100
relatedness <- 0.05
noReplicates <- 250

simulateFunctionPK <- function(h2, N, blockSize, relatedness, K, P) {
    NTilde <- 2*round(N * (1-P)/(1-K), -log10(blockSize)) + blockSize

    threshold <- qnorm(1-K)
    g <- c()
    itCount <- NTilde / blockSize
    covMat <- diag(blockSize) * h2 + ((diag(blockSize) * 0 + 1) - diag(blockSize)) * h2 * relatedness
    mu <- rep(0,blockSize)
    for (i in 1:itCount) {
      g_samp <- MASS::mvrnorm(1, mu = mu, Sigma = covMat)
      g <- c(g, g_samp)
    }
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
  g_final <- c(g[final1], g[final2])

  return(list(y_final,g_final))
}


h2_grid <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55,0.65,0.75)
K_grid <- c(0.2, 0.1, 0.01)
P_grid <- c(0.5,0.25) # how many percent is P
K_imprecision <- c(0.5,0.2,0)
resultsDat <- NULL

for (h2_val in h2_grid) {
  for (K_val in K_grid) {
    for (P_val in P_grid) {
        P_use <- K_val * P_val
        for (imp_val in K_imprecision) {

      for (n in 1:noReplicates) {
        print(n)
        samp_val <- simulateFunctionPK(h2_val, N, blockSize, relatedness, K_val, P_use)
        g_use <- samp_val[[2]]
        y_use <- samp_val[[1]]
        m1 <- lm(y_use ~ g_use)
        h2Obs <- summary(m1)$adj.r.squared

        if (imp_val != 0) {
          intSize <- 100
          h2Original_vec <- rep(NA, intSize)
          h2New_vec <- rep(NA, intSize)

          for (jj in 1:intSize) {
            K_use <- 0
            while (K_use == 0) {
              K_use <- max(rnorm(1, mean = K_val, sd = 0.5 * K_val), 0)
            }

            h2Original_vec[jj] <- classicFormula(h2Obs, K_use, P_use)
            h2New_vec[jj] <- newFormula(h2Obs, K_use)
          }

          h2Original <- mean(h2Original_vec)
          h2New <- mean(h2New_vec)
        } else {
          h2Original <- classicFormula(h2Obs, K_val, P_use)
          h2New <- newFormula(h2Obs, P_use)
        }

        res_vec <- c(h2_val, K_val,P_val,imp_val, h2Original, h2New)
        resultsDat <- rbind(resultsDat, res_vec)
      }
    }

          
      }
    
    }
}






h2_grid <- c( 0.15, 0.25, 0.35, 0.45, 0.55, 0.65)
P_grid <- c(0.5,0.25) # how many percent is P from K
K_grid <- c(0.2, 0.1, 0.01)

cl <- makeCluster(64, outfile = "") 

registerDoSNOW(cl)

#Write the same version in parallel for the replicates
res <- foreach(n = 1:noReplicates, .combine = 'rbind') %dopar% {
resultsDat <- NULL
  for (h2_val in h2_grid) {
    print(h2_val)
    for (K_val in K_grid) {
      for (P_val in P_grid) {
        P_use <- K_val * P_val

          samp_val <- simulateFunctionPK(h2_val, N, blockSize, relatedness, K_val, P_use)
          g_use <- samp_val[[2]]
          y_use <- samp_val[[1]]
          m1 <- lm(y_use ~ g_use)
          h2Obs <- summary(m1)$adj.r.squared

          h2Original <- classicFormula(h2Obs, K_val, P_use)
          h2New <- newFormula(h2Obs, P_use)

          res_vec <- c(h2_val, K_val, P_val, h2Original, h2New)
          resultsDat <- rbind(resultsDat, res_vec)
        
      }
    }
}
resultsDat

    
}

stopCluster(cl)



resFinal <- as.data.table(res)
names(resFinal) <- c("h2", "K","PComp", "origH2", "newH2")

resultsDatLong <- melt(resFinal,id.vars=c("h2","K","PComp"))

resultsDatLong[,mse:=mean((value-h2)**2),by=c("h2","K","PComp","variable")]
resultsDatLong[,varH2:=var(value),by=c("h2","K","PComp","variable")]
resultsDatLong[,bias2:=(mean(value)-h2)**2,by=c("h2","K","PComp","variable")]

mseDat <- unique(resultsDatLong[,list(h2,  K,PComp ,variable,mse,varH2,bias2)])

fwrite(mseDat,file="/home/svenerik.ojavee/mseDat_PK.txt")



