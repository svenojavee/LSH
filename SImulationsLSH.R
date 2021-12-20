library(data.table)
library(MASS)
library(doParallel)
library(doSNOW)


originalFormula <- function(h2Obs, K) {
    return(h2Obs * K * (1-K) /(dnorm(qnorm(K)))**2)
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


#1. Simulate weakly related individuals

N <- 10000
blockSize <- 100
relatedness <- 0.05
noReplicates <- 250

simulateFunction <- function(h2, N, blockSize, relatedness, K) {
    threshold <- qnorm(1-K)
    g <- c()
  itCount <- N / blockSize
  covMat <- diag(blockSize) * h2 + ((diag(blockSize) * 0 + 1) - diag(blockSize)) * h2 * relatedness
    mu <- rep(0,blockSize)
    for (i in 1:itCount) {
      g_samp <- MASS::mvrnorm(1, mu = mu, Sigma = covMat)
          g <- c(g, g_samp)
    }
  varG <- var(g)
    e <- rnorm(N,mean=0,sd=sqrt(1-varG))
  l <- g + e
  y <- rep(0, N)
  y[l > threshold] <- 1

  return(list(y,g))
}

h2_grid <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55,0.65,0.75)
K_grid <- c(0.001,0.005,0.01,0.05,0.1,0.5)
K_imprecision <- c(0.5,0.2,0)
resultsDat <- NULL

for (h2_val in h2_grid) {
  for (K_val in K_grid) {
    for (imp_val in K_imprecision) {

      for (n in 1:noReplicates) {
        print(n)
        samp_val <- simulateFunction(h2_val, N, blockSize, relatedness, K_val)
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

            h2Original_vec[jj] <- originalFormula(h2Obs, K_use)
            h2New_vec[jj] <- newFormula(h2Obs, K_use)
          }

          h2Original <- mean(h2Original_vec)
          h2New <- mean(h2New_vec)
        } else {
          h2Original <- originalFormula(h2Obs, K_val)
          h2New <- newFormula(h2Obs, K_val)
        }

        res_vec <- c(h2_val, K_val,imp_val, h2Original, h2New)
        resultsDat <- rbind(resultsDat, res_vec)
      }
    }
    }
}



h2_grid <- c( 0.15, 0.25, 0.35, 0.45, 0.55, 0.65)
#h2_grid <- c(0.05)

K_grid <- c(0.001,0.005,0.01,0.05,0.1,0.5)
K_imprecision <- c(0.2,0.1,0)


cl <- makeCluster(64, outfile = "") 

registerDoSNOW(cl)

#Write the same version in parallel for the replicates
res <- foreach(n = 1:noReplicates, .combine = 'rbind') %dopar% {
resultsDat <- NULL
  for (h2_val in h2_grid) {
    print(h2_val)
    for (K_val in K_grid) {
    for (imp_val in K_imprecision) {
        samp_val <- simulateFunction(h2_val, N, blockSize, relatedness, K_val)
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
            while (K_use <= 0.0001) {
              K_use <- max(rnorm(1, mean = K_val, sd = imp_val * K_val), 0)
            }

            h2Original_vec[jj] <- originalFormula(h2Obs, K_use)
            h2New_vec[jj] <- newFormula(h2Obs, K_use)
          }

          h2Original <- mean(h2Original_vec)
          h2New <- mean(h2New_vec)
        } else {
          h2Original <- originalFormula(h2Obs, K_val)
          h2New <- newFormula(h2Obs, K_val)
        }

        res_vec <- c(h2_val, K_val,imp_val, h2Original, h2New)
        resultsDat <- rbind(resultsDat, res_vec)
      
    }
    }
}
resultsDat

    
}

stopCluster(cl)



resFinal <- as.data.table(res)
names(resFinal) <- c("h2", "K","Imprecision", "origH2", "newH2")

resultsDatLong <- melt(resFinal,id.vars=c("h2","K","Imprecision"))

resultsDatLong[,mse:=mean((value-h2)**2),by=c("h2","K","Imprecision","variable")]
resultsDatLong[,varH2:=var(value),by=c("h2","K","Imprecision","variable")]
resultsDatLong[,bias2:=(mean(value)-h2)**2,by=c("h2","K","Imprecision","variable")]

mseDat <- unique(resultsDatLong[,list(h2,  K, Imprecision ,variable,mse,varH2,bias2)])

fwrite(mseDat,file="/home/svenerik.ojavee/mseDat.txt")



summary(resultsDatLong$value[resultsDatLong$variable=="origH2"])
summary(resultsDatLong$value[resultsDatLong$variable=="newH2"])
