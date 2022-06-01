#We define the functions that can be used to calculate LSH estimate, SE of LSH estimate and two types of CI
#We need the library pbivnorm


###################################
#Functions for ascertained versions

LSH_ascertain <- function(h2Obs,P,K, precision = 0.0000001, convergenceSteps = 100){
  tP <- qnorm(P)
  tK <- qnorm(K)
  z <- dnorm(tK)
  xx <- h2Obs
  a1 <- 0
  a2 <- 1
  it <- 1
  while(abs(a2-a1) > precision){
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
    if(it > convergenceSteps){
      print(paste0("Has not converged in ",convergenceSteps," steps"))
      return(xx)
    }
  }

  return(xx)
}

LSH_ascertain_SE <- function(h2Obs,P,K,SEh2Obs){
  require(pbivnorm)
  corDer <- function(t,h2Liab){
    return(exp(-t**2/(1 + h2Liab))/(2*pi*sqrt(1 - h2Liab ** 2)))
  }
  derFun <- function(h2Liab,K,P){
    minT <- qnorm(K)
    zz2 <- dnorm(minT) ** 2
    numerat <- 2 * h2Liab * corDer(minT,h2Liab) - zz2 * h2Liab + K - pbivnorm(minT,minT,h2Liab)
    denomin <- (zz2 * h2Liab + K - pbivnorm(minT,minT,h2Liab))**3
    return(zz2 * P * (1-P) * numerat / denomin)
  }
  h2Liab <- LSH_ascertain(h2Obs,P,K) #That is the existing function to calculate LSH
  return(SEh2Obs / derFun(h2Liab,K,P))
}

#Convenience function to output the classical confidence interval 
CI_classical_ascertain <- function(h2Obs,P,K,SEh2Obs,alphaThreshold=0.05){
    LSHEst <- LSH_ascertain(h2Obs,P,K)
    LSHSEEst <- LSH_ascertain_SE(h2Obs,P,K,SEh2Obs)
    print(paste0("(",LSHEst - LSHSEEst*qnorm(1-alphaThreshold/2),",",LSHEst + LSHSEEst*qnorm(1-alphaThreshold/2),")"))
}
#Convenience function to output a non-symmetric CI that is guaranteed to be bounded between 0 and 1
#(using clog-log transform)
CI_nonSymmetric_cloglog_ascertain <- function(h2Obs,P,K,SEh2Obs,alphaThreshold=0.05){
    LSHEst <- LSH_ascertain(h2Obs,P,K)
    LSHSEEst <- LSH_ascertain_SE(h2Obs,P,K,SEh2Obs)
    Q <- LSHSEEst / abs(LSHEst * log(LSHEst))
    z <- qnorm(1-alphaThreshold/2)
    print(paste0("(", exp(log(LSHEst) * exp(Q * z) ) ,",",exp(log(LSHEst) * exp(-Q * z) ),")"))
}



######################################
#Functions for unascertained versions


LSH <- function(h2Obs,K, precision = 0.0000001, convergenceSteps = 100){
    t <- qnorm(1-K)
    z2 <- dnorm(t) **2
    xx <- h2Obs
    a1 <- 0
    a2 <- 1
    it <- 1
    while(abs(a2-a1) > precision){
        it <- it+1
        res <- (xx * z2)/(xx * z2 + 1-K - pbivnorm::pbivnorm(t,t,xx)) - h2Obs
        if(res < 0){ #Find the solution from between xx and a2
            a1 <- xx
            xx <- 0.5*(a2 + xx)
        }
        if(res > 0){
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

LSH_SE <- function(h2Obs,K,SEh2Obs){
  require(pbivnorm)
  corDer <- function(t,h2Liab){
    return(exp(-t**2/(1 + h2Liab))/(2*pi*sqrt(1 - h2Liab ** 2)))
  }
  derFun <- function(h2Liab,K){
    minT <- qnorm(K)
    zz2 <- dnorm(minT) ** 2
    numerat <- corDer(minT,h2Liab) * h2Liab * zz2 + zz2 * (K - corDer(minT,h2Liab))
    denomin <- (h2Liab * zz2 + K - pbivnorm(minT,minT,h2Liab))**2
    return(numerat/denomin)
  }
    h2Liab <- LSH(h2Obs,K) #That is the existing function to calculate LSH
    return(SEh2Obs / derFun(h2Liab,K))
}

#Convenience function to output the classical confidence interval 
CI_classical <- function(h2Obs,K,SEh2Obs,alphaThreshold=0.05){
    LSHEst <- LSH(h2Obs,K)
    LSHSEEst <- LSH_SE(h2Obs,K,SEh2Obs)
    print(paste0("(",LSHEst - LSHSEEst*qnorm(1-alphaThreshold/2),",",LSHEst + LSHSEEst*qnorm(1-alphaThreshold/2),")"))
}
#Convenience function to output a non-symmetric CI that is guaranteed to be bounded between 0 and 1
#(using clog-log transform)
CI_nonSymmetric_cloglog <- function(h2Obs,K,SEh2Obs,alphaThreshold=0.05){
    LSHEst <- LSH(h2Obs,K)
    LSHSEEst <- LSH_SE(h2Obs,K,SEh2Obs)
    Q <- LSHSEEst / abs(LSHEst * log(LSHEst))
    z <- qnorm(1-alphaThreshold/2)
    print(paste0("(", exp(log(LSHEst) * exp(Q * z) ) ,",",exp(log(LSHEst) * exp(-Q * z) ),")"))
}




