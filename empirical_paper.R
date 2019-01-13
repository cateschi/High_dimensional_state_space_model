setwd("C:/Users/Utente/Dropbox/Relevant Literature PhD Project/Application_google_trends")


# Settings
avg <- T # if TRUE, then aggregation is done on the aggregated Google trends by taking the average
fixed <- F # if TRUE, then the standard deviation of the RGB is kept as fixed
load.w <- T # if TRUE, matrix of factor loadings is estimated on a weekly basis
CCGT <- F # if TRUE, model includes auxiliary series of both CC and GT
fixed_sigma_u <- T
PLS <- F # if TRUE, uses the loadings from the PLS in the state space model
EN <- T # if TRUE, Elastic Net to target the predictors 
desGT <- F # if TRUE, deseasonalize Google Trends
lagGT <- F # if TRUE, include the 12 lags of the GTs as aditional variables
seas_fact <- F # if TRUE, estimate the seasonality of the factor in the state space model

monthlyGT <- F # if TRUE, use directly the monthly GT in the estimation from 2004 until 2017
until2012 <- F
until2013 <- F
until2014 <- F
until2015 <- F
until2016 <- F
until2017 <- F # if TRUE, download data from 01(02)/2010 until the end of 2017, otherwise from 01(02)/2004 until the end of 2011
f2007u2017 <- F
f2004u2017 <- T



#### Packages ####

library(timsac)
library(magic)
library(readxl)
library(data.table)
library(foreach)
library(glmnet)
library(devtools)
install_github("gabrielrvsc/HDeconometrics")
library(HDeconometrics)
library(xts)
library(fUnitRoots)
library(lubridate)
library(zoo)
library(cointReg)
library(qpcR)
library(openxlsx)
library(MVN)
library(plsdepot)
library(pls)
library(fpp)
library(mvnTest)
library(portes)
library(extRemes)
library(matrixcalc)



#### Functions ####

#Functions for Kalman filter
if (fixed==T){
  KF_slopes_univ <- function(par,y,opti,k,delta,outofsample,parP10,nstates){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- log(1155)
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8))
    Q <- D%*%R%*%D
    
    
    #Bulid T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tmatrix <- Ty
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){
      
      #Bulid Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Z <- Zy
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z)
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  }
  
  KF_slopes <- function(par,y,k,delta,opti,outofsample,parP10,nstates){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- log(1203)
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- par[9]
    sigma_omegax <- par[10]
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), 0, exp(sigma_Rx), exp(sigma_omegax)*diag(11))
    R[32,2] <- tanh(par[11])
    R[2,32] <- tanh(par[11])
    Q <- D%*%R%*%D
    H <- adiag(diag(0,5,5), exp(2*par[12]))
    
    #Bulid T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- adiag(Tymu, Tyomega)
    Tmatrix <- adiag(Ty, Tx)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){ 
      
      #Bulid Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- c(1,0)
      Zx <- rep(Zx,6)
      Zx <- c(Zx,1)
      Zx <- rbind(Zx)
      Z <- adiag(Zy,Zx)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + H
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,(ncol(waves)+1)] <- rep(0,nstates)
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  }
} else {
  KF_slopes_univ <- function(par,y,opti,k,delta,outofsample,parP10,nstates){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8))
    Q <- D%*%R%*%D

    
    #Bulid T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tmatrix <- Ty
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){
      
      #Bulid Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Z <- Zy
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z)
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl <- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  }
  
  KF_slopes <- function(par,y,k,delta,opti,outofsample,parP10,nstates){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- par[9]
    sigma_omegax <- par[10]
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    st.for <- matrix(0,nrow(y),len)
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), 0, exp(sigma_Rx), exp(sigma_omegax)*diag(11))
    R[32,2] <- tanh(par[11])
    R[2,32] <- tanh(par[11])
    Q <- D%*%R%*%D
    H <- adiag(diag(0,5,5), exp(2*par[12]))
    
    #Bulid T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- adiag(Tymu, Tyomega)
    Tmatrix <- adiag(Ty, Tx)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){ 
      
      #Bulid Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- c(1,0)
      Zx <- rep(Zx,6)
      Zx <- c(Zx,1)
      Zx <- rbind(Zx)
      Z <- adiag(Zy,Zx)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + H
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      for (j in 1:nrow(y)){
        st.for[j,i] <- epshatoutofsample[j]/sqrt(Fmatrix[j,j])
      }
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,(ncol(waves)+1)] <- rep(0,nstates)
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt,st.for=st.for))
    }
  }
}

if (fixed_sigma_u==T){
  KF_slopes_mixed_factor <- function(par,y,opti,k,delta,outofsample,parP10,nstates,lambda,H){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- log(1)
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), exp(sigma_Rx))
    gamma <- par[9]
    R[31,2] <- tanh(gamma)
    R[2,31] <- tanh(gamma)
    Q <- D%*%R%*%D
    
    #Build T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- as.matrix(1)
    Tmatrix <- adiag(Ty, Tx)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){
      
      #Bulild Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- as.matrix(lambda,length(lambda),1)
      Z <- adiag(Zy,Zx)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + adiag(diag(0,ncol(waves),ncol(waves)),diag(diag(H),nrow(y)-ncol(waves),nrow(y)-ncol(waves)))
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,c((ncol(waves)+1):ncol(Kg))] <- matrix(0,nstates,ncol(Kg)-(ncol(waves)))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  } #y1 has smooth trend and y2 has local level
  
  KF_slopes_mixed_factor_CC <- function(par,y,opti,k,delta,outofsample,parP10,nstates,lambda,H){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- par[9]
    sigma_omegax <- par[10]
    sigma_Rz <- log(1)
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    st.for <- matrix(0,nrow(y),len)
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), 0, exp(sigma_Rx), exp(sigma_omegax)*diag(11), exp(sigma_Rz))
    R[32,2] <- tanh(par[11])
    R[2,32] <- tanh(par[11])
    R[44,2] <- tanh(par[12])
    R[2,44] <- tanh(par[12])
    Q <- D%*%R%*%D
    
    #Bulid T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- adiag(Tymu, Tyomega)
    Tz <- as.matrix(1)
    Tmatrix <- adiag(Ty, Tx, Tz)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){ 
      
      #Bulid Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- c(1,0)
      Zx <- rep(Zx,6)
      Zx <- c(Zx,1)
      Zx <- rbind(Zx)
      Zz <- as.matrix(lambda,length(lambda),1)
      Z <- adiag(Zy,Zx,Zz)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + adiag(diag(0,ncol(waves),ncol(waves)),exp(2*par[13]),diag(diag(H),nrow(y)-ncol(waves)-1,nrow(y)-ncol(waves)-1))
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      for (j in 1:nrow(y)){
        st.for[j,i] <- epshatoutofsample[j]/sqrt(Fmatrix[j,j])
      }
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,(ncol(waves)+1)] <- rep(0,nstates)
        }
        if (is.na(epshatoutofsample[ncol(waves)+2,])){
          Kg[,c((ncol(waves)+2):ncol(Kg))] <- matrix(0,nstates,ncol(Kg)-(ncol(waves)+1))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
      
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt,st.for=st.for))
    }
  } #y1 and y2 have smooth trend and y3 has local level
} else {
  KF_slopes_mixed_factor <- function(par,y,opti,k,delta,outofsample,parP10,nstates,lambda,H){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- par[9]
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), exp(sigma_Rx))
    gamma <- par[10]
    R[31,2] <- tanh(gamma)
    R[2,31] <- tanh(gamma)
    Q <- D%*%R%*%D
    
    #Build T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- as.matrix(1)
    Tmatrix <- adiag(Ty, Tx)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){
      
      #Bulild Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- as.matrix(lambda,length(lambda),1)
      Z <- adiag(Zy,Zx)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + adiag(diag(0,ncol(waves),ncol(waves)),diag(diag(H),nrow(y)-ncol(waves),nrow(y)-ncol(waves)))
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,c((ncol(waves)+1):ncol(Kg))] <- matrix(0,nstates,ncol(Kg)-(ncol(waves)))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  } #y1 has smooth trend and y2 has local level
  
  KF_slopes_mixed_factor_CC <- function(par,y,opti,k,delta,outofsample,parP10,nstates,lambda,H){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- par[9]
    sigma_omegax <- par[10]
    sigma_Rz <- par[11]
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), 0, exp(sigma_Rx), exp(sigma_omegax)*diag(11), exp(sigma_Rz))
    R[32,2] <- tanh(par[12])
    R[2,32] <- tanh(par[12])
    R[44,2] <- tanh(par[13])
    R[2,44] <- tanh(par[13])
    Q <- D%*%R%*%D
    
    #Bulid T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- adiag(Tymu, Tyomega)
    Tz <- as.matrix(1)
    Tmatrix <- adiag(Ty, Tx, Tz)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){ 
      
      #Bulid Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- c(1,0)
      Zx <- rep(Zx,6)
      Zx <- c(Zx,1)
      Zx <- rbind(Zx)
      Zz <- as.matrix(lambda,length(lambda),1)
      Z <- adiag(Zy,Zx,Zz)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + adiag(diag(0,ncol(waves),ncol(waves)),exp(2*par[14]),diag(diag(H),nrow(y)-ncol(waves)-1,nrow(y)-ncol(waves)-1))
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,(ncol(waves)+1)] <- rep(0,nstates)
        }
        if (is.na(epshatoutofsample[ncol(waves)+2,])){
          Kg[,c((ncol(waves)+2):ncol(Kg))] <- matrix(0,nstates,ncol(Kg)-(ncol(waves)+1))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  } #y1 and y2 have smooth trend and y3 has local level
}

if (seas_fact==T && fixed_sigma_u==T){
  KF_slopes_mixed_factor_seas <- function(par,y,opti,k,delta,outofsample,parP10,nstates,lambda,H){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- log(1)
    sigma_omegax <- par[9]
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), exp(sigma_Rx), exp(sigma_omegax)*diag(11))
    gamma <- par[10]
    R[31,2] <- tanh(gamma)
    R[2,31] <- tanh(gamma)
    Q <- D%*%R%*%D
    
    #Build T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- adiag(as.matrix(1),Tyomega)
    Tmatrix <- adiag(Ty, Tx)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){
      
      
      #Bulild Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- as.matrix(lambda,length(lambda),1)
      Zx <- cbind(Zx,Zx%*%c(rep(c(1,0),5),1))
      Z <- adiag(Zy,Zx)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + adiag(diag(0,ncol(waves),ncol(waves)),diag(diag(H),nrow(y)-ncol(waves),nrow(y)-ncol(waves)))
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,c((ncol(waves)+1):ncol(Kg))] <- matrix(0,nstates,ncol(Kg)-(ncol(waves)))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  } #y1 has smooth trend and y2 has local level

  KF_slopes_mixed_factor_seas_CC <- function(par,y,opti,k,delta,outofsample,parP10,nstates,lambda,H){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- par[9]
    sigma_omegax <- par[10]
    sigma_Rz <- log(1)
    sigma_omegaz <- par[11]
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), 0, exp(sigma_Rx), exp(sigma_omegax)*diag(11), exp(sigma_Rz), exp(sigma_omegaz)*diag(11))
    R[32,2] <- tanh(par[12])
    R[2,32] <- tanh(par[12])
    R[44,2] <- tanh(par[13])
    R[2,44] <- tanh(par[13])
    Q <- D%*%R%*%D
    
    #Bulid T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- adiag(Tymu, Tyomega)
    Tz <- adiag(as.matrix(1),Tyomega)
    Tmatrix <- adiag(Ty, Tx, Tz)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){ 
      
      #Bulid Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- c(1,0)
      Zx <- rep(Zx,6)
      Zx <- c(Zx,1)
      Zx <- rbind(Zx)
      Zz <- as.matrix(lambda,length(lambda),1)
      Zz <- cbind(Zz,Zz%*%c(rep(c(1,0),5),1))
      Z <- adiag(Zy,Zx,Zz)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + adiag(diag(0,ncol(waves),ncol(waves)),exp(2*par[14]),diag(diag(H),nrow(y)-ncol(waves)-1,nrow(y)-ncol(waves)-1))
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,(ncol(waves)+1)] <- rep(0,nstates)
        }
        if (is.na(epshatoutofsample[ncol(waves)+2,])){
          Kg[,c((ncol(waves)+2):ncol(Kg))] <- matrix(0,nstates,ncol(Kg)-(ncol(waves)+1))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  } #y1 and y2 have smooth trend and y3 has local level
  
} else if (seas_fact==T && fixed_sigma_u==F){
  KF_slopes_mixed_factor_seas <- function(par,y,opti,k,delta,outofsample,parP10,nstates,lambda,H){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- par[9]
    sigma_omegax <- par[10]
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), exp(sigma_Rx), exp(sigma_omegax)*diag(11))
    gamma <- par[11]
    R[31,2] <- tanh(gamma)
    R[2,31] <- tanh(gamma)
    Q <- D%*%R%*%D
    
    #Build T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- adiag(as.matrix(1),Tyomega)
    Tmatrix <- adiag(Ty, Tx)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){
      
      #Bulild Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- as.matrix(lambda,length(lambda),1)
      Zx <- cbind(Zx,Zx%*%c(rep(c(1,0),5),1))
      Z <- adiag(Zy,Zx)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + adiag(diag(0,ncol(waves),ncol(waves)),diag(diag(H),nrow(y)-ncol(waves),nrow(y)-ncol(waves)))
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,c((ncol(waves)+1):ncol(Kg))] <- matrix(0,nstates,ncol(Kg)-(ncol(waves)))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  } #y1 has smooth trend and y2 has local level
  
  KF_slopes_mixed_factor_seas_CC <- function(par,y,opti,k,delta,outofsample,parP10,nstates,lambda,H){
    len <- length(y[1,])
    sigma_Ry <- par[1]
    sigma_omegay <- par[2]
    sigma_lambda <- par[3]
    sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
    sigma_Rx <- par[9]
    sigma_omegax <- par[10]
    sigma_Rz <- par[11]
    sigma_omegaz <- par[12]
    x10 <- rep(0,nstates)
    Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
    Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
    P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
    Pttm1[[1]] <- P10
    xtt <- matrix(0,nstates,(len))
    xttm1 <- matrix(0,nstates,(len+1))
    xttm1[,1] <- x10
    R <- diag(1,nstates,nstates)
    D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), 0, exp(sigma_Rx), exp(sigma_omegax)*diag(11), exp(sigma_Rz), exp(sigma_omegaz)*diag(11))
    R[32,2] <- tanh(par[13])
    R[2,32] <- tanh(par[13])
    R[44,2] <- tanh(par[14])
    R[2,44] <- tanh(par[14])
    Q <- D%*%R%*%D
    
    #Bulid T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- adiag(Tymu, Tyomega)
    Tz <- adiag(as.matrix(1),Tyomega)
    Tmatrix <- adiag(Ty, Tx, Tz)
    
    #initialization of loglikelihood
    logl <- 0
    
    #Start of KF recursions
    for (i in 1:len){ 
      
      #Bulid Z:
      Zy <- c(1,0)
      Zy <- rep(Zy,6)
      Zy <- c(Zy,1)
      Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
      Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
      Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
      Zx <- c(1,0)
      Zx <- rep(Zx,6)
      Zx <- c(Zx,1)
      Zx <- rbind(Zx)
      Zz <- as.matrix(lambda,length(lambda),1)
      Zz <- cbind(Zz,Zz%*%c(rep(c(1,0),5),1))
      Z <- adiag(Zy,Zx,Zz)
      ncol(Z)
      nrow(Z)
      
      epshatoutofsample <- y[,i] - Z%*%xttm1[,i]
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + adiag(diag(0,ncol(waves),ncol(waves)),exp(2*par[15]),diag(diag(H),nrow(y)-ncol(waves)-1,nrow(y)-ncol(waves)-1))
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if ((NaN %in% Fmatrix)==T){
        logl<- -P10[1]
      } else {
        svdFmatrix <- svd(Fmatrix)
        Kg <- Pttm1[[i]]%*%t(Z)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u) #kalman gain
        if (is.na(epshatoutofsample[1,])){
          Kg[,c(1:ncol(waves))] <- matrix(0,nstates,ncol(waves))
        }
        if (is.na(epshatoutofsample[ncol(waves)+1,])){
          Kg[,(ncol(waves)+1)] <- rep(0,nstates)
        }
        if (is.na(epshatoutofsample[ncol(waves)+2,])){
          Kg[,c((ncol(waves)+2):ncol(Kg))] <- matrix(0,nstates,ncol(Kg)-(ncol(waves)+1))
        }
        epshatoutofsample <- ifelse(is.na(epshatoutofsample), 0, epshatoutofsample)
        xtt[,i] <- xttm1[,i]+Kg%*%epshatoutofsample #compute x_{t|t}
        epshatinsample <- y[,i]-Z%*%xtt[,i] #in-sample forecast error (after y_t has been observed)
        epshatinsample <- ifelse(is.na(epshatinsample), 0, epshatinsample)
        Ptt[[i]] <- Pttm1[[i]]-Kg%*%Z%*%Pttm1[[i]] #compute P_{t|t}
        Pttm1[[i+1]] <- Tmatrix%*%Ptt[[i]]%*%t(Tmatrix)+Q #compute P_{t+1|t}
        xttm1[,i+1] <- Tmatrix%*%xtt[,i] #compute x_{t+1|t}
        
        #The optimization criterion
        if (outofsample) {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatoutofsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        } else {
          if (i <= (nstates-13) ){
            logl <- logl - nrow(y)/2*log(2*pi)
          } else if (i > (nstates-13) ){ # diffuse log likelihood
            logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)%*%epshatinsample
            if ((NaN %in% logl)==T){
              logl<- -P10[1]
            }
          }
        }
      }
    }
    if (opti) {
      return(-logl)
    }
    else {
      return(list(logl=-logl, xtt=xtt,xttm1=xttm1,Pttm1=Pttm1,Ptt=Ptt))
    }
  } #y1 and y2 have smooth trend and y3 has local level
  
}


# function to aggregate differences of flow variables
flow_diff <- function(x,k,per){  #x should be vector
  s <- nrow(aggregate(x,by=ym(per), sum))
  a <- matrix(NA,nrow=s, ncol=k)
  for (i in k:1){
    a[,i] <- aggregate(x,by=ym(per), function(x){rev(x)[i]})$x
    a[is.na(a)] <- 0
  }
  agg <- a%*%c(1:k) + rbind(0,a[-nrow(a),])%*%c((k-1):0)
  return(agg)
}

# function to aggregate flow variables
flow <- function(x,k,per){  #x should be vector
  s <- nrow(aggregate(x,by=ym(per), sum))
  a <- matrix(NA,nrow=s, ncol=k)
  for (i in k:1){
    a[,i] <- aggregate(x,by=ym(per), function(x){rev(x)[i]})$x
    a[is.na(a)] <- 0
  }
  agg <- a%*%rep(1,k)
  return(agg)
}

# function to aggregate flow variables by taking the average
flow_avg <- function(x,k,per){  #x should be vector
  s <- nrow(aggregate(x,by=ym(per), sum))
  a <- matrix(NA,nrow=s, ncol=k)
  for (i in k:1){
    a[,i] <- aggregate(x,by=ym(per), function(x){rev(x)[i]})$x
    a[is.na(a)] <- 0
  }
  agg <- a%*%rep(1,k)/k
  return(agg)
}

ym <- function(x) {list(format(x,format="%m"),format(x,format="%Y"))}



########## Import Data ###############

# Labour force on a monthly basis
if (until2017==TRUE){
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:109),]
  waves <- LFS[,c(1,43,85,127,169)]
} else if (f2007u2017==T){
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:73),]
  waves <- LFS[,c(1,43,85,127,169)]
} else if (until2012==TRUE){
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:49,145:204),]
  waves <- LFS[,c(1,43,85,127,169)]
} else if (until2013==TRUE){
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:61,157:204),]
  waves <- LFS[,c(1,43,85,127,169)]
} else if (until2014==TRUE){
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:73,169:204),]
  waves <- LFS[,c(1,43,85,127,169)]
} else if (until2015==TRUE){
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:85,181:204),]
  waves <- LFS[,c(1,43,85,127,169)]
} else if (until2016==TRUE){
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:97,193:204),]
  waves <- LFS[,c(1,43,85,127,169)]
} else if (monthlyGT == T){
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:37),]
  waves <- LFS[,c(1,43,85,127,169)]
} else if (f2004u2017 == T){
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:37),]
  waves <- LFS[,c(1,43,85,127,169)]
} else {
  LFS <- read_excel("DATA_LFS/lfs2017.xlsx", col_names = FALSE)
  LFS <- LFS[-c(1:37,133:204),]
  waves <- LFS[,c(1,43,85,127,169)]
}
waves <- as.data.frame(waves)
len.m <- nrow(waves)

# Claimant Counts on a monthly basis
if (until2017==TRUE){
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:109),3]*1000)
  CC <- as.data.frame(CC)
} else if (f2007u2017==T){
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:73),3]*1000)
  CC <- as.data.frame(CC)
} else if (until2012==TRUE){
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:49,145:204),3]*1000)
  CC <- as.data.frame(CC)
} else if (until2013==TRUE){
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:61,157:204),3]*1000)
  CC <- as.data.frame(CC)
} else if (until2014==TRUE){
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:73,169:204),3]*1000)
  CC <- as.data.frame(CC)
} else if (until2015==TRUE){
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:85,181:204),3]*1000)
  CC <- as.data.frame(CC)
} else if (until2016==TRUE){
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:97,193:204),3]*1000)
  CC <- as.data.frame(CC)
} else if (monthlyGT == T){
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:37),3]*1000)
  CC <- as.data.frame(CC)
} else if (f2004u2017 == T){
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:37),3]*1000)
  CC <- as.data.frame(CC)
} else {
  CC <- read_excel("DATA_LFS/cc2017.xlsx", col_names = FALSE)
  CC <- rbind(CC[-c(1:37,133:204),3]*1000)
  CC <- as.data.frame(CC)
}

# Google trends both on a weekly and on a monthly basis
if (until2017==TRUE){
  dataset <- read_excel("Data_Google_trends/2010_2017/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2010_2017/dataset_month.xlsx", col_names = T) # monthly
} else if (f2007u2017==T){
  dataset <- read_excel("Data_Google_trends/2007_2017/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2007_2017/dataset_month.xlsx", col_names = T) # monthly
} else if (until2012==TRUE){
  dataset <- read_excel("Data_Google_trends/2005_2012/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2005_2012/dataset_month.xlsx", col_names = T) # monthly
} else if (until2013==TRUE){
  dataset <- read_excel("Data_Google_trends/2006_2013/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2006_2013/dataset_month.xlsx", col_names = T) # monthly
} else if (until2014==TRUE){
  dataset <- read_excel("Data_Google_trends/2007_2014/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2007_2014/dataset_month.xlsx", col_names = T) # monthly
} else if (until2015==TRUE){
  dataset <- read_excel("Data_Google_trends/2008_2015/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2008_2015/dataset_month.xlsx", col_names = T) # monthly
} else if (until2016==TRUE){
  dataset <- read_excel("Data_Google_trends/2009_2016/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2009_2016/dataset_month.xlsx", col_names = T) # monthly
} else if (monthlyGT == T){
  dataset <- read_excel("Data_Google_trends/2004_2017/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2004_2017/dataset_month.xlsx", col_names = T) # monthly
} else if (f2004u2017 == T){
  dataset <- read_excel("Data_Google_trends/2004_2017/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2004_2017/dataset_month.xlsx", col_names = T) # monthly
} else {
  dataset <- read_excel("Data_Google_trends/2004_2011/dataset.xlsx", col_names = T) # weekly
  dataset_month <- read_excel("Data_Google_trends/2004_2011/dataset_month.xlsx", col_names = T) # monthly
}
x.w <- dataset[,c(2:ncol(dataset))]
x.w <- as.matrix(x.w)
x <- dataset_month[,c(2:ncol(dataset_month))]
x <- as.matrix(x)
# for (i in 1:ncol(x)){
#   plot(x.w[,i],type="l")
# }

# Aggregate weekly Google trends
dataset <- as.data.frame(dataset)
week <- as.factor(dataset[,1])
week <- as.Date(week,"%Y-%m-%d")
agg.gt <- matrix(NA,nrow(x),ncol(x))
agg.gt5 <- matrix(NA,nrow(x),ncol(x))
agg.gt.avg <- matrix(NA,nrow(x),ncol(x))
agg.gt5.avg <- matrix(NA,nrow(x),ncol(x))

# Aggregate as flow variables
for (i in 1:ncol(x.w)){
  agg.gt[,i] <- flow(x=x.w[,i], k=4, per=week) # flow variables
  agg.gt5 <- flow(x=x.w[,i], k=5, per=week)
  five.obs <- which(!is.na(aggregate(x.w[,i],by=ym(week),function(x){rev(x)[5]})$x))
  agg.gt[five.obs,i] <- agg.gt5[five.obs]
}
# Rescale according to largest value on the entire dataset
for (i in 1:ncol(agg.gt)){
  if (max(agg.gt[,i]) != 0){
    agg.gt[,i] <- agg.gt[,i]*100/max(agg.gt[,i])
  } else {
    agg.gt[,i] <- agg.gt[,i]
  }
}
for (i in 1:ncol(agg.gt)){
  print(max(agg.gt[,i]))
}
# Aggregate as flow variables by taking the average
for (i in 1:ncol(x.w)){
  agg.gt.avg[,i] <- flow_avg(x=x.w[,i], k=4, per=week) # flow variables
  agg.gt5.avg <- flow_avg(x=x.w[,i], k=5, per=week)
  five.obs.avg <- which(!is.na(aggregate(x.w[,i],by=ym(week),function(x){rev(x)[5]})$x))
  agg.gt.avg[five.obs.avg,i] <- agg.gt5.avg[five.obs.avg]
}
# Compare to monthly Google trends
for (i in ncol(x):1){
  ts.plot(cbind(x[,i],agg.gt[,i],agg.gt.avg[,i]), col=c("black","red","green"))
}

# Remove constant (are almost always constant) monthly and weekly Google trends
if (until2017==TRUE){
  x <- x[,-c(7,9,43)]
  x.w <- x.w[,-c(7,9,43)]
  agg.gt <- agg.gt[,-c(7,9,43)]
  agg.gt.avg <- agg.gt.avg[,-c(7,9,43)] 
} else if (f2007u2017==TRUE){
  x <- x[,-c(1,6,7,8,9,12,20,24,34,40,43,55,57,58,100)]
  x.w <- x.w[,-c(1,6,7,8,9,12,20,24,34,40,43,55,57,58,100)]
  agg.gt <- agg.gt[,-c(1,6,7,8,9,12,20,24,34,40,43,55,57,58,100)]
  agg.gt.avg <- agg.gt.avg[,-c(1,6,7,8,9,12,20,24,34,40,43,55,57,58,100)] 
} else if (until2013==TRUE){
  x <- x[,-c(1,2,6,7,8,9,12,20,34,40,57,58)]
  x.w <- x.w[,-c(1,2,6,7,8,9,12,20,34,40,57,58)]
  agg.gt <- agg.gt[,-c(1,2,6,7,8,9,12,20,34,40,57,58)]
  agg.gt.avg <- agg.gt.avg[,-c(1,2,6,7,8,9,12,20,34,40,57,58)]
} else if (until2014==TRUE){
  x <- x[,-c(1,6,7,8,9,12,20,24,34,40,57,58,96)]
  x.w <- x.w[,-c(1,6,7,8,9,12,20,24,34,40,57,58,96)]
  agg.gt <- agg.gt[,-c(1,6,7,8,9,12,20,24,34,40,57,58,96)]
  agg.gt.avg <- agg.gt.avg[,-c(1,6,7,8,9,12,20,24,34,40,57,58,96)] 
} else if (until2015==TRUE){
  x <- x[,-c(6,7,8,9,12,40,57)]
  x.w <- x.w[,-c(6,7,8,9,12,40,57)]
  agg.gt <- agg.gt[,-c(6,7,8,9,12,40,57)]
  agg.gt.avg <- agg.gt.avg[,-c(6,7,8,9,12,40,57)]
} else if (until2016==TRUE){
  x <- x[,-c(6,7,8,9,12,49,57)]
  x.w <- x.w[,-c(6,7,8,9,12,49,57)]
  agg.gt <- agg.gt[,-c(6,7,8,9,12,49,57)]
  agg.gt.avg <- agg.gt.avg[,-c(6,7,8,9,12,49,57)]
} else if (monthlyGT == T){
  x <- x[,-c(6,7,8,9,12,13,20,34,40,57)]
  x.w <- x.w[,-c(6,7,8,9,12,13,20,34,40,57)]
} else if (f2004u2017 == T){
  x <- x[,-c(6,7,8,9,12,13,20,23,24,20,30,34,40,41,57,58)]
  x.w <- x.w[,-c(6,7,8,9,12,13,20,23,24,20,30,34,40,41,57,58)]
  agg.gt <- agg.gt[,-c(6,7,8,9,12,13,20,23,24,20,30,34,40,41,57,58)]
  agg.gt.avg <- agg.gt.avg[,-c(6,7,8,9,12,13,20,23,24,20,30,34,40,41,57,58)]
} else {
  x <- x[,-c(1,6,7,8,9,12,20,34,40,57,96)]
  x.w <- x.w[,-c(1,6,7,8,9,12,20,34,40,57,96)]
  agg.gt <- agg.gt[,-c(1,6,7,8,9,12,20,34,40,57,96)]
  agg.gt.avg <- agg.gt.avg[,-c(1,6,7,8,9,12,20,34,40,57,96)]  
}

if (avg==T){
  agg.gt <- agg.gt.avg
}
if (monthlyGT==T){
  agg.gt <- x
}



#### Estimation of the LFS model, without auxiliary series ####
len.m <- nrow(waves)
y <- matrix(0,5,len.m)
for (j in 1:5){
  y[j,] <- waves[1:len.m,j]
}
#ts.plot(cbind(y[1,],y[2,],y[3,],y[4,],y[5,]),col=c("black","black","black","black","black"))

se <- as.data.frame(LFS[1:len.m,c(2,44,86,128,170)])
k <- se

objopt <-  optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                       log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2))),
                 KF_slopes_univ,y=y[1:5,],k=k,delta=0.21,opti=T,outofsample=T,parP10=1000000000000,nstates=30,  hessian=F, method="L-BFGS-B" )
par <- objopt$par
#exp(par)/(1-0.21^2) #to get the estimates for the standard deviation of the sampling errors
obj <- KF_slopes_univ(par=objopt$par,y=y[1:5,],k=k,delta=0.21,opti=F,outofsample=T,parP10=1000000000000,nstates=30)

#ts.plot(unlist(lapply(obj$Ptt, function(x) x[1,1]))[31:ncol(y)])
#ts.plot(cbind(y[1,],y[2,],y[3,],y[4,],y[5,],obj$xtt[1,]),col=c("black","black","black","black","black","red"))
#ts.plot(cbind(obj$xttm1[2,]),col=c("green"))

P.L <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]))[(31-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))])
P.R <- mean(unlist(lapply(obj$Ptt, function(x) x[2,2]))[(31-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
P.theta <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(31-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])

par.univ <- par

slope <- obj$xtt[2,(31-13):length(obj$xtt[2,])]




#### Estimation of the LFS model, with auxiliary series of CC ####
len.m <- nrow(waves)
y <- matrix(0,6,len.m)
for (j in 1:5){
  y[j,] <- waves[1:len.m,j]
}
y[6,] <- CC[1:len.m,1]
#ts.plot(cbind(y[1,],y[2,],y[3,],y[4,],y[5,],y[6,]),col=c("black","black","black","black","black","red"))

# Order of integration
# adfTest(diff(y[j,],2), lag=trunc((len.m-1)^(1/3)), type = "nc")

se <- as.data.frame(LFS[1:len.m,c(2,44,86,128,170)])
k <- se

objopt <-  optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                       log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),
                       log(3000),log(0.02),0,log(1000)),
                 KF_slopes,y=y,k=k,delta=0.21,opti=T,outofsample=T,parP10=1000000000000,nstates=43,  hessian=F, method="L-BFGS-B")
par <- objopt$par
#exp(par)/(1-0.21^2) #to get the estimates for the standard deviation of the sampling errors
obj <- KF_slopes(par=objopt$par,y=y,k=k,delta=0.21,opti=F,outofsample=T,parP10=1000000000000,nstates=43)

#ts.plot(unlist(lapply(obj$Ptt, function(x) x[31,31]))[44:ncol(y)])
#ts.plot(cbind(y[1,],y[2,],y[3,],y[4,],y[5,],obj$xtt[1,]),col=c("black","black","black","black","black","red"))
#ts.plot(cbind(obj$xttm1[2,]),col=c("green"))
#ts.plot(cbind(y[6,],obj$xtt[31,]),col=c("black","red"))
#ts.plot(cbind(obj$xttm1[32,]),col=c("green"))

#estimated seasonal components
season_y_est_tt <- obj$xtt[3,]+obj$xtt[5,]+obj$xtt[7,]+obj$xtt[9,]+obj$xtt[11,]+obj$xtt[13,]
season_y_est_ttm1 <- obj$xttm1[3,]+obj$xttm1[5,]+obj$xttm1[7,]+obj$xttm1[9,]+obj$xttm1[11,]+obj$xttm1[13,]
season_x_est_tt <- obj$xtt[33,]+obj$xtt[35,]+obj$xtt[37,]+obj$xtt[39,]+obj$xtt[41,]+obj$xtt[43,]
season_x_est_ttm1 <- obj$xttm1[33,]+obj$xttm1[35,]+obj$xttm1[37,]+obj$xttm1[39,]+obj$xttm1[41,]+obj$xttm1[43,]
#ts.plot(cbind(season_y_est_tt,season_y_est_ttm1[1:len.m]),lty=c(1,2), col=c("red","red"))
#ts.plot(cbind(season_x_est_tt,season_x_est_ttm1[1:len.m]),lty=c(1,2), col=c("red","red"))

#estimated theta
theta_y_est_tt <- obj$xtt[1,] + season_y_est_tt 
theta_y_est_ttm1 <- obj$xttm1[1,] + season_y_est_ttm1 
theta_x_est_tt <- obj$xtt[31,] + season_y_est_tt 
theta_x_est_ttm1 <- obj$xttm1[31,] + season_y_est_ttm1 
#ts.plot(cbind(y[1,],y[2,],y[3,],y[4,],y[5,],theta_y_est_tt,theta_y_est_ttm1[1:len.m]),lty=c(1,1,1,1,1,1,2), col=c("black","black","black","black","black","red","red"))
#ts.plot(cbind(y[6,],theta_x_est_tt,theta_x_est_ttm1[1:len.m]),lty=c(1,1,2), col=c("black","red","red"))

P.L.biv <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]))[(44-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))])
P.R.biv <- mean(unlist(lapply(obj$Ptt, function(x) x[2,2]))[(44-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
P.theta.biv <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                    2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(44-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])

par.biv <- par


# Diagnostic tests for the forecast errors

eta <- obj$st.for[,(44-13):ncol(obj$st.for)]

# serial correlation
p.serial.Ljung <- rep(NA,6)
for (i in 1:6){
  if (i==1) {
    p.serial.Ljung[i] <- Box.test(eta[i,], lag =4, type = "Ljung-Box", fitdf = 2)$p.value
  } else {
    p.serial.Ljung[i] <- Box.test(eta[i,], lag =4, type = "Ljung-Box", fitdf = 3)$p.value
    if (i==6){
      p.serial.Ljung[i] <- Box.test(eta[i,], lag =4, type = "Ljung-Box", fitdf = 1)$p.value
    }
  }
}


# univariate normality
p.normal.SW <- rep(NA,6)
p.normal.CVM <- rep(NA,6)
p.normal.Lillie <- rep(NA,6)
p.normal.SF <- rep(NA,6)
p.normal.AD <- rep(NA,6)
p.normal.JB <- rep(NA,6)
for (i in 1:6){
  p.normal.SW[i] <- mvn(data = eta[i,], univariateTest = "SW")$univariateNormality$`p value`
  p.normal.CVM[i] <- mvn(data = eta[i,], univariateTest = "CVM")$univariateNormality$`p value`
  p.normal.Lillie[i] <- mvn(data = eta[i,], univariateTest = "Lillie")$univariateNormality$`p value`
  p.normal.SF[i] <- mvn(data = eta[i,], univariateTest = "SF")$univariateNormality$`p value`
  p.normal.AD[i] <- mvn(data = eta[i,], univariateTest = "AD")$univariateNormality$`p value`
  p.normal.JB[i] <- jarque.bera.test(eta[i,])$p.value
}

p.normal <- cbind(p.normal.SW, p.normal.CVM, p.normal.Lillie, p.normal.SF, p.normal.AD,p.normal.JB)

#heteroscedasticity
half <- trunc(ncol(eta)/2)
p.heter <- rep(NA,6)
F.cv <- rep(NA,6)
for (i in 1:6){
  F.cv[i] <- sum((eta[i,(ncol(eta)-half+1):ncol(eta)])^2)/sum((eta[i,1:half])^2)
  p.heter[i] <- pf(F.cv[i], half, half, lower.tail = F, log.p = FALSE)
}



#### Test for non-stationarity in Google trends ####

# Differenciate data
diff.x <- diff(agg.gt,1) 
diff.x.w <- diff(x.w,1)

# Control for multiple hypotheses testing
n <- nrow(diff.x)
s <- ncol(diff.x) #number of t statistics

p.adf <- rep(NA,s)
test.adf <- rep(NA,s)
for (i in 1:s) { 
  test <- adfTest(agg.gt[2:nrow(agg.gt),i], lag=trunc((n-1)^(1/3)), type = "ct")
  p.adf[i] <- test@test$p.value
  test.adf[i] <- test@test$statistic
}

test.adf <- -test.adf
tis <- sort(test.adf, decreasing=F) # sort the statistics from lowest significant to highest significant, and append a vector that saves their original index
tis <- cbind(tis,sort(test.adf, decreasing=F, index.return=T)$ix)

block <- rep(NA, s)
nr.bl <- rep(NA, s)
for (i in 1:s){
  block[i] <- trunc(getBandwidthAnd(diff.x[,i], kernel = "ba", check = TRUE)) # optimally obtain block length
  if (block[i] == 0){
    block[i] <- 1
  }
  nr.bl[i] <- ceiling(n/block[i])
}
B <- 999
alpha <- 0.05
A <- diff.x

store_t <- matrix(0,nrow=B,ncol=s)
store_p <- matrix(0,nrow=B,ncol=s)

# Bootstrap loop
set.seed(43)
for (b in 1:B) {
  A.star <- matrix(0,nrow=n,ncol=s)
  for (i in 1:s){
    v.star <- rep(0,n)
    start.b <- sample.int(n-block[i]+1, size=nr.bl[i], replace=T) # draw the starting points of the blocks
    for (j in 1:nr.bl[i]) {
      if (j < nr.bl[i] || block[i]==1 || n %% block[i]==0) {
        index.j <- 1:block[i] + start.b[j] - 1
        v.star <- c(v.star,A[index.j,i]) # draw the bootstrap sample
      } else {
        index.j <- 1:(n - trunc(n/block[i])*block[i]) + start.b[j] - 1
        v.star <- c(v.star,A[index.j,i]) # draw the bootstrap sample
      }
    }
    v.star <- v.star[-1:-n] # eliminate elements which are 0
    A.star[,i] <- v.star
  }
  for (j in 1:ncol(A.star)){ # bootstrap sample of level variables
    A.star[,j] <- cumsum(A.star[,j]) 
  }
  for (j in 1:s) {
    test <- adfTest(A.star[,j], lag=trunc((n-1)^(1/3)), type = "ct")
    store_t[b,j] <- test@test$statistic # store bootstrap t statistiscs
    store_p[b,j] <- test@test$p.value # store bootstrap p-values
  }
}

# FDR loop
c <- rep(0,s) # initialize vector of critical values c_1, ..., c_s
for (j in 1:s){
  rej_j <- matrix(NA,B,(B+1))
  rej_j_only <- matrix(NA,B,(B+1))
  fdr <- rep(0,(B+1)) # initialize vector to store FDR for all possible (B+1) critical values
  Tb_j <- matrix(NA,nrow=B,ncol=j)
  Tb_j <- -store_t[,tis[1:j,2]] # select only the bootstrap test statistics for those units corresponding to the set of j least significant tests
  if (j==1){ # Romano et al. (2008) equation (8)
    Tb_sort_j <- Tb_j  
    cv <- c(sort(Tb_sort_j, decreasing=T),-99999) # take the most significant statistics across all bootstrap replications, and sort these from high to low - these are all possible critical values sorted in descending order; -99999 represents minus infinity
    for (i in 1:(B+1)){
      for (l in 1:B){
        rej_j[l,i] <- Tb_sort_j[l] >= cv[i]
      }
    }
    fdr <- apply(rej_j,2,mean)/s # calculate the bootstrap FDR for each potential critical value as the probability of rejection of T^*_j:j (average over all replications) divided by the number of rejections
    c[j] <- cv[sum(fdr <= alpha)] # take c_1 as the minimal cv that assures FDR <= alpha
  } else { # Romano et al. (2008) equations (9) and (10)
    Tb_sort_j <- t(apply(Tb_j, 1, sort, decreasing=F)) # for each bootstrap replication, sort the statistics from low to high
    cv <- c(sort(Tb_sort_j[,j], decreasing=T),-99999) # take the most significant statistics across all bootstrap replications, and sort these from high to low - these are all possible critical values sorted in descending order; -99999 represents minus infinity
    for (i in 1:(B+1)){
      for (l in 1:B){
        rej_j[l,i] <- Tb_sort_j[l,j] >= cv[i]
      }
    }
    for (r in (s-j+1):s){ # sum over r, equation (10)
      if (r==(s-j+1)){
        rej_j_only <- (Tb_sort_j[,(j-1)] < c[j-1])*rej_j
        fdr <- (r-s+j)/r*apply(rej_j_only,2,mean) # calculate the bootstrap FDR for each potential critical value as the probability of rejection of T^*_j:j while not rejecting T^*_j-1:j (average over all replications) divided by the number of rejections
      } else if (r > (s-j+1) && r < s) {
        selb <- as.matrix(Tb_sort_j[,(s-r+1):(j-1)])
        selcv <- c[(s-r+1):(j-1)]
        rej_j_to_r <- rep(1,(B))
        for (i in 1:ncol(selb)){
          selb[,i] <- selb[,i] >= selcv[i]
          rej_j_to_r <- rej_j_to_r*selb[,i]
        }
        rej_j_to_r <- rej_j_to_r*(Tb_sort_j[,(s-r)] < c[s-r])*rej_j
        fdr <- fdr + (r-s+j)/r*apply(rej_j_to_r,2,mean) 
      } else {
        selb <- as.matrix(Tb_sort_j[,(s-r+1):(j-1)])
        selcv <- c[(s-r+1):(j-1)]
        rej_j_to_r <- rep(1,(B))
        for (i in 1:ncol(selb)){
          selb[,i] <- selb[,i] >= selcv[i]
          rej_j_to_r <- rej_j_to_r*selb[,i]
        }
        rej_j_to_r <- rej_j_to_r*rej_j
        fdr <- fdr + (r-s+j)/r*apply(rej_j_to_r,2,mean)
      }
    }
    c[j] <- cv[sum(fdr <= alpha)]
  }
}

-rev(tis[,1]) < -rev(c) #test sequentially from most significant to lowest (stepdown procedure)
for (i in 1:s){
  if(-rev(tis[,1])[i] >= -rev(c)[i]) break
  j.star <- i+1
}
sort(-test.adf) < -2.89 # compare with ADF test without controlling for MHT (5% cl for sample size of 100)
nonstat <- x[,rev(tis[,2])[j.star:s]]
stat <- x[,-rev(tis[,2])[j.star:s]]
# plot non-stationary series
for (i in ncol(nonstat):1){
  ts.plot(cbind(nonstat[,i]))
  #acf(nonstat[,i])
}
# plot stationary series
for (i in ncol(stat):1){
  ts.plot(cbind(stat[,i]))
  #acf(stat[,i])
}

# Keep non-stationary Google trends:
agg.gt <- agg.gt[,rev(tis[,2])[j.star:s]] 
x.w <- x.w[,rev(tis[,2])[j.star:s]]
x <- x[,rev(tis[,2])[j.star:s]]



#### Target the predictors ####

# Keep only GTs highly correlated with \hat{R}^y_t
keep <- rep(NA, ncol(agg.gt))
for (i in 1:ncol(agg.gt)){
  if (abs(cor(agg.gt[(32-13):nrow(agg.gt),i],slope)) >= 0.4){
    keep[i] <- i
    ts.plot(cbind(agg.gt[(32-13):nrow(agg.gt),i]),col="red")
    par(new=len.m)
    ts.plot(cbind(slope),col="black")
    axis(side=4)
  }
}
length(keep[which(!is.na(keep))])
# agg.gt <- agg.gt[,keep[which(!is.na(keep))]]
# x <- x[,keep[which(!is.na(keep))]]
# x.w <- x.w[,keep[which(!is.na(keep))]]

# Deseasonalize Google trends
if (desGT==T){
gt.des <- agg.gt
  for (i in 1:ncol(agg.gt)){
    # gt.des[,i] <- ds(agg.gt[,i], type = "monthly", lag.max=12, ic="AIC", standardizeQ=F)$z
    ts_beer = ts(agg.gt[,i], frequency = 12, start = 2004)
    decompose_beer = decompose(ts_beer, "additive")
    gt.des[,i] = ts_beer - decompose_beer$seasonal
    ts.plot(cbind(agg.gt[,i],gt.des[,i]),col=c("black","red"))
  }
  agg.gt <- gt.des
}

if (EN ==T){ 
  # Elastic Net
  diff.x <- diff(agg.gt[(32-13):nrow(agg.gt),],1)
  diff.y <- diff(slope,1)
  for (i in 1:ncol(diff.x)){
    diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
  }
  diff.y <- (diff.y-mean(diff.y))/sd(diff.y)
  alpha_grid <- seq(0.05, 0.95, 0.01) # grid values for alpha
  AIC <- rep(NA,length(alpha_grid))
  tune_lambda <- rep(NA,length(alpha_grid))
  nselect_var <- rep(NA,length(alpha_grid))
  for (i in 1:length(alpha_grid)) {
    en <- ic.glmnet(diff.x, diff.y, crit = "aic", alpha=alpha_grid[i], lambda=seq(0.1, 0, -0.001)) #alpha=1 lasso, alpha=0 ridge
    # plot(en$glmnet,"lambda")
    # plot(en)
    tune_lambda[i] <- en$lambda
    AIC[i] <- en$ic[2]
    nselect_var[i] <-en$nvar-1
    # coef(en)
    # en$glmnet$lambda
  }
  nselect_var[which.min(AIC)]
  fit <- glmnet(diff.x, diff.y, lambda = tune_lambda[which.min(AIC)], alpha = alpha_grid[which.min(AIC)])
  coef(fit)
  sel <- fit$beta@i[-1]+1
  for (i in sel){
    ts.plot(cbind(agg.gt[(32-13):nrow(agg.gt),i]),col="red")
    par(new=len.m)
    ts.plot(cbind(slope),col="black")
    axis(side=4)
  }
  hhh <- agg.gt
  agg.gt <- agg.gt[,sel]
  #x <- x[,sel]
  #x.w <- x.w[,sel]
}


if (PLS==T){
  # Partial Least Squares
  diff.x <- diff(agg.gt[32:nrow(agg.gt),],1)
  diff.y <- diff(slope,1)
  for (i in 1:ncol(diff.x)){
    diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
  }
  diff.y <- (diff.y-mean(diff.y))/sd(diff.y)
  
  # plsdepot package
  pls1 <- plsreg1(diff.x, diff.y, comps = 10, crosval = F)
  pls1$x.scores%*%t(pls1$x.loads)
  
  for (i in 1:ncol(pls1$x.scores)){
    ts.plot(cbind(cumsum(pls1$x.scores[,i])),col="red")
    par(new=len.m)
    ts.plot(cbind(slope[-1]),col="black")
    axis(side=4)
  }
  
  # pls package
  pls2 <- plsr(diff.y ~ diff.x, ncomp = 8, validation = "none")
  pls2$scores%*%t(pls2$loadings)
  plot(RMSEP(pls2))
  
  for (i in 1:ncol(pls2$scores)){
    ts.plot(cbind(cumsum(pls2$scores[,i])),col="red")
    par(new=len.m)
    ts.plot(cbind(slope),col="black")
    axis(side=4)
  }
  
  for (i in 1:ncol(pls2$scores)){
    print(cor(cumsum(pls2$scores[,i]),slope[-1]))
  }
}


#### Get the factor for the monthly selected predictors with PCA ####

if (lagGT==T){
  lag.agg.gt <- agg.gt[13:nrow(agg.gt),]
  for (i in 1:12){
    lag.agg.gt <- cbind(lag.agg.gt,agg.gt[(13-i):(nrow(agg.gt)-i),])
  }
  agg.gt <- lag.agg.gt
}

# Differenciate data
diff.x <- diff(agg.gt,1) 
diff.x.w <- diff(x.w,1)

if (EN==T){
  diff.x.w <- diff(x.w[,sel],1)
}

if (load.w==T){
  diff.x <- diff.x.w
}

# Standardize the monthly Google trends
for (i in 1:ncol(diff.x)){
  diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
}
diff.x <- as.matrix(diff.x)
kmax <- 10 # max number of factors
cov <- diff.x%*%t(diff.x)
len.m <- nrow(diff.x)
nvar <- ncol(diff.x)
kmax.eigenvalues <- eigen(cov)$values[1:kmax]
kmax.eigenvectors <- eigen(cov)$vectors[,1:kmax] 
F.hat.max <- kmax.eigenvectors*sqrt(len.m-1) # facotrs
t(F.hat.max)%*%F.hat.max/(len.m-1)
Lambda.hat.max <- t(diff.x)%*%F.hat.max/(len.m-1) # loadings
t(Lambda.hat.max)%*%Lambda.hat.max
Comps.hat.max <- F.hat.max%*%t(Lambda.hat.max) # common components
e.max <- diff.x-Comps.hat.max # estimated hydiosincratic components (hc)

# We need the kmax variance of the hc for the information criteria (IC)
sigma.hat.maxj <- c(1:nvar)
Vj <- sigma.hat.maxj
for (j in 1:nvar){
  sigma.hat.maxj[j] <- 1/len.m*t(e.max[,j])%*%(e.max[,j])
}
sigma.hat.max <- 1/nvar*sum(sigma.hat.maxj)

# Define vectors to store the values for the IC (Bai and Ng)
PC1 <- c(1:kmax)
PC2 <- PC1
PC3 <- PC1
IC1 <- PC1
IC2 <- PC1
IC3 <- PC1

# Loop to choose the number of factors based on IC
for (i in 1:kmax){
  k.eigenvalues <- eigen(cov)$values[1:i]
  k.eigenvectors <- eigen(cov)$vectors[,1:i]
  F.hat <- k.eigenvectors*sqrt(len.m-1)
  Lambda.hat <- t(diff.x)%*%F.hat/(len.m-1)
  Comps.hat <- F.hat%*%t(Lambda.hat)
  e <- diff.x-Comps.hat
  # Try with differenced data
  # e <- diff(e)
  for (j in 1:nvar){
    Vj[j] <- 1/len.m*t(e[,j])%*%(e[,j])
  }
  V <- 1/nvar*sum(Vj)
  # Information Criteria for stationary factors
  PC1[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
  PC2[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
  PC3[i] <- V + i*sigma.hat.max*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
  IC1[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
  IC2[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
  IC3[i] <- log(V) + i*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
}

# Information Criteria
k.PC1 <- which.min(PC1)
k.PC2 <- which.min(PC2)
k.PC3 <- which.min(PC3)
k.IC1 <- which.min(IC1)
k.IC2 <- which.min(IC2)
k.IC3 <- which.min(IC3)

# Double-check with a built-in function
a <- prcomp(diff.x, rank.=kmax)
screeplot(a)
plot(sqrt(kmax.eigenvalues))

#Non-stationary standardized Google trends
stand.x <- diff.x
for (i in 1:ncol(diff.x)){
  stand.x[,i] <- cumsum(diff.x[,i])
}
stand.x <- as.matrix(stand.x)

# Non-stationary monthly factors
F.ns.hat <- F.hat.max
for (i in 1:kmax){
  F.ns.hat[,i] <- cumsum(F.hat.max[,i])
  plot(F.ns.hat[,i], type="l")
}
# covariance matrix of the idiosyncratic components:
idio <- stand.x - F.ns.hat[,1]%*%t(Lambda.hat.max[,1])
H <- 1/(len.m-1)*t(idio)%*%idio

ts.plot(cbind(F.ns.hat[,1]), col=c("red"))
par(new=len.m)
ts.plot(cbind(y[1,1:len.m],y[2,1:len.m],y[3,1:len.m],y[4,1:len.m],y[5,1:len.m]),col=c("black","black","black","black","black"))
axis(side = 4)

ts.plot(cbind(F.ns.hat[,1]), col=c("red"))
par(new=len.m)
ts.plot(cbind(y[1,(nrow(waves)-nrow(stand.x)+1):nrow(waves)]),col=c("black"))
axis(side = 4)

ts.plot(cbind(F.ns.hat[2:len.m,1]), col=c("red"))
par(new=len.m)
ts.plot(cbind(diff(y[1,(nrow(waves)-nrow(stand.x)+1):nrow(waves)],1)),col=c("black"))
axis(side = 4)

for (i in 1:kmax){
  print(paste(cor(slope,F.ns.hat[31:length(F.ns.hat[,i]),i])))
  ts.plot(cbind(F.ns.hat[31:length(F.ns.hat[,i]),i]), col=c("red"))
  par(new=len.m)
  ts.plot(cbind(slope),col=c("black"))
  axis(side = 4)
}

if (load.w==T){
  diff.x <- diff(agg.gt,1)
  # Standardize the monthly Google trends
  for (i in 1:ncol(diff.x)){
    diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
  }
  diff.x <- as.matrix(diff.x)
  stand.x <- diff.x
  for (i in 1:ncol(diff.x)){
    stand.x[,i] <- cumsum(diff.x[,i])
  }
  stand.x <- as.matrix(stand.x)
}

if (PLS==T){
  Lambda.hat.max[,1] <- pls2$loadings[,1]
  idio <- stand.x[32:nrow(stand.x),] - pls2$scores[,1]%*%t(pls2$loadings[,1])
  H <- 1/(nrow(idio)-1)*t(idio)%*%idio
  diff.x <- diff(agg.gt[32:nrow(agg.gt),],1)
  for (i in 1:ncol(diff.x)){
    diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
  }
  stand.x <- diff.x
  for (i in 1:ncol(diff.x)){
    stand.x[,i] <- cumsum(diff.x[,i])
  }
  stand.x <- as.matrix(stand.x)
}



#### Estimation of the LFS model, with auxiliary series of Google trends ####

if (CCGT == F) {
  # Weekly Google trends
  # Temporally disaggregate the waves
  # merged <- matrix(rep(NA, length(F.hat.max[,1])*5), ncol=5, byrow=TRUE)
  # merged <- cbind(F.hat.max[,1], merged)
  # merged <- as.data.frame(merged)
  # merged <- cbind(merged, week)
  # 
  # for (i in 1:length(week)){
  #   for (j in 1:length(month)){
  #     if (merged[i,7]==month[j]){
  #       merged[i,2] <- waves[j,1]
  #       merged[i,3] <- waves[j,2]
  #       merged[i,4] <- waves[j,3]
  #       merged[i,5] <- waves[j,4]
  #       merged[i,6] <- waves[j,5]
  #     } else {
  #       merged[i,2] <- merged[i,2]
  #       merged[i,3] <- merged[i,3]
  #       merged[i,4] <- merged[i,4]
  #       merged[i,5] <- merged[i,5]
  #       merged[i,6] <- merged[i,6]
  #     }
  #   }
  # }
  # for (i in 1:(ncol(merged)-1)){
  #   merged[,i] <- as.numeric(na.locf(merged[,i], fromLast = TRUE))
  # }
  # y <- t(cbind(merged[,2:6],sw))
  
  # Monthly google trends
  len.m <- nrow(waves)
  merged <- matrix(0,len.m,5)
  for (j in 1:5){
    merged[,j] <- waves[1:len.m,j]
  }
  y <- t(cbind(merged,rbind(matrix(NA,(nrow(merged)-nrow(stand.x)),ncol(stand.x)),stand.x)))
  
  # Kalman filter estimation
  se <- as.data.frame(LFS[,c(2,44,86,128,170)])
  k <- se
  if (seas_fact==F){
    objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                          log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),0), 
                    KF_slopes_mixed_factor,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                    parP10=1000000000000,nstates=31,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
    par <- objopt$par
    #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
    obj <- KF_slopes_mixed_factor(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                  nstates=31,lambda=Lambda.hat.max[,1],H=H)
    
    #ts.plot(unlist(lapply(obj$Ptt, function(x) x[1,1]))[(ncol(y)-nrow(stand.x)+31):ncol(y)])
    #ts.plot(cbind(y[1,],y[2,],y[3,],y[4,],y[5,],obj$xtt[1,]),col=c("black","black","black","black","black","red"))
    #ts.plot(cbind(obj$xttm1[2,-c(1:19)]),col=c("green"))
    if (PLS == T){
      ts.plot(cbind(cumsum(pls2$scores[,1]),obj$xtt[31,(32-13):(len.m)]),col=c("black","red"))
    } else {
      ts.plot(cbind(F.ns.hat[(31-13):len.m,1],obj$xtt[31,(ncol(y)-nrow(stand.x)+31-13):(len.m)]),col=c("black","red"))
    }
    
    ts.plot(cbind(obj$xtt[2,(ncol(y)-nrow(stand.x)+31-13):(len.m-1)]), col=c("green"))
    par(new=len.m)
    ts.plot(cbind(obj$xtt[31,(ncol(y)-nrow(stand.x)+31-13):(len.m-1)]),col=c("red"))
    axis(side = 4)
    
    P.L.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]))[(ncol(y)-nrow(stand.x)+31-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))])
    P.R.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[2,2]))[(ncol(y)-nrow(stand.x)+31-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    P.theta.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                       2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(ncol(y)-nrow(stand.x)+31-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    
    par.gt <- par
    
    # Test for normaility and independence in the residuals
    #Build T:
    Tymu <- matrix(c(1,1,0,1),2,2, byrow=T)
    C <- array(0,dim=c(2,2,5))
    for (l in 1:5){
      C[,,l] <- matrix(c(cos((pi*l)/6),  sin((pi*l)/6), -sin((pi*l)/6), cos((pi*l)/6)),2,2,byrow=TRUE)
    }
    Tyomega <- adiag(C[,,1],C[,,2],C[,,3],C[,,4],C[,,5],-1)
    ncol(Tyomega)
    nrow(Tyomega)
    Tylambda <- diag(4)
    TyE <- rbind(matrix(0,9,5), cbind(diag(4), c(0,0,0,0)))
    delta <- delta
    TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
    TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
    Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
    Tx <- as.matrix(1)
    Tmatrix <- adiag(Ty, Tx)
    
    eta <- obj$xtt[,(33-13):ncol(obj$xtt)] - Tmatrix%*%obj$xtt[,(32-13):(ncol(obj$xtt)-1)] # estimated residuals
    eta <- eta[-c(1,23:30),] # remove residuals that should be zero
    for (i in nrow(eta):1){
      if (Box.test(eta[i,], lag = 12, type = "Box-Pierce", fitdf = 0)$p.value < 0.05)
        print(i) # print if the null hypothesis of independence is rejected
    }
    for (i in nrow(eta):1){
      if (mvn(data = eta[i,], univariateTest = "Lillie")$univariateNormality$`p value` < 0.05)
        print(i) # print if the null hypothesis of univariate normality is rejected
    }
    mvn(data = t(eta[c(1,22),]), mvnTest = "hz")$multivariateNormality$`p value` < 0.05 # TRUE if the null hypothesis of multivariate normality is rejected
    mvn(data = t(eta[c(1,22),]), mvnTest = "royston")$multivariateNormality$`p value` < 0.05 # TRUE if the null hypothesis of multivariate normality is rejected
    mvn(data = t(eta[c(1,22),]), mvnTest = "dh")$multivariateNormality$`p value` < 0.05 # TRUE if the null hypothesis of multivariate normality is rejected
    mvn(data = t(eta[c(1,22),]), mvnTest = "energy")$multivariateNormality$`p value` < 0.05 # TRUE if the null hypothesis of multivariate normality is rejected
    
  } else if (seas_fact==T && fixed_sigma_u==T) {
    objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                          log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(0.02),0), 
                    KF_slopes_mixed_factor_seas,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                    parP10=1000000000000,nstates=42,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
    par.gt <- objopt$par
    #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
    obj <- KF_slopes_mixed_factor_seas(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                  nstates=42,lambda=Lambda.hat.max[,1],H=H)
  
    P.L.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]))[(ncol(y)-nrow(stand.x)+43-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))])
    P.R.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[2,2]))[(ncol(y)-nrow(stand.x)+43-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    P.theta.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                       2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(ncol(y)-nrow(stand.x)+43-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    
  } else if (seas_fact==T && fixed_sigma_u==F) {
    objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                          log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(1),log(1),0), 
                    KF_slopes_mixed_factor_seas,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                    parP10=1000000000000,nstates=42,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
    par.gt <- objopt$par
    #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
    obj <- KF_slopes_mixed_factor_seas(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                       nstates=42,lambda=Lambda.hat.max[,1],H=H)
    
    P.L.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]))[(ncol(y)-nrow(stand.x)+43-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))])
    P.R.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[2,2]))[(ncol(y)-nrow(stand.x)+43-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    P.theta.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                       2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(ncol(y)-nrow(stand.x)+43-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    
  }
}

#### Estimation of the LFS model, with auxiliary series of Claimant Counts and Google trends ####

if (CCGT == TRUE) {
  # Monthly google trends
  len.m <- nrow(waves)
  merged <- matrix(0,len.m,6)
  for (j in 1:5){
    merged[,j] <- waves[1:len.m,j]
  }
  merged[,6] <- CC[1:len.m,1]
  y <- t(cbind(merged,rbind(matrix(NA,(nrow(merged)-nrow(stand.x)),ncol(stand.x)),stand.x)))
  
  # Kalman filter estimation
  se <- as.data.frame(LFS[,c(2,44,86,128,170)])
  k <- se
  
  if (seas_fact==F){
    
    objopt <-  optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                           log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(3000),log(0.02),0,0,log(1000)),
                     KF_slopes_mixed_factor_CC,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                     parP10=1000000000000,nstates=44,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
    par <- objopt$par
    #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
    obj <- KF_slopes_mixed_factor_CC(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                     nstates=44,lambda=Lambda.hat.max[,1],H=H)
    
    #ts.plot(unlist(lapply(obj$Ptt, function(x) x[1,1]))[45:ncol(y)])
    #ts.plot(cbind(y[1,],y[2,],y[3,],y[4,],y[5,],obj$xtt[1,]),col=c("black","black","black","black","black","red"))
    #ts.plot(cbind(obj$xttm1[2,]),col=c("green"))
    ts.plot(cbind(F.ns.hat[(45-13):len.m,1],obj$xtt[44,(ncol(y)-nrow(stand.x)+45-13):(len.m)]),col=c("black","red"))
    
    ts.plot(cbind(obj$xtt[2,(45-13):(len.m-1)]), col=c("green"))
    par(new=len.m)
    ts.plot(cbind(obj$xtt[44,(45-13):(len.m-1)]),col=c("red"))
    axis(side = 4)
    
    P.L.cc.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]))[(45-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))])
    P.R.cc.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[2,2]))[(45-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    P.theta.cc.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                          2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(45-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    
    par.cc.gt <- par
  
    
    # Diagnostic tests for the forecast errors
    
    eta <- obj$st.for[,(45-13):ncol(obj$st.for)]
    
    # serial correlation
    p.serial.Ljung <- rep(NA,6)
    for (i in 1:6){
      if (i==1) {
        p.serial.Ljung[i] <- Box.test(eta[i,], lag =16, type = "Ljung-Box", fitdf = 2)$p.value
      } else {
        p.serial.Ljung[i] <- Box.test(eta[i,], lag =16, type = "Ljung-Box", fitdf = 3)$p.value
        if (i==6){
        p.serial.Ljung[i] <- Box.test(eta[i,], lag =16, type = "Ljung-Box", fitdf = 1)$p.value
        }
      }
    }
    
    
    # univariate normality
    p.normal.SW <- rep(NA,6)
    p.normal.CVM <- rep(NA,6)
    p.normal.Lillie <- rep(NA,6)
    p.normal.SF <- rep(NA,6)
    p.normal.AD <- rep(NA,6)
    p.normal.JB <- rep(NA,6)
    for (i in 1:6){
      p.normal.SW[i] <- mvn(data = eta[i,], univariateTest = "SW")$univariateNormality$`p value`
      p.normal.CVM[i] <- mvn(data = eta[i,], univariateTest = "CVM")$univariateNormality$`p value`
      p.normal.Lillie[i] <- mvn(data = eta[i,], univariateTest = "Lillie")$univariateNormality$`p value`
      p.normal.SF[i] <- mvn(data = eta[i,], univariateTest = "SF")$univariateNormality$`p value`
      p.normal.AD[i] <- mvn(data = eta[i,], univariateTest = "AD")$univariateNormality$`p value`
      p.normal.JB[i] <- jarque.bera.test(eta[i,])$p.value
    }
    
    p.normal <- cbind(p.normal.SW, p.normal.CVM, p.normal.Lillie, p.normal.SF, p.normal.AD,p.normal.JB)
    
    #heteroscedasticity
    half <- trunc(ncol(eta)/2)
    p.heter <- rep(NA,6)
    F.cv <- rep(NA,6)
    for (i in 1:6){
      F.cv[i] <- sum((eta[i,(ncol(eta)-half+1):ncol(eta)])^2)/sum((eta[i,1:half])^2)
      p.heter[i] <- pf(F.cv[i], half, half, lower.tail = F, log.p = FALSE)
    }
    
    #correct for multiple hypotheses testing
    p.serial <- c(p.serial.Box, p.serial.Ljung)
    p.correct.serial <- p.adjust(p.serial.Ljung, method = "BY", n = length(p.serial.Ljung))
    for (i in 1:length(p.correct.serial)){
      if(p.correct.serial[i] <= 0.05 )
        print(p.correct.serial[i])
    }
    which(p.correct.serial <= 0.05)
    
    p.correct.normal <- p.adjust(p.normal[,1], method = "BY", n = length(p.normal[,1]))
    for (i in 1:length(p.correct.normal)){
      if(p.correct.normal[i] <= 0.05 )
        print(p.correct.normal[i])
    }
    which(p.correct.normal <= 0.05)
    
    p.correct.heter <- p.adjust(p.heter, method = "BY", n = length(p.heter))
    for (i in 1:length(p.correct.heter)){
      if(p.correct.heter[i] <= 0.05 )
        print(p.correct.heter[i])
    }
    which(p.correct.heter <= 0.05)
    
    
    for (i in 1:6){
      plot(eta[i,],type="l")
    }
    for (i in 1:6){
      acf(eta[i,])
    }
  
    
  } else if (seas_fact==T && fixed_sigma_u==T) {
    objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                          log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(3000),log(0.02),log(0.02),0,0,log(1000)), 
                    KF_slopes_mixed_factor_seas_CC,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                    parP10=1000000000000,nstates=55,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
    par.cc.gt <- objopt$par
    #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
    obj <- KF_slopes_mixed_factor_seas_CC(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                       nstates=55,lambda=Lambda.hat.max[,1],H=H)
    
    P.L.cc.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]))[(ncol(y)-nrow(stand.x)+56-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))])
    P.R.cc.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[2,2]))[(ncol(y)-nrow(stand.x)+56-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    P.theta.cc.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                       2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(ncol(y)-nrow(stand.x)+56-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    
  } else if (seas_fact==T && fixed_sigma_u==F) {
    objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                          log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(3000),log(0.02),log(1),log(1),0,0,log(1000)), 
                    KF_slopes_mixed_factor_seas_CC,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                    parP10=1000000000000,nstates=55,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
    par.cc.gt <- objopt$par
    #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
    obj <- KF_slopes_mixed_factor_seas_CC(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                          nstates=55,lambda=Lambda.hat.max[,1],H=H)
    
    P.L.cc.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]))[(ncol(y)-nrow(stand.x)+56-13):length(unlist(lapply(obj$Ptt, function(x) x[1,1])))])
    P.R.cc.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[2,2]))[(ncol(y)-nrow(stand.x)+56-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
    P.theta.cc.gt <- mean(unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                       2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[(ncol(y)-nrow(stand.x)+56-13):length(unlist(lapply(obj$Ptt, function(x) x[2,2])))])
  }
}


#### Comparison of the estimates ####

P.L.biv/P.L
P.R.biv/P.R
P.theta.biv/P.theta

P.L.gt/P.L
P.R.gt/P.R
P.theta.gt/P.theta

P.L.cc.gt/P.L
P.R.cc.gt/P.R
P.theta.cc.gt/P.theta

lr.test(10260.19, 11983.11, alpha = 0.01, df = 4)


#### Save estimation results ####

if (seas_fact==F){
  
  if (CCGT == FALSE){
    if (avg==F && fixed==T && load.w==F){ #if TRUE, then append results
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.biv[1:4]),exp(par.biv[5:8])/(1-0.21^2),exp(par.biv[9:10]),exp(par.biv[12]),tanh(par.biv[11])),P.L.biv/P.L,P.R.biv/P.R,P.theta.biv/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_CC_fixed.xlsx")
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.gt[1:4]),exp(par.gt[5:8])/(1-0.21^2),tanh(par.gt[9])),P.L.gt/P.L,P.R.gt/P.R,P.theta.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_gt_fixed.xlsx")
    } else if (avg==T && fixed==F && load.w==F) {
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.gt[1:4]),exp(par.gt[5:8])/(1-0.21^2),tanh(par.gt[9])),P.L.gt/P.L,P.R.gt/P.R,P.theta.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_gt_avg.xlsx")
    } else if (avg==F && fixed==F && load.w==T) {
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.gt[1:4]),exp(par.gt[5:8])/(1-0.21^2),tanh(par.gt[9])),P.L.gt/P.L,P.R.gt/P.R,P.theta.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_gt_load.w.xlsx")
    } else if (avg==T && fixed==F && load.w==T) {
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.gt[1:4]),exp(par.gt[5:8])/(1-0.21^2),tanh(par.gt[9])),P.L.gt/P.L,P.R.gt/P.R,P.theta.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_gt_avg_load.w.xlsx")
    } else { # all FALSE
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.biv[1:4]),exp(par.biv[5:8])/(1-0.21^2),exp(par.biv[9:10]),exp(par.biv[12]),tanh(par.biv[11])),P.L.biv/P.L,P.R.biv/P.R,P.theta.biv/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_CC.xlsx")
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.gt[1:4]),exp(par.gt[5:8])/(1-0.21^2),tanh(par.gt[9])),P.L.gt/P.L,P.R.gt/P.R,P.theta.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_gt.xlsx")
    }
  } else {
    if (avg==F && fixed==T && load.w==F){ #if TRUE, then append results
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.biv[1:4]),exp(par.biv[5:8])/(1-0.21^2),exp(par.biv[9:10]),exp(par.biv[12]),tanh(par.biv[11])),P.L.biv/P.L,P.R.biv/P.R,P.theta.biv/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_CC_fixed.xlsx")
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.cc.gt[1:4]),exp(par.cc.gt[5:8])/(1-0.21^2),exp(par.cc.gt[9:10]),exp(par.cc.gt[13]),tanh(par.cc.gt[11:12])),P.L.cc.gt/P.L,P.R.cc.gt/P.R,P.theta.cc.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_cc_gt_fixed.xlsx")
    } else if (avg==T && fixed==F && load.w==F) {
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.cc.gt[1:4]),exp(par.cc.gt[5:8])/(1-0.21^2),exp(par.cc.gt[9:10]),exp(par.cc.gt[13]),tanh(par.cc.gt[11:12])),P.L.cc.gt/P.L,P.R.cc.gt/P.R,P.theta.cc.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_cc_gt_avg.xlsx")
    } else if (avg==F && fixed==F && load.w==T) {
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.cc.gt[1:4]),exp(par.cc.gt[5:8])/(1-0.21^2),exp(par.cc.gt[9:10]),exp(par.cc.gt[13]),tanh(par.cc.gt[11:12])),P.L.cc.gt/P.L,P.R.cc.gt/P.R,P.theta.cc.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_cc_gt_load.w.xlsx")
    } else if (avg==T && fixed==F && load.w==T) {
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.cc.gt[1:4]),exp(par.cc.gt[5:8])/(1-0.21^2),exp(par.cc.gt[9:10]),exp(par.cc.gt[13]),tanh(par.cc.gt[11:12])),P.L.cc.gt/P.L,P.R.cc.gt/P.R,P.theta.cc.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_cc_gt_avg_load.w.xlsx")
    } else { # all FALSE
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.biv[1:4]),exp(par.biv[5:8])/(1-0.21^2),exp(par.biv[9:10]),exp(par.biv[12]),tanh(par.biv[11])),P.L.biv/P.L,P.R.biv/P.R,P.theta.biv/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_CC.xlsx")
      write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                                 c(c(exp(par.cc.gt[1:4]),exp(par.cc.gt[5:8])/(1-0.21^2),exp(par.cc.gt[9:10]),exp(par.cc.gt[13]),tanh(par.cc.gt[11:12])),P.L.cc.gt/P.L,P.R.cc.gt/P.R,P.theta.cc.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_cc_gt.xlsx")
    }
  }
} else if (seas_fact==T) {
    if (avg==T & fixed==F & load.w==T & CCGT==F){
    write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                              c(c(exp(par.gt[1:4]),exp(par.gt[5:8])/(1-0.21^2),exp(par.gt[9]),tanh(par.gt[10])),P.L.gt/P.L,P.R.gt/P.R,P.theta.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_gt_avg_load.w.xlsx")
  } else if (avg==T & fixed==F & load.w==T & CCGT==T){ 
    write.xlsx(qpcR:::cbind.na(c(exp(par.univ[1:4]),exp(par.univ[5:8])/(1-0.21^2)),
                              c(c(exp(par.cc.gt[1:4]),exp(par.cc.gt[5:8])/(1-0.21^2),exp(par.cc.gt[9:10]),exp(par.cc.gt[14]),exp(par.cc.gt[11]),tanh(par.cc.gt[12:13])),P.L.cc.gt/P.L,P.R.cc.gt/P.R,P.theta.cc.gt/P.theta)), "Paper/Empirical analysis R/Estimation results/2004_2017/est_cc_gt_avg_load.w.xlsx")
  }
}


#### For the prediction ####

## Univariate ##
len.m <- nrow(waves)
h <- len.m - trunc((len.m)/3) # forecast in third half of the sample 
P.L <- rep(NA,(len.m-h))
P.R <- rep(NA,(len.m-h))
P.theta <- rep(NA,(len.m-h))
par.univ <- matrix(NA, (len.m-h), 8)
R_y_est_tt <- rep(NA,(len.m-h))
R_y_est_ttm1 <- rep(NA,(len.m-h))
L_y_est_tt <- rep(NA,(len.m-h))
L_y_est_ttm1 <- rep(NA,(len.m-h))
theta_y_est_tt <- rep(NA,(len.m-h))
theta_y_est_ttm1 <- rep(NA,(len.m-h))
slope <- lapply(seq_len(len.m-1-h), function(X) rep(NA,len.m))

for (r in h:(len.m-1)){ # r is the number of observed values
  y <- matrix(0,5,(r+1))
  for (j in 1:5){
    y[j,1:r] <- waves[1:r,j]
    y[j,(r+1)] <- NA
  }
  
  se <- as.data.frame(LFS[1:r,c(2,44,86,128,170)])
  k <- rbind(se,se[nrow(se),])
  
  objopt <-  optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                         log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2))),
                   KF_slopes_univ,y=y[1:5,],k=k,delta=0.21,opti=T,outofsample=T,parP10=1000000000000,nstates=30,  hessian=F, method="L-BFGS-B" )
  par.univ[(r-h+1),] <- c(objopt$par)
  
  obj <- KF_slopes_univ(par=objopt$par,y=y[1:5,],k=k,delta=0.21,opti=F,outofsample=T,parP10=1000000000000,nstates=30)
  
  P.L[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]))[length(unlist(lapply(obj$Ptt, function(x) x[1,1])))]
  P.R[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[2,2]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
  P.theta[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                    2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
  #ts.plot(unlist(lapply(obj$Ptt, function(x) x[1,1]))[31:ncol(y)])
  
  #estimated R
  R_y_est_tt[r-h+1] <- obj$xtt[2,ncol(y)] 
  R_y_est_ttm1[r-h+1] <- obj$xttm1[2,ncol(y)]
  
  #estimated L
  L_y_est_tt[r-h+1] <- obj$xtt[1,ncol(y)] 
  L_y_est_ttm1[r-h+1] <- obj$xttm1[1,ncol(y)]
  
  #estimated seasonal components
  season_y_est_tt <- obj$xtt[3,ncol(y)]+obj$xtt[5,ncol(y)]+obj$xtt[7,ncol(y)]+obj$xtt[9,ncol(y)]+obj$xtt[11,ncol(y)]+obj$xtt[13,ncol(y)]
  season_y_est_ttm1 <- obj$xttm1[3,ncol(y)]+obj$xttm1[5,ncol(y)]+obj$xttm1[7,ncol(y)]+obj$xttm1[9,ncol(y)]+obj$xttm1[11,ncol(y)]+obj$xttm1[13,ncol(y)]
  
  #estimated theta
  theta_y_est_tt[r-h+1] <- obj$xtt[1,ncol(y)] + season_y_est_tt 
  theta_y_est_ttm1[r-h+1] <- obj$xttm1[1,ncol(y)] + season_y_est_ttm1 
  
  slope[[r-h+1]] <- obj$xtt[2,(31-13):length(obj$xtt[2,])]
  
  print(paste(P.L[r-h+1],P.R[r-h+1],P.theta[r-h+1]))
}

now_df_univ <- cbind(P.L,P.R,P.theta,R_y_est_tt,R_y_est_ttm1,L_y_est_tt,L_y_est_ttm1,theta_y_est_tt,theta_y_est_ttm1)
if (fixed == T) {
  saveRDS(now_df_univ, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_univ_fixed.csv"))
  saveRDS(par.univ, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.univ_fixed.csv"))
} else {
  saveRDS(now_df_univ, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_univ.csv"))
  saveRDS(par.univ, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.univ.csv"))
}


## Bivariate with CC ##
len.m <- nrow(waves)
h <- len.m - trunc((len.m)/3)  # forecast in third half of the sample 
P.L.biv <- rep(NA,(len.m-h))
P.R.biv <- rep(NA,(len.m-h))
P.theta.biv <- rep(NA,(len.m-h))
par.biv <- matrix(NA, (len.m-h), 12)
R_y_est_tt.biv <- rep(NA,(len.m-h))
R_y_est_ttm1.biv <- rep(NA,(len.m-h))
L_y_est_tt.biv <- rep(NA,(len.m-h))
L_y_est_ttm1.biv <- rep(NA,(len.m-h))
theta_y_est_tt.biv <- rep(NA,(len.m-h))
theta_y_est_ttm1.biv <- rep(NA,(len.m-h))

for (r in h:(len.m-1)){ # r is the number of observed values
  y <- matrix(0,6,(r+1))
  for (j in 1:5){
    y[j,1:r] <- waves[1:r,j]
    y[j,(r+1)] <- NA
  }
  y[6,1:r] <- CC[1:r,1]
  y[6,(r+1)] <- NA
  
  se <- as.data.frame(LFS[1:r,c(2,44,86,128,170)])
  k <- rbind(se,se[nrow(se),])
  
  objopt <-  optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                         log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),
                         log(3000),log(0.02),0,log(1000)),
                   KF_slopes,y=y,k=k,delta=0.21,opti=T,outofsample=T,parP10=1000000000000,nstates=43,  hessian=F, method="L-BFGS-B" )
  par.biv[(r-h+1),] <- objopt$par
  
  obj <- KF_slopes(par=objopt$par,y=y,k=k,delta=0.21,opti=F,outofsample=T,parP10=1000000000000,nstates=43)
  
  P.L.biv[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]))[length(unlist(lapply(obj$Ptt, function(x) x[1,1])))]
  P.R.biv[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[2,2]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
  P.theta.biv[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                        2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
  #ts.plot(unlist(lapply(obj$Ptt, function(x) x[32,32]))[31:ncol(y)])
  
  #estimated R
  R_y_est_tt.biv[r-h+1] <- obj$xtt[2,ncol(y)] 
  R_y_est_ttm1.biv[r-h+1] <- obj$xttm1[2,ncol(y)]
  
  #estimated L
  L_y_est_tt.biv[r-h+1] <- obj$xtt[1,ncol(y)] 
  L_y_est_ttm1.biv[r-h+1] <- obj$xttm1[1,ncol(y)]
  
  #estimated seasonal components
  season_y_est_tt.biv <- obj$xtt[3,ncol(y)]+obj$xtt[5,ncol(y)]+obj$xtt[7,ncol(y)]+obj$xtt[9,ncol(y)]+obj$xtt[11,ncol(y)]+obj$xtt[13,ncol(y)]
  season_y_est_ttm1.biv <- obj$xttm1[3,ncol(y)]+obj$xttm1[5,ncol(y)]+obj$xttm1[7,ncol(y)]+obj$xttm1[9,ncol(y)]+obj$xttm1[11,ncol(y)]+obj$xttm1[13,ncol(y)]
  
  #estimated theta
  theta_y_est_tt.biv[r-h+1] <- obj$xtt[1,ncol(y)] + season_y_est_tt.biv 
  theta_y_est_ttm1.biv[r-h+1] <- obj$xttm1[1,ncol(y)] + season_y_est_ttm1.biv
  
  print(paste(P.L.biv[r-h+1],P.R.biv[r-h+1],P.theta.biv[r-h+1]))
}

now_df_biv <- cbind(P.L.biv,P.R.biv,P.theta.biv,R_y_est_tt.biv,R_y_est_ttm1.biv,L_y_est_tt.biv,L_y_est_ttm1.biv,theta_y_est_tt.biv,theta_y_est_ttm1.biv)
if (fixed == T) {
  saveRDS(now_df_biv, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_biv_fixed.csv"))
  saveRDS(par.biv, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.biv_fixed.csv"))
} else {
  saveRDS(now_df_biv, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_biv.csv"))
  saveRDS(par.biv, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.biv.csv"))
}


## Bivariate with Google trends ##
if (monthlyGT==T && CCGT==F){
  len.m <- nrow(waves)
  h <- len.m - trunc((len.m)/3)  # forecast in third half of the sample 
  P.L.gt <- rep(NA,(len.m-h))
  P.R.gt <- rep(NA,(len.m-h))
  P.theta.gt <- rep(NA,(len.m-h))
  if (seas_fact==F){
    par.gt <- matrix(NA, (len.m-h), 9)
  } else if (seas_fact==T && fixed_sigma_u==T) {
    par.gt <- matrix(NA, (len.m-h), 10)
  } else if (seas_fact==T && fixed_sigma_u==F) {
    par.gt <- matrix(NA, (len.m-h), 11)
  }
  R_y_est_tt.gt <- rep(NA,(len.m-h))
  R_y_est_ttm1.gt <- rep(NA,(len.m-h))
  L_y_est_tt.gt <- rep(NA,(len.m-h))
  L_y_est_ttm1.gt <- rep(NA,(len.m-h))
  theta_y_est_tt.gt <- rep(NA,(len.m-h))
  theta_y_est_ttm1.gt <- rep(NA,(len.m-h))
  sel.var <- lapply(seq_len(len.m-h), function(X) rep(NA,ncol(dataset)))
  
  for (r in h:(len.m-1)){ # r is the number of observed values
    
    len.m <- nrow(waves)
    agg.gt <- x[1:(r+2),]
    
    if (EN ==T){ 
      agg.gt <- hhh[1:(r+2),]
      # Elastic Net
      diff.x <- diff(agg.gt[(32-13):nrow(agg.gt),],1)
      diff.y <- diff(slope[[r-h+1]],1)
      for (i in 1:ncol(diff.x)){
        diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
      }
      diff.y <- (diff.y-mean(diff.y))/sd(diff.y)
      alpha_grid <- seq(0.05, 0.95, 0.01) # grid values for alpha
      AIC <- rep(NA,length(alpha_grid))
      tune_lambda <- rep(NA,length(alpha_grid))
      nselect_var <- rep(NA,length(alpha_grid))
      for (i in 1:length(alpha_grid)) {
        en <- ic.glmnet(diff.x, diff.y, crit = "aic", alpha=alpha_grid[i], lambda=seq(0.1, 0, -0.001)) #alpha=1 lasso, alpha=0 ridge
        # plot(en$glmnet,"lambda")
        # plot(en)
        tune_lambda[i] <- en$lambda
        AIC[i] <- en$ic[2]
        nselect_var[i] <-en$nvar-1
        # coef(en)
        # en$glmnet$lambda
      }
      nselect_var[which.min(AIC)]
      fit <- glmnet(diff.x, diff.y, lambda = tune_lambda[which.min(AIC)], alpha = alpha_grid[which.min(AIC)])
      coef(fit)
      sel <- fit$beta@i[-1]+1
      agg.gt <- agg.gt[,sel]
      sel.var[[r-h+1]] <- sel
    }
    
    # Differenciate data
    diff.x <- diff(agg.gt,1) 
    diff.x.w <- diff(x.w,1)
    
    if(EN==T){
      diff.x.w <- diff(x.w[,sel],1)
    }
    
    if (load.w==T){
      diff.x <- diff.x.w
    }
    
    # Standardize the monthly Google trends
    for (i in 1:ncol(diff.x)){
      diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
    }
    diff.x <- as.matrix(diff.x)
    kmax <- 8 # max number of factors
    cov <- diff.x%*%t(diff.x)
    len.m <- nrow(diff.x)
    nvar <- ncol(diff.x)
    kmax.eigenvalues <- eigen(cov)$values[1:kmax]
    kmax.eigenvectors <- eigen(cov)$vectors[,1:kmax] 
    F.hat.max <- kmax.eigenvectors*sqrt(len.m-1) # facotrs
    t(F.hat.max)%*%F.hat.max
    Lambda.hat.max <- t(diff.x)%*%F.hat.max/(len.m-1) # loadings
    Comps.hat.max <- F.hat.max%*%t(Lambda.hat.max) # common components
    e.max <- diff.x-Comps.hat.max # estimated hydiosincratic components (hc)
    
    # We need the kmax variance of the hc for the information criteria (IC)
    sigma.hat.maxj <- c(1:nvar)
    Vj <- sigma.hat.maxj
    for (j in 1:nvar){
      sigma.hat.maxj[j] <- 1/len.m*t(e.max[,j])%*%(e.max[,j])
    }
    sigma.hat.max <- 1/nvar*sum(sigma.hat.maxj)
    
    # Define vectors to store the values for the IC (Bai and Ng)
    PC1 <- c(1:kmax)
    PC2 <- PC1
    PC3 <- PC1
    IC1 <- PC1
    IC2 <- PC1
    IC3 <- PC1
    
    # Loop to choose the number of factors based on IC
    for (i in 1:kmax){
      k.eigenvalues <- eigen(cov)$values[1:i]
      k.eigenvectors <- eigen(cov)$vectors[,1:i]
      F.hat <- k.eigenvectors*sqrt(len.m-1)
      Lambda.hat <- t(diff.x)%*%F.hat/(len.m-1)
      Comps.hat <- F.hat%*%t(Lambda.hat)
      e <- diff.x-Comps.hat
      # Try with differenced data
      # e <- diff(e)
      for (j in 1:nvar){
        Vj[j] <- 1/len.m*t(e[,j])%*%(e[,j])
      }
      V <- 1/nvar*sum(Vj)
      # Information Criteria for stationary factors
      PC1[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
      PC2[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
      PC3[i] <- V + i*sigma.hat.max*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
      IC1[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
      IC2[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
      IC3[i] <- log(V) + i*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
    }
    
    # Information Criteria
    k.PC1 <- which.min(PC1)
    k.PC2 <- which.min(PC2)
    k.PC3 <- which.min(PC3)
    k.IC1 <- which.min(IC1)
    k.IC2 <- which.min(IC2)
    k.IC3 <- which.min(IC3)
    
    # Double-check with a built-in function
    a <- prcomp(diff.x, rank.=kmax)
    screeplot(a)
    plot(sqrt(kmax.eigenvalues))
    
    #Non-stationary standardized Google trends
    stand.x <- diff.x
    for (i in 1:ncol(diff.x)){
      stand.x[,i] <- cumsum(diff.x[,i])
    }
    stand.x <- as.matrix(stand.x)
    
    # Non-stationary monthly factors
    F.ns.hat <- F.hat.max
    for (i in 1:kmax){
      F.ns.hat[,i] <- cumsum(F.hat.max[,i])
      plot(F.ns.hat[,i], type="l")
    }
    # covariance matrix of the idiosyncratic components:
    idio <- stand.x - F.ns.hat[,1]%*%t(Lambda.hat.max[,1])
    H <- 1/(len.m-1)*t(idio)%*%idio
    
    merged <- matrix(0,(r+1),5)
    for (j in 1:5){
      merged[1:r,j] <- waves[1:r,j]
      merged[(r+1),j] <- NA
    }
    y <- t(cbind(merged,stand.x))
    
    se <- as.data.frame(LFS[1:r,c(2,44,86,128,170)])
    k <- rbind(se,se[nrow(se),])
    
    if (seas_fact==F){
      
      objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                            log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),0), 
                      KF_slopes_mixed_factor,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                      parP10=1000000000000,nstates=31,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
      par.gt[(r-h+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                    nstates=31,lambda=Lambda.hat.max[,1],H=H)
    } else if (seas_fact==T && fixed_sigma_u==T){
      objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                            log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(0.02),0), 
                      KF_slopes_mixed_factor_seas,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                      parP10=1000000000000,nstates=42,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
      par.gt <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_seas(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                         nstates=42,lambda=Lambda.hat.max[,1],H=H)
      
    } else if (seas_fact==T && fixed_sigma_u==F) {
      objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                            log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(1),log(1),0), 
                      KF_slopes_mixed_factor_seas,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                      parP10=1000000000000,nstates=42,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
      par.gt <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_seas(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                         nstates=42,lambda=Lambda.hat.max[,1],H=H)
    }
    
    P.L.gt[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]))[length(unlist(lapply(obj$Ptt, function(x) x[1,1])))]
    P.R.gt[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[2,2]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
    P.theta.gt[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                              2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
    #ts.plot(unlist(lapply(obj$Ptt, function(x) x[31,31]))[31:ncol(y)])
    
    #estimated R
    R_y_est_tt.gt[r-h+1] <- obj$xtt[2,ncol(y)] 
    R_y_est_ttm1.gt[r-h+1] <- obj$xttm1[2,ncol(y)]
    
    #estimated L
    L_y_est_tt.gt[r-h+1] <- obj$xtt[1,ncol(y)] 
    L_y_est_ttm1.gt[r-h+1] <- obj$xttm1[1,ncol(y)]
    
    #estimated seasonal components
    season_y_est_tt.gt <- obj$xtt[3,ncol(y)]+obj$xtt[5,ncol(y)]+obj$xtt[7,ncol(y)]+obj$xtt[9,ncol(y)]+obj$xtt[11,ncol(y)]+obj$xtt[13,ncol(y)]
    season_y_est_ttm1.gt <- obj$xttm1[3,ncol(y)]+obj$xttm1[5,ncol(y)]+obj$xttm1[7,ncol(y)]+obj$xttm1[9,ncol(y)]+obj$xttm1[11,ncol(y)]+obj$xttm1[13,ncol(y)]
    
    #estimated theta
    theta_y_est_tt.gt[r-h+1] <- obj$xtt[1,ncol(y)] + season_y_est_tt.gt 
    theta_y_est_ttm1.gt[r-h+1] <- obj$xttm1[1,ncol(y)] + season_y_est_ttm1.gt
    
    print(paste(P.L.gt[r-h+1],P.R.gt[r-h+1],P.theta.gt[r-h+1]))
  }
  
  now_df_gt <- cbind(P.L.gt,P.R.gt,P.theta.gt,R_y_est_tt.gt,R_y_est_ttm1.gt,L_y_est_tt.gt,L_y_est_ttm1.gt,theta_y_est_tt.gt,theta_y_est_ttm1.gt)
  if (avg == F && fixed == T && load.w == F) {
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt.csv"))
  } else if (avg == T && fixed == F && load.w == F) {
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt_avg.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt_avg.csv"))
  } else if (avg == F && fixed == F && load.w == T) {
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt_load.w.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt_load.w.csv"))
  } else if (avg == T && fixed == F && load.w == T) {
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt_avg_load.w.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt_avg_load.w.csv"))
  } else { # all FALSE
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt.csv"))
  }
  
  if (EN==T){
    saveRDS(sel.var, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/sel.var.gt.csv"))
  }
  
  
} else if (monthlyGT==F && CCGT==F){
  len.m <- nrow(x)
  len <- nrow(x.w)
  dataset <- as.data.frame(dataset)
  week <- as.factor(dataset[1:nrow(dataset),1])
  week <- as.Date(week,"%Y-%m-%d")
  dataset_month <- as.data.frame(dataset_month)
  month <- as.factor(dataset_month[1:nrow(dataset_month),1])
  month <- as.Date(month,"%Y-%m-%d")
  h <- len.m-1 - trunc((len.m-1)/3) # forecast in third part of the sample 
  for (j in 1:len){
    if (strftime(month[h+2], "%Y-%m") == strftime(week[j], "%Y-%m")) break
    h.week <- j+1
  }
  h.month <- h+2
  
  P.L.gt <- rep(NA,(len-h.week+1))
  P.R.gt <- rep(NA,(len-h.week+1))
  P.theta.gt <- rep(NA,(len-h.week+1))
  if (seas_fact==F){
    par.gt <- matrix(NA, (len-h.week+1), 9)
  } else if (seas_fact==T && fixed_sigma_u==T) {
    par.gt <- matrix(NA, (len-h.week+1), 10)
  } else if (seas_fact==T && fixed_sigma_u==F) {
    par.gt <- matrix(NA, (len-h.week+1), 11)
  }
  R_y_est_tt.gt <- rep(NA,(len-h.week+1))
  R_y_est_ttm1.gt <- rep(NA,(len-h.week+1))
  L_y_est_tt.gt <- rep(NA,(len-h.week+1))
  L_y_est_ttm1.gt <- rep(NA,(len-h.week+1))
  theta_y_est_tt.gt <- rep(NA,(len-h.week+1))
  theta_y_est_ttm1.gt <- rep(NA,(len-h.week+1))
  sel.var <- lapply(seq_len(len-h.week+1), function(X) rep(NA,ncol(dataset)))
  
  for (r in h.week:len){
    
    len.m <- nrow(x)
    len <- nrow(x.w)
    dataset <- as.data.frame(dataset)
    week <- as.factor(dataset[1:nrow(dataset),1])
    week <- as.Date(week,"%Y-%m-%d")
    dataset_month <- as.data.frame(dataset_month)
    month <- as.factor(dataset_month[1:nrow(dataset_month),1])
    month <- as.Date(month,"%Y-%m-%d")
    
    for (j in 1:len.m){
      if (strftime(month[j], "%Y-%m") == strftime(week[r], "%Y-%m")) break
      h <- j+1
    }
    
    # Aggregate weekly Google trends
    dataset <- as.data.frame(dataset)
    week <- as.factor(dataset[1:r,1])
    week <- as.Date(week,"%Y-%m-%d")
    agg.gt <- matrix(NA,nrow(x[1:h,]),ncol(x))
    agg.gt2 <- matrix(NA,nrow(x[1:h,]),ncol(x))
    agg.gt3 <- matrix(NA,nrow(x[1:h,]),ncol(x))
    agg.gt4 <- matrix(NA,nrow(x[1:h,]),ncol(x))
    agg.gt5 <- matrix(NA,nrow(x[1:h,]),ncol(x))
    
    if (avg==T){
      
      # Aggregate as flow variables
      for (i in 1:ncol(x.w)){
        agg.gt[,i] <- flow_avg(x=x.w[1:r,i], k=1, per=week) # flow variables
        agg.gt2 <- flow_avg(x=x.w[1:r,i], k=2, per=week) # flow variables
        agg.gt3 <- flow_avg(x=x.w[1:r,i], k=3, per=week) # flow variables
        agg.gt4 <- flow_avg(x=x.w[1:r,i], k=4, per=week) # flow variables
        agg.gt5 <- flow_avg(x=x.w[1:r,i], k=5, per=week)
        two.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[2]})$x))
        three.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[3]})$x))
        four.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[4]})$x))
        five.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[5]})$x))
        agg.gt[two.obs,i] <- agg.gt2[two.obs]
        agg.gt[three.obs,i] <- agg.gt3[three.obs]
        agg.gt[four.obs,i] <- agg.gt4[four.obs]
        agg.gt[five.obs,i] <- agg.gt5[five.obs]
      }
    } else {
      
      # Aggregate as flow variables
      for (i in 1:ncol(x.w)){
        agg.gt[,i] <- flow_avg(x=x.w[1:r,i], k=1, per=week) # flow variables
        agg.gt2 <- flow(x=x.w[1:r,i], k=2, per=week) # flow variables
        agg.gt3 <- flow(x=x.w[1:r,i], k=3, per=week) # flow variables
        agg.gt4 <- flow(x=x.w[1:r,i], k=4, per=week) # flow variables
        agg.gt5 <- flow(x=x.w[1:r,i], k=5, per=week)
        two.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[2]})$x))
        three.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[3]})$x))
        four.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[4]})$x))
        five.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[5]})$x))
        agg.gt[two.obs,i] <- agg.gt2[two.obs]
        agg.gt[three.obs,i] <- agg.gt3[three.obs]
        agg.gt[four.obs,i] <- agg.gt4[four.obs]
        agg.gt[five.obs,i] <- agg.gt5[five.obs]
      }
      
      # Rescale according to largest value on the entire dataset
      for (i in 1:ncol(agg.gt)){
        if (max(agg.gt[,i]) != 0){
          agg.gt[,i] <- agg.gt[,i]*100/max(agg.gt[,i])
        } else {
          agg.gt[,i] <- agg.gt[,i]
        }
      }
      
    }
    
    if (EN ==T){ 
      # Elastic Net
      diff.x <- diff(agg.gt[(32-13):nrow(agg.gt),],1)
      diff.y <- diff(slope[[h-h.month+1]],1)
      for (i in 1:ncol(diff.x)){
        diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
      }
      diff.y <- (diff.y-mean(diff.y))/sd(diff.y)
      alpha_grid <- seq(0.05, 0.95, 0.01) # grid values for alpha
      AIC <- rep(NA,length(alpha_grid))
      tune_lambda <- rep(NA,length(alpha_grid))
      nselect_var <- rep(NA,length(alpha_grid))
      for (i in 1:length(alpha_grid)) {
        en <- ic.glmnet(diff.x, diff.y, crit = "aic", alpha=alpha_grid[i], lambda=seq(0.1, 0, -0.001)) #alpha=1 lasso, alpha=0 ridge
        # plot(en$glmnet,"lambda")
        # plot(en)
        tune_lambda[i] <- en$lambda
        AIC[i] <- en$ic[2]
        nselect_var[i] <-en$nvar-1
        # coef(en)
        # en$glmnet$lambda
      }
      nselect_var[which.min(AIC)]
      fit <- glmnet(diff.x, diff.y, lambda = tune_lambda[which.min(AIC)], alpha = alpha_grid[which.min(AIC)])
      coef(fit)
      sel <- fit$beta@i[-1]+1
      agg.gt <- agg.gt[,sel]
      sel.var[[r-h.week+1]] <- sel
    }
    
    # Differenciate data
    diff.x <- diff(agg.gt,1)
    diff.x.w <- diff(x.w[1:r,],1)
    
    if (EN==T){
      diff.x.w <- diff(x.w[1:r,sel],1)
    }
    
    if (load.w==T){
      diff.x <- diff.x.w
    }
    
    # Standardize the monthly Google trends
    for (i in 1:ncol(diff.x)){
      diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
    }
    diff.x <- as.matrix(diff.x)
    kmax <- 4 # max number of factors
    cov <- diff.x%*%t(diff.x)
    len.m <- nrow(diff.x)
    nvar <- ncol(diff.x)
    kmax.eigenvalues <- eigen(cov)$values[1:kmax]
    kmax.eigenvectors <- eigen(cov)$vectors[,1:kmax] 
    F.hat.max <- kmax.eigenvectors*sqrt(len.m-1) # facotrs
    t(F.hat.max)%*%F.hat.max
    Lambda.hat.max <- t(diff.x)%*%F.hat.max/(len.m-1) # loadings
    Comps.hat.max <- F.hat.max%*%t(Lambda.hat.max) # common components
    e.max <- diff.x-Comps.hat.max # estimated hydiosincratic components (hc)
    
    # We need the kmax variance of the hc for the information criteria (IC)
    sigma.hat.maxj <- c(1:nvar)
    Vj <- sigma.hat.maxj
    for (j in 1:nvar){
      sigma.hat.maxj[j] <- 1/len.m*t(e.max[,j])%*%(e.max[,j])
    }
    sigma.hat.max <- 1/nvar*sum(sigma.hat.maxj)
    
    # Define vectors to store the values for the IC (Bai and Ng)
    PC1 <- c(1:kmax)
    PC2 <- PC1
    PC3 <- PC1
    IC1 <- PC1
    IC2 <- PC1
    IC3 <- PC1
    
    # Loop to choose the number of factors based on IC
    for (i in 1:kmax){
      k.eigenvalues <- eigen(cov)$values[1:i]
      k.eigenvectors <- eigen(cov)$vectors[,1:i]
      F.hat <- k.eigenvectors*sqrt(len.m-1)
      Lambda.hat <- t(diff.x)%*%F.hat/(len.m-1)
      Comps.hat <- F.hat%*%t(Lambda.hat)
      e <- diff.x-Comps.hat
      # Try with differenced data
      # e <- diff(e)
      for (j in 1:nvar){
        Vj[j] <- 1/len.m*t(e[,j])%*%(e[,j])
      }
      V <- 1/nvar*sum(Vj)
      # Information Criteria for stationary factors
      PC1[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
      PC2[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
      PC3[i] <- V + i*sigma.hat.max*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
      IC1[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
      IC2[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
      IC3[i] <- log(V) + i*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
    }
    
    # Information Criteria
    k.PC1 <- which.min(PC1)
    k.PC2 <- which.min(PC2)
    k.PC3 <- which.min(PC3)
    k.IC1 <- which.min(IC1)
    k.IC2 <- which.min(IC2)
    k.IC3 <- which.min(IC3)
    
    #Non-stationary standardized Google trends
    stand.x <- diff.x
    for (i in 1:ncol(diff.x)){
      stand.x[,i] <- cumsum(diff.x[,i])
    }
    stand.x <- as.matrix(stand.x)
    
    # Non-stationary monthly factors
    F.ns.hat <- F.hat.max
    for (i in 1:kmax){
      F.ns.hat[,i] <- cumsum(F.hat.max[,i])
      # plot(F.ns.hat[,i], type="l")
    }
    # covariance matrix of the idiosyncratic components:
    idio <- stand.x - F.ns.hat[,1]%*%t(Lambda.hat.max[,1])
    H <- 1/(len.m-1)*t(idio)%*%idio
    
    if (load.w==T){
      diff.x <- diff(agg.gt,1)
      # Standardize the monthly Google trends
      for (i in 1:ncol(diff.x)){
        diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
      }
      diff.x <- as.matrix(diff.x)
      stand.x <- diff.x
      for (i in 1:ncol(diff.x)){
        stand.x[,i] <- cumsum(diff.x[,i])
      }
      stand.x <- as.matrix(stand.x)
    }
    
    # Monthly google trends
    merged <- matrix(0,(h-1),5)
    for (j in 1:5){
      merged[1:(h-1),j] <- waves[1:(h-1),j]
      merged[(h-1),j] <- NA
    }
    y <- t(cbind(merged,rbind(matrix(NA,(nrow(merged)-nrow(stand.x)),ncol(stand.x)),stand.x)))
    
    # Kalman filter estimation
    se <- as.data.frame(LFS[1:(h-1),c(2,44,86,128,170)])
    k <- rbind(se,se[nrow(se),])
    
    if (seas_fact==F){
      possibleError <- tryCatch(objopt <-  optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                                 log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),0),
                       KF_slopes_mixed_factor,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                       parP10=1000000000000,nstates=31,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B" )
      , error = function(e) e)
      if(inherits(possibleError, "error")) next
      par.gt[(r-h.week+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                  nstates=31,lambda=Lambda.hat.max[,1],H=H)
    } else if (seas_fact==T && fixed_sigma_u==T) {
      possibleError <- tryCatch(objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                                                      log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(0.02),0), 
                                                KF_slopes_mixed_factor_seas,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                                                parP10=1000000000000,nstates=42,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
                                , error = function(e) e)
      if(inherits(possibleError, "error")) next
      par.gt[(r-h.week+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_seas(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                         nstates=42,lambda=Lambda.hat.max[,1],H=H)
      
    } else if (seas_fact==T && fixed_sigma_u==F) {
      possibleError <- tryCatch(objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                            log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(1),log(1),0), 
                      KF_slopes_mixed_factor_seas,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                      parP10=1000000000000,nstates=42,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
                      , error = function(e) e)
      if(inherits(possibleError, "error")) next
      par.gt[(r-h.week+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_seas(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                         nstates=42,lambda=Lambda.hat.max[,1],H=H)
      }
    
    P.L.gt[r-h.week+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]))[length(unlist(lapply(obj$Ptt, function(x) x[1,1])))]
    P.R.gt[r-h.week+1] <- unlist(lapply(obj$Ptt, function(x) x[2,2]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
    P.theta.gt[r-h.week+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                              2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
    #ts.plot(unlist(lapply(obj$Ptt, function(x) x[31,31]))[31:ncol(y)])
    
    #estimated R
    R_y_est_tt.gt[r-h.week+1] <- obj$xtt[2,ncol(y)] 
    R_y_est_ttm1.gt[r-h.week+1] <- obj$xttm1[2,ncol(y)]
    
    #estimated L
    L_y_est_tt.gt[r-h.week+1] <- obj$xtt[1,ncol(y)] 
    L_y_est_ttm1.gt[r-h.week+1] <- obj$xttm1[1,ncol(y)]
    
    #estimated seasonal components
    season_y_est_tt.gt <- obj$xtt[3,ncol(y)]+obj$xtt[5,ncol(y)]+obj$xtt[7,ncol(y)]+obj$xtt[9,ncol(y)]+obj$xtt[11,ncol(y)]+obj$xtt[13,ncol(y)]
    season_y_est_ttm1.gt <- obj$xttm1[3,ncol(y)]+obj$xttm1[5,ncol(y)]+obj$xttm1[7,ncol(y)]+obj$xttm1[9,ncol(y)]+obj$xttm1[11,ncol(y)]+obj$xttm1[13,ncol(y)]
    
    #estimated theta
    theta_y_est_tt.gt[r-h.week+1] <- obj$xtt[1,ncol(y)] + season_y_est_tt.gt 
    theta_y_est_ttm1.gt[r-h.week+1] <- obj$xttm1[1,ncol(y)] + season_y_est_ttm1.gt
    
    print(paste(P.L.gt[r-h.week+1],P.R.gt[r-h.week+1],P.theta.gt[r-h.week+1]))
  }
  
  now_df_gt <- cbind(P.L.gt,P.R.gt,P.theta.gt,R_y_est_tt.gt,R_y_est_ttm1.gt,L_y_est_tt.gt,L_y_est_ttm1.gt,theta_y_est_tt.gt,theta_y_est_ttm1.gt)
  if (avg == F && fixed == T && load.w == F) {
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt.csv"))
  } else if (avg == T && fixed == F && load.w == F) {
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt_avg.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt_avg.csv"))
  } else if (avg == F && fixed == F && load.w == T) {
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt_load.w.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt_load.w.csv"))
  } else if (avg == T && fixed == F && load.w == T) {
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt_avg_load.w.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt_avg_load.w.csv"))
  } else { # all FALSE
    saveRDS(now_df_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt.csv"))
    saveRDS(par.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt.csv"))
  }
  
  if (EN==T){
    saveRDS(sel.var, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/sel.var.gt.csv"))
  }
}


## Bivariate with Claimant Counts and Google trends ##
if (monthlyGT==T && CCGT==T){
  len.m <- nrow(waves)
  h <- len.m - trunc((len.m)/3)  # forecast in third half of the sample 
  P.L.cc.gt <- rep(NA,(len.m-h))
  P.R.cc.gt <- rep(NA,(len.m-h))
  P.theta.cc.gt <- rep(NA,(len.m-h))
  if (seas_fact==F){
    par.cc.gt <- matrix(NA, (len.m-h), 13)
  } else if (seas_fact==T && fixed_sigma_u==T) {
    par.cc.gt <- matrix(NA, (len.m-h), 14)
  } else if (seas_fact==T && fixed_sigma_u==F) {
    par.cc.gt <- matrix(NA, (len.m-h), 15)
  }
  R_y_est_tt.cc.gt <- rep(NA,(len.m-h))
  R_y_est_ttm1.cc.gt <- rep(NA,(len.m-h))
  L_y_est_tt.cc.gt <- rep(NA,(len.m-h))
  L_y_est_ttm1.cc.gt <- rep(NA,(len.m-h))
  theta_y_est_tt.cc.gt <- rep(NA,(len.m-h))
  theta_y_est_ttm1.cc.gt <- rep(NA,(len.m-h))
  sel.var <- lapply(seq_len(len.m-h), function(X) rep(NA,ncol(dataset)))
  
  for (r in h:(len.m-1)){ # r is the number of observed values
    
    len.m <- nrow(waves)
    agg.gt <- x[1:(r+2),]
    
    if (EN ==T){ 
      agg.gt <- hhh[1:(r+2),]
      # Elastic Net
      diff.x <- diff(agg.gt[(32-13):nrow(agg.gt),],1)
      diff.y <- diff(slope[[r-h+1]],1)
      for (i in 1:ncol(diff.x)){
        diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
      }
      diff.y <- (diff.y-mean(diff.y))/sd(diff.y)
      alpha_grid <- seq(0.05, 0.95, 0.01) # grid values for alpha
      AIC <- rep(NA,length(alpha_grid))
      tune_lambda <- rep(NA,length(alpha_grid))
      nselect_var <- rep(NA,length(alpha_grid))
      for (i in 1:length(alpha_grid)) {
        en <- ic.glmnet(diff.x, diff.y, crit = "aic", alpha=alpha_grid[i], lambda=seq(0.1, 0, -0.001)) #alpha=1 lasso, alpha=0 ridge
        # plot(en$glmnet,"lambda")
        # plot(en)
        tune_lambda[i] <- en$lambda
        AIC[i] <- en$ic[2]
        nselect_var[i] <-en$nvar-1
        # coef(en)
        # en$glmnet$lambda
      }
      nselect_var[which.min(AIC)]
      fit <- glmnet(diff.x, diff.y, lambda = tune_lambda[which.min(AIC)], alpha = alpha_grid[which.min(AIC)])
      coef(fit)
      sel <- fit$beta@i[-1]+1
      agg.gt <- agg.gt[,sel]
      sel.var[[r-h+1]] <- sel
    }
    
    # Differenciate data
    diff.x <- diff(agg.gt,1) 
    diff.x.w <- diff(x.w,1)
    
    if(EN==T){
      diff.x.w <- diff(x.w[,sel],1)
    }
    
    if (load.w==T){
      diff.x <- diff.x.w
    }
    
    # Standardize the monthly Google trends
    for (i in 1:ncol(diff.x)){
      diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
    }
    diff.x <- as.matrix(diff.x)
    kmax <- 8 # max number of factors
    cov <- diff.x%*%t(diff.x)
    len.m <- nrow(diff.x)
    nvar <- ncol(diff.x)
    kmax.eigenvalues <- eigen(cov)$values[1:kmax]
    kmax.eigenvectors <- eigen(cov)$vectors[,1:kmax] 
    F.hat.max <- kmax.eigenvectors*sqrt(len.m-1) # facotrs
    t(F.hat.max)%*%F.hat.max
    Lambda.hat.max <- t(diff.x)%*%F.hat.max/(len.m-1) # loadings
    Comps.hat.max <- F.hat.max%*%t(Lambda.hat.max) # common components
    e.max <- diff.x-Comps.hat.max # estimated hydiosincratic components (hc)
    
    # We need the kmax variance of the hc for the information criteria (IC)
    sigma.hat.maxj <- c(1:nvar)
    Vj <- sigma.hat.maxj
    for (j in 1:nvar){
      sigma.hat.maxj[j] <- 1/len.m*t(e.max[,j])%*%(e.max[,j])
    }
    sigma.hat.max <- 1/nvar*sum(sigma.hat.maxj)
    
    # Define vectors to store the values for the IC (Bai and Ng)
    PC1 <- c(1:kmax)
    PC2 <- PC1
    PC3 <- PC1
    IC1 <- PC1
    IC2 <- PC1
    IC3 <- PC1
    
    # Loop to choose the number of factors based on IC
    for (i in 1:kmax){
      k.eigenvalues <- eigen(cov)$values[1:i]
      k.eigenvectors <- eigen(cov)$vectors[,1:i]
      F.hat <- k.eigenvectors*sqrt(len.m-1)
      Lambda.hat <- t(diff.x)%*%F.hat/(len.m-1)
      Comps.hat <- F.hat%*%t(Lambda.hat)
      e <- diff.x-Comps.hat
      # Try with differenced data
      # e <- diff(e)
      for (j in 1:nvar){
        Vj[j] <- 1/len.m*t(e[,j])%*%(e[,j])
      }
      V <- 1/nvar*sum(Vj)
      # Information Criteria for stationary factors
      PC1[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
      PC2[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
      PC3[i] <- V + i*sigma.hat.max*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
      IC1[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
      IC2[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
      IC3[i] <- log(V) + i*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
    }
    
    # Information Criteria
    k.PC1 <- which.min(PC1)
    k.PC2 <- which.min(PC2)
    k.PC3 <- which.min(PC3)
    k.IC1 <- which.min(IC1)
    k.IC2 <- which.min(IC2)
    k.IC3 <- which.min(IC3)
    
    # Double-check with a built-in function
    a <- prcomp(diff.x, rank.=kmax)
    screeplot(a)
    plot(sqrt(kmax.eigenvalues))
    
    #Non-stationary standardized Google trends
    stand.x <- diff.x
    for (i in 1:ncol(diff.x)){
      stand.x[,i] <- cumsum(diff.x[,i])
    }
    stand.x <- as.matrix(stand.x)
    
    # Non-stationary monthly factors
    F.ns.hat <- F.hat.max
    for (i in 1:kmax){
      F.ns.hat[,i] <- cumsum(F.hat.max[,i])
      plot(F.ns.hat[,i], type="l")
    }
    # covariance matrix of the idiosyncratic components:
    idio <- stand.x - F.ns.hat[,1]%*%t(Lambda.hat.max[,1])
    H <- 1/(len.m-1)*t(idio)%*%idio
    
    merged <- matrix(0,(r+1),6)
    for (j in 1:5){
      merged[1:r,j] <- waves[1:r,j]
      merged[(r+1),j] <- NA
    }
    merged[1:(r+1),6] <- CC[1:(r+1),1] # assume you also don't know current value of CC
    merged[(r+1),6] <- NA
    y <- t(cbind(merged,stand.x))
    
    se <- as.data.frame(LFS[1:r,c(2,44,86,128,170)])
    k <- rbind(se,se[nrow(se),])
    
    if (seas_fact==F){
    
      objopt <-  optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                             log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(3000),log(0.02),0,0,log(1000)),
                       KF_slopes_mixed_factor_CC,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                       parP10=1000000000000,nstates=44,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B" )
      par.cc.gt[(r-h+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_CC(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                       nstates=44,lambda=Lambda.hat.max[,1],H=H)
    } else if (seas_fact==T && fixed_sigma_u==T) {
      objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                            log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(3000),log(0.02),log(0.02),0,0,log(1000)), 
                      KF_slopes_mixed_factor_seas_CC,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                      parP10=1000000000000,nstates=55,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
      par.cc.gt[(r-h+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_seas_CC(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                            nstates=55,lambda=Lambda.hat.max[,1],H=H)
    
    } else if (seas_fact==T && fixed_sigma_u==F) {
      objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                            log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(3000),log(0.02),log(1),log(1),0,0,log(1000)), 
                      KF_slopes_mixed_factor_seas_CC,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                      parP10=1000000000000,nstates=55,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
      par.cc.gt[(r-h+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_seas_CC(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                            nstates=55,lambda=Lambda.hat.max[,1],H=H)
    }
    
    P.L.cc.gt[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]))[length(unlist(lapply(obj$Ptt, function(x) x[1,1])))]
    P.R.cc.gt[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[2,2]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
    P.theta.cc.gt[r-h+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                                 2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
    #ts.plot(unlist(lapply(obj$Ptt, function(x) x[1,1]))[31:ncol(y)])
    
    #estimated R
    R_y_est_tt.cc.gt[r-h+1] <- obj$xtt[2,ncol(y)] 
    R_y_est_ttm1.cc.gt[r-h+1] <- obj$xttm1[2,ncol(y)]
    
    #estimated L
    L_y_est_tt.cc.gt[r-h+1] <- obj$xtt[1,ncol(y)] 
    L_y_est_ttm1.cc.gt[r-h+1] <- obj$xttm1[1,ncol(y)]
    
    #estimated seasonal components
    season_y_est_tt.cc.gt <- obj$xtt[3,ncol(y)]+obj$xtt[5,ncol(y)]+obj$xtt[7,ncol(y)]+obj$xtt[9,ncol(y)]+obj$xtt[11,ncol(y)]+obj$xtt[13,ncol(y)]
    season_y_est_ttm1.cc.gt <- obj$xttm1[3,ncol(y)]+obj$xttm1[5,ncol(y)]+obj$xttm1[7,ncol(y)]+obj$xttm1[9,ncol(y)]+obj$xttm1[11,ncol(y)]+obj$xttm1[13,ncol(y)]
    
    #estimated theta
    theta_y_est_tt.cc.gt[r-h+1] <- obj$xtt[1,ncol(y)] + season_y_est_tt.cc.gt 
    theta_y_est_ttm1.cc.gt[r-h+1] <- obj$xttm1[1,ncol(y)] + season_y_est_ttm1.cc.gt
    
    print(paste(P.L.cc.gt[r-h+1],P.R.cc.gt[r-h+1],P.theta.cc.gt[r-h+1]))
  }
  
  now_df_cc_gt <- cbind(P.L.cc.gt,P.R.cc.gt,P.theta.cc.gt,R_y_est_tt.cc.gt,R_y_est_ttm1.cc.gt,L_y_est_tt.cc.gt,L_y_est_ttm1.cc.gt,theta_y_est_tt.cc.gt,theta_y_est_ttm1.cc.gt)
  if (avg == F && fixed == T && load.w == F) {
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt.csv"))
  } else if (avg == T && fixed == F && load.w == F) {
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt_avg.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt_avg.csv"))
  } else if (avg == F && fixed == F && load.w == T) {
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt_load.w.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt_load.w.csv"))
  } else if (avg == T && fixed == F && load.w == T) {
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt_avg_load.w.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt_avg_load.w.csv"))
  } else { # all FALSE
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt.csv"))
  }
  
  if (EN==T){
    saveRDS(sel.var, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/sel.var.gt.cc.csv"))
  }
  
  
  
} else if (monthlyGT==F && CCGT==T){
  len.m <- nrow(x)
  len <- nrow(x.w)
  dataset <- as.data.frame(dataset)
  week <- as.factor(dataset[1:nrow(dataset),1])
  week <- as.Date(week,"%Y-%m-%d")
  dataset_month <- as.data.frame(dataset_month)
  month <- as.factor(dataset_month[1:nrow(dataset_month),1])
  month <- as.Date(month,"%Y-%m-%d")
  h <- len.m-1 - trunc((len.m-1)/3) # forecast in third part of the sample 
  for (j in 1:len){
    if (strftime(month[h+2], "%Y-%m") == strftime(week[j], "%Y-%m")) break
    h.week <- j+1
  }
  h.month <- h+2
  
  P.L.cc.gt <- rep(NA,(len-h.week+1))
  P.R.cc.gt <- rep(NA,(len-h.week+1))
  P.theta.cc.gt <- rep(NA,(len-h.week+1))
  if (seas_fact==F){
    par.cc.gt <- matrix(NA, (len-h.week+1), 13)
  } else if (seas_fact==T && fixed_sigma_u==T) {
    par.cc.gt <- matrix(NA, (len-h.week+1), 14)
  } else if (seas_fact==T && fixed_sigma_u==F) {
    par.cc.gt <- matrix(NA, (len-h.week+1), 15)
  }
  R_y_est_tt.cc.gt <- rep(NA,(len-h.week+1))
  R_y_est_ttm1.cc.gt <- rep(NA,(len-h.week+1))
  L_y_est_tt.cc.gt <- rep(NA,(len-h.week+1))
  L_y_est_ttm1.cc.gt <- rep(NA,(len-h.week+1))
  theta_y_est_tt.cc.gt <- rep(NA,(len-h.week+1))
  theta_y_est_ttm1.cc.gt <- rep(NA,(len-h.week+1))
  sel.var <- lapply(seq_len(len-h.week+1), function(X) rep(NA,ncol(dataset)))
  
  for (r in h.week:len){
    
    len.m <- nrow(x)
    len <- nrow(x.w)
    dataset <- as.data.frame(dataset)
    week <- as.factor(dataset[1:nrow(dataset),1])
    week <- as.Date(week,"%Y-%m-%d")
    dataset_month <- as.data.frame(dataset_month)
    month <- as.factor(dataset_month[1:nrow(dataset_month),1])
    month <- as.Date(month,"%Y-%m-%d")
    
    for (j in 1:len.m){
      if (strftime(month[j], "%Y-%m") == strftime(week[r], "%Y-%m")) break
      h <- j+1
    }
    
    # Aggregate weekly Google trends
    dataset <- as.data.frame(dataset)
    week <- as.factor(dataset[1:r,1])
    week <- as.Date(week,"%Y-%m-%d")
    agg.gt <- matrix(NA,nrow(x[1:h,]),ncol(x))
    agg.gt2 <- matrix(NA,nrow(x[1:h,]),ncol(x))
    agg.gt3 <- matrix(NA,nrow(x[1:h,]),ncol(x))
    agg.gt4 <- matrix(NA,nrow(x[1:h,]),ncol(x))
    agg.gt5 <- matrix(NA,nrow(x[1:h,]),ncol(x))
    
    if (avg==T){
      
      # Aggregate as flow variables
      for (i in 1:ncol(x.w)){
        agg.gt[,i] <- flow_avg(x=x.w[1:r,i], k=1, per=week) # flow variables
        agg.gt2 <- flow_avg(x=x.w[1:r,i], k=2, per=week) # flow variables
        agg.gt3 <- flow_avg(x=x.w[1:r,i], k=3, per=week) # flow variables
        agg.gt4 <- flow_avg(x=x.w[1:r,i], k=4, per=week) # flow variables
        agg.gt5 <- flow_avg(x=x.w[1:r,i], k=5, per=week)
        two.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[2]})$x))
        three.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[3]})$x))
        four.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[4]})$x))
        five.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[5]})$x))
        agg.gt[two.obs,i] <- agg.gt2[two.obs]
        agg.gt[three.obs,i] <- agg.gt3[three.obs]
        agg.gt[four.obs,i] <- agg.gt4[four.obs]
        agg.gt[five.obs,i] <- agg.gt5[five.obs]
      }
    } else {
      
      # Aggregate as flow variables
      for (i in 1:ncol(x.w)){
        agg.gt[,i] <- flow_avg(x=x.w[1:r,i], k=1, per=week) # flow variables
        agg.gt2 <- flow(x=x.w[1:r,i], k=2, per=week) # flow variables
        agg.gt3 <- flow(x=x.w[1:r,i], k=3, per=week) # flow variables
        agg.gt4 <- flow(x=x.w[1:r,i], k=4, per=week) # flow variables
        agg.gt5 <- flow(x=x.w[1:r,i], k=5, per=week)
        two.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[2]})$x))
        three.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[3]})$x))
        four.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[4]})$x))
        five.obs <- which(!is.na(aggregate(x.w[1:r,i],by=ym(week),function(x){rev(x)[5]})$x))
        agg.gt[two.obs,i] <- agg.gt2[two.obs]
        agg.gt[three.obs,i] <- agg.gt3[three.obs]
        agg.gt[four.obs,i] <- agg.gt4[four.obs]
        agg.gt[five.obs,i] <- agg.gt5[five.obs]
      }
      
      # Rescale according to largest value on the entire dataset
      for (i in 1:ncol(agg.gt)){
        if (max(agg.gt[,i]) != 0){
          agg.gt[,i] <- agg.gt[,i]*100/max(agg.gt[,i])
        } else {
          agg.gt[,i] <- agg.gt[,i]
        }
      }
      
    }
    
    if (EN ==T){ 
      # Elastic Net
      diff.x <- diff(agg.gt[(32-13):nrow(agg.gt),],1)
      diff.y <- diff(slope[[h-h.month+1]],1)
      for (i in 1:ncol(diff.x)){
        diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
      }
      diff.y <- (diff.y-mean(diff.y))/sd(diff.y)
      alpha_grid <- seq(0.05, 0.95, 0.01) # grid values for alpha
      AIC <- rep(NA,length(alpha_grid))
      tune_lambda <- rep(NA,length(alpha_grid))
      nselect_var <- rep(NA,length(alpha_grid))
      for (i in 1:length(alpha_grid)) {
        en <- ic.glmnet(diff.x, diff.y, crit = "aic", alpha=alpha_grid[i], lambda=seq(0.1, 0, -0.001)) #alpha=1 lasso, alpha=0 ridge
        # plot(en$glmnet,"lambda")
        # plot(en)
        tune_lambda[i] <- en$lambda
        AIC[i] <- en$ic[2]
        nselect_var[i] <-en$nvar-1
        # coef(en)
        # en$glmnet$lambda
      }
      nselect_var[which.min(AIC)]
      fit <- glmnet(diff.x, diff.y, lambda = tune_lambda[which.min(AIC)], alpha = alpha_grid[which.min(AIC)])
      coef(fit)
      sel <- fit$beta@i[-1]+1
      agg.gt <- agg.gt[,sel]
      sel.var[[r-h.week+1]] <- sel
    }
    
    # Differenciate data
    diff.x <- diff(agg.gt,1)
    diff.x.w <- diff(x.w[1:r,],1)
    
    if (EN==T){
      diff.x.w <- diff(x.w[1:r,sel],1)
    }
    
    if (load.w==T){
      diff.x <- diff.x.w
    }
    
    # Standardize the monthly Google trends
    for (i in 1:ncol(diff.x)){
      diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
    }
    diff.x <- as.matrix(diff.x)
    kmax <- 4 # max number of factors
    cov <- diff.x%*%t(diff.x)
    len.m <- nrow(diff.x)
    nvar <- ncol(diff.x)
    kmax.eigenvalues <- eigen(cov)$values[1:kmax]
    kmax.eigenvectors <- eigen(cov)$vectors[,1:kmax] 
    F.hat.max <- kmax.eigenvectors*sqrt(len.m-1) # facotrs
    t(F.hat.max)%*%F.hat.max
    Lambda.hat.max <- t(diff.x)%*%F.hat.max/(len.m-1) # loadings
    Comps.hat.max <- F.hat.max%*%t(Lambda.hat.max) # common components
    e.max <- diff.x-Comps.hat.max # estimated hydiosincratic components (hc)
    
    # We need the kmax variance of the hc for the information criteria (IC)
    sigma.hat.maxj <- c(1:nvar)
    Vj <- sigma.hat.maxj
    for (j in 1:nvar){
      sigma.hat.maxj[j] <- 1/len.m*t(e.max[,j])%*%(e.max[,j])
    }
    sigma.hat.max <- 1/nvar*sum(sigma.hat.maxj)
    
    # Define vectors to store the values for the IC (Bai and Ng)
    PC1 <- c(1:kmax)
    PC2 <- PC1
    PC3 <- PC1
    IC1 <- PC1
    IC2 <- PC1
    IC3 <- PC1
    
    # Loop to choose the number of factors based on IC
    for (i in 1:kmax){
      k.eigenvalues <- eigen(cov)$values[1:i]
      k.eigenvectors <- eigen(cov)$vectors[,1:i]
      F.hat <- k.eigenvectors*sqrt(len.m-1)
      Lambda.hat <- t(diff.x)%*%F.hat/(len.m-1)
      Comps.hat <- F.hat%*%t(Lambda.hat)
      e <- diff.x-Comps.hat
      # Try with differenced data
      # e <- diff(e)
      for (j in 1:nvar){
        Vj[j] <- 1/len.m*t(e[,j])%*%(e[,j])
      }
      V <- 1/nvar*sum(Vj)
      # Information Criteria for stationary factors
      PC1[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
      PC2[i] <- V + i*sigma.hat.max*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
      PC3[i] <- V + i*sigma.hat.max*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
      IC1[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log((len.m-1)*nvar/((len.m-1)+nvar))
      IC2[i] <- log(V) + i*((len.m-1)+nvar)/((len.m-1)*nvar)*log(min((len.m-1),nvar))
      IC3[i] <- log(V) + i*(log(min((len.m-1),nvar))/min((len.m-1),nvar))
    }
    
    # Information Criteria
    k.PC1 <- which.min(PC1)
    k.PC2 <- which.min(PC2)
    k.PC3 <- which.min(PC3)
    k.IC1 <- which.min(IC1)
    k.IC2 <- which.min(IC2)
    k.IC3 <- which.min(IC3)
    
    #Non-stationary standardized Google trends
    stand.x <- diff.x
    for (i in 1:ncol(diff.x)){
      stand.x[,i] <- cumsum(diff.x[,i])
    }
    stand.x <- as.matrix(stand.x)
    
    # Non-stationary monthly factors
    F.ns.hat <- F.hat.max
    for (i in 1:kmax){
      F.ns.hat[,i] <- cumsum(F.hat.max[,i])
      # plot(F.ns.hat[,i], type="l")
    }
    # covariance matrix of the idiosyncratic components:
    idio <- stand.x - F.ns.hat[,1]%*%t(Lambda.hat.max[,1])
    H <- 1/(len.m-1)*t(idio)%*%idio
    
    if (load.w==T){
      diff.x <- diff(agg.gt,1)
      # Standardize the monthly Google trends
      for (i in 1:ncol(diff.x)){
        diff.x[,i] <- (diff.x[,i]-mean(diff.x[,i]))/sd(diff.x[,i])
      }
      diff.x <- as.matrix(diff.x)
      stand.x <- diff.x
      for (i in 1:ncol(diff.x)){
        stand.x[,i] <- cumsum(diff.x[,i])
      }
      stand.x <- as.matrix(stand.x)
    }
    
    # Monthly google trends
    merged <- matrix(0,(h-1),6)
    for (j in 1:5){
      merged[1:(h-1),j] <- waves[1:(h-1),j]
      merged[(h-1),j] <- NA
    }
    merged[1:(h-1),6] <- CC[1:(h-1),1] # assume you also don't know current value of CC
    merged[(h-1),6] <- NA
    y <- t(cbind(merged,rbind(matrix(NA,(nrow(merged)-nrow(stand.x)),ncol(stand.x)),stand.x)))
    
    # Kalman filter estimation
    se <- as.data.frame(LFS[1:(h-1),c(2,44,86,128,170)])
    k <- rbind(se,se[nrow(se),])
    
    if (seas_fact==F){
      possibleError <- tryCatch(objopt <-  optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                                 log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(3000),log(0.02),0,0,log(1000)),
                       KF_slopes_mixed_factor_CC,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                       parP10=1000000000000,nstates=44,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B" )
                       , error = function(e) e)
      if(inherits(possibleError, "error")) next
      par.cc.gt[(r-h.week+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_CC(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                       nstates=44,lambda=Lambda.hat.max[,1],H=H)
    } else if (seas_fact==T && fixed_sigma_u==T) {
      possibleError <- tryCatch(objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                                                      log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(3000),log(0.02),log(0.02),0,0,log(1000)), 
                                                KF_slopes_mixed_factor_seas_CC,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                                                parP10=1000000000000,nstates=55,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
                                , error = function(e) e)
      if(inherits(possibleError, "error")) next
      par.cc.gt[(r-h.week+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_seas_CC(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                            nstates=55,lambda=Lambda.hat.max[,1],H=H)
    } else if (seas_fact==T && fixed_sigma_u==F) {
      possibleError <- tryCatch(objopt <- optim(par=c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                            log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)),log(3000),log(0.02),log(1),log(1),0,0,log(1000)), 
                      KF_slopes_mixed_factor_seas_CC,y=y,opti=T,k=k,delta=0.21,outofsample=T,
                      parP10=1000000000000,nstates=55,lambda=Lambda.hat.max[,1],H=H,  hessian=F, method="L-BFGS-B")
                      , error = function(e) e)
      if(inherits(possibleError, "error")) next
      par.cc.gt[(r-h.week+1),] <- objopt$par
      #exp(par)/(1-0.21^2) # to get the estimates for the standard deviation of the sampling errors
      obj <- KF_slopes_mixed_factor_seas_CC(par=objopt$par,y=y,opti=F,k=k,delta=0.21,outofsample=T,parP10=1000000000000,
                                            nstates=55,lambda=Lambda.hat.max[,1],H=H)
    }
    
    P.L.cc.gt[r-h.week+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]))[length(unlist(lapply(obj$Ptt, function(x) x[1,1])))]
    P.R.cc.gt[r-h.week+1] <- unlist(lapply(obj$Ptt, function(x) x[2,2]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
    P.theta.cc.gt[r-h.week+1] <- unlist(lapply(obj$Ptt, function(x) x[1,1]+x[3,3]+x[5,5]+x[7,7]+x[9,9]+x[11,11]+x[13,13]+
                                                 2*x[1,3]+2*x[3,5]+2*x[5,7]+2*x[7,9]+2*x[9,11]+2*x[11,13]))[length(unlist(lapply(obj$Ptt, function(x) x[2,2])))]
    #ts.plot(unlist(lapply(obj$Ptt, function(x) x[1,1]))[31:ncol(y)])
    
    #estimated R
    R_y_est_tt.cc.gt[r-h.week+1] <- obj$xtt[2,ncol(y)] 
    R_y_est_ttm1.cc.gt[r-h.week+1] <- obj$xttm1[2,ncol(y)]
    
    #estimated L
    L_y_est_tt.cc.gt[r-h.week+1] <- obj$xtt[1,ncol(y)] 
    L_y_est_ttm1.cc.gt[r-h.week+1] <- obj$xttm1[1,ncol(y)]
    
    #estimated seasonal components
    season_y_est_tt.cc.gt <- obj$xtt[3,ncol(y)]+obj$xtt[5,ncol(y)]+obj$xtt[7,ncol(y)]+obj$xtt[9,ncol(y)]+obj$xtt[11,ncol(y)]+obj$xtt[13,ncol(y)]
    season_y_est_ttm1.cc.gt <- obj$xttm1[3,ncol(y)]+obj$xttm1[5,ncol(y)]+obj$xttm1[7,ncol(y)]+obj$xttm1[9,ncol(y)]+obj$xttm1[11,ncol(y)]+obj$xttm1[13,ncol(y)]
    
    #estimated theta
    theta_y_est_tt.cc.gt[r-h.week+1] <- obj$xtt[1,ncol(y)] + season_y_est_tt.cc.gt 
    theta_y_est_ttm1.cc.gt[r-h.week+1] <- obj$xttm1[1,ncol(y)] + season_y_est_ttm1.cc.gt
  
    print(paste(P.L.cc.gt[r-h.week+1],P.R.cc.gt[r-h.week+1],P.theta.cc.gt[r-h.week+1]))
  }
  
  now_df_cc_gt <- cbind(P.L.cc.gt,P.R.cc.gt,P.theta.cc.gt,R_y_est_tt.cc.gt,R_y_est_ttm1.cc.gt,L_y_est_tt.cc.gt,L_y_est_ttm1.cc.gt,theta_y_est_tt.cc.gt,theta_y_est_ttm1.cc.gt)
  if (avg == F && fixed == T && load.w == F) {
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt.csv"))
  } else if (avg == T && fixed == F && load.w == F) {
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt_avg.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt_avg.csv"))
  } else if (avg == F && fixed == F && load.w == T) {
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt_load.w.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt_load.w.csv"))
  } else if (avg == T && fixed == F && load.w == T) {
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt_avg_load.w.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt_avg_load.w.csv"))
  } else { # all FALSE
    saveRDS(now_df_cc_gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt.csv"))
    saveRDS(par.cc.gt, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt.csv"))
  }
  
  if (EN==T){
    saveRDS(sel.var, file.path("Paper/Empirical analysis R/Nowcast results/2004_2017/sel.var.gt.cc.csv"))
  }
}



#### Comparison of nowcast performance #####

if (fixed == T) {
  now_df_univ <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_univ_fixed.csv")
  par.univ <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.univ_fixed.csv")
} else {
  now_df_univ <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_univ.csv")
  par.univ <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.univ.csv")
}

if (fixed == T) {
  now_df_biv <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_biv_fixed.csv")
  par.biv <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.biv_fixed.csv")
} else {
  now_df_biv <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_biv.csv")
  par.biv <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.biv.csv")
}

if (CCGT == FALSE){
  if (avg == F && fixed == T && load.w == F) {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt.csv")
  } else if (avg == T && fixed == F && load.w == F) {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt_avg.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt_avg.csv")
  } else if (avg == F && fixed == F && load.w == T) {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt_load.w.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt_load.w.csv")
  } else if (avg == T && fixed == F && load.w == T) {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt_avg_load.w.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt_avg_load.w.csv")
  } else {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_gt.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.gt.csv")
  }
} else {
  if (avg == F && fixed == T && load.w == F) {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt.csv")
  } else if (avg == T && fixed == F && load.w == F) {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt_avg.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt_avg.csv")
  } else if (avg == F && fixed == F && load.w == T) {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt_load.w.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt_load.w.csv")
  } else if (avg == T && fixed == F && load.w == T) {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt_avg_load.w.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt_avg_load.w.csv")
  } else {
    now_df_gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/now_df_cc_gt.csv")
    par.gt <- readRDS("Paper/Empirical analysis R/Nowcast results/2004_2017/par.cc.gt.csv")
  }
}

P.L <- now_df_univ[,1]
P.R <- now_df_univ[,2]
P.theta <- now_df_univ[,3]

P.L.biv <- now_df_biv[,1]
P.R.biv <- now_df_biv[,2]
P.theta.biv <- now_df_biv[,3]

P.L.gt <- now_df_gt[,1]
P.R.gt <- now_df_gt[,2]
P.theta.gt <- now_df_gt[,3]

mean(P.L.biv)/mean(P.L)
mean(P.R.biv)/mean(P.R)
mean(P.theta.biv)/mean(P.theta)

mean(P.L.gt,na.rm=TRUE)/mean(P.L)
mean(P.R.gt,na.rm=TRUE)/mean(P.R)
mean(P.theta.gt,na.rm=TRUE)/mean(P.theta)

len <- nrow(x.w)
len.m <- nrow(diff.x)

P.L.gt.week <- matrix(NA, nrow=5, ncol=length(month[(h.month-1):len.m]))
for (i in 1:5){
  P.L.gt.week[i,] <- aggregate(P.L.gt,by=ym(week[h.week:len]),function(x){x[i]})$x
}
apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE)

P.R.gt.week <- matrix(NA, nrow=5, ncol=length(month[(h.month-1):len.m]))
for (i in 1:5){
  P.R.gt.week[i,] <- aggregate(P.R.gt,by=ym(week[h.week:len]),function(x){x[i]})$x
}
apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE)

P.theta.gt.week <- matrix(NA, nrow=5, ncol=length(month[(h.month-1):len.m]))
for (i in 1:5){
  P.theta.gt.week[i,] <- aggregate(P.theta.gt,by=ym(week[h.week:len]),function(x){x[i]})$x
}
apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)

now_df_gt_last <- matrix(NA, length(month[(h.month-1):len.m]), 3)
now_df_gt_last[,1] <- aggregate(now_df_gt[,4],by=ym(week[h.week:len]),function(x){rev(x)[1]})$x
now_df_gt_last[,2] <- aggregate(now_df_gt[,6],by=ym(week[h.week:len]),function(x){rev(x)[1]})$x
now_df_gt_last[,3] <- aggregate(now_df_gt[,8],by=ym(week[h.week:len]),function(x){rev(x)[1]})$x



#### Save nowcast results ####

if (CCGT == FALSE){
  if (avg==F && fixed==T && load.w==F){
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast_fixed.csv", row.names=FALSE)
  } else if (avg==T && fixed==F && load.w==F){
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast_avg.csv", row.names=FALSE)
  } else if (avg==F && fixed==F && load.w==T){
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast_load.w.csv", row.names=FALSE)
  } else if (avg==T && fixed==F && load.w==T){
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast_avg_load.w.csv", row.names=FALSE)  
  } else { #all false
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast.csv", row.names=FALSE)
  }
} else {
  if (avg==F && fixed==T && load.w==F){
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast_fixed_cc_gt.csv", row.names=FALSE)
  } else if (avg==T && fixed==F && load.w==F){
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast_avg_cc_gt.csv", row.names=FALSE)
  } else if (avg==F && fixed==F && load.w==T){
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast_load.w_cc_gt.csv", row.names=FALSE)
  } else if (avg==T && fixed==F && load.w==T){
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast_avg_load.w_cc_gt.csv", row.names=FALSE)  
  } else { #all false
    write.csv(data.frame(month[h.month:(len.m+1)],cbind(now_df_univ[,c(4,6,8)],now_df_biv[,c(4,6,8)],now_df_gt_last,waves[(h.month+nrow(waves)-nrow(diff.x)-1):(len.m+nrow(waves)-nrow(diff.x)),],diff(as.matrix(waves),1)[(h.month+nrow(waves)-nrow(diff.x)-2):(len.m+nrow(waves)-nrow(diff.x)-1),])), 
              file = "Paper/Empirical analysis R/Nowcast results/2004_2017/nowcast_cc_gt.csv", row.names=FALSE)
  }
}


if (CCGT == FALSE){
  if (avg==F && fixed==T && load.w==F){
    write.xlsx(c(mean(P.L.biv)/mean(P.L),mean(P.R.biv)/mean(P.R),mean(P.theta.biv)/mean(P.theta)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_CC_fixed.xlsx")
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_gt_fixed.xlsx")
  } else if (avg==T && fixed==F && load.w==F){
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_gt_avg.xlsx")
  } else if (avg==F && fixed==F && load.w==T){
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_gt_load.w.xlsx")
  } else if (avg==T && fixed==F && load.w==T){
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_gt_avg_load.w.xlsx")
  } else { # all FALSE
    write.xlsx(c(mean(P.L.biv)/mean(P.L),mean(P.R.biv)/mean(P.R),mean(P.theta.biv)/mean(P.theta)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_CC.xlsx")
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_gt.xlsx")
  }
} else {
  if (avg==F && fixed==T && load.w==F){
    write.xlsx(c(mean(P.L.biv)/mean(P.L),mean(P.R.biv)/mean(P.R),mean(P.theta.biv)/mean(P.theta)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_CC_fixed.xlsx")
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_cc_gt_fixed.xlsx")
  } else if (avg==T && fixed==F && load.w==F){
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_cc_gt_avg.xlsx")
  } else if (avg==F && fixed==F && load.w==T){
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_cc_gt_load.w.xlsx")
  } else if (avg==T && fixed==F && load.w==T){
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_cc_gt_avg_load.w.xlsx")
  } else { # all FALSE
    write.xlsx(c(mean(P.L.biv)/mean(P.L),mean(P.R.biv)/mean(P.R),mean(P.theta.biv)/mean(P.theta)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_CC.xlsx")
    write.xlsx(c(mean(P.L.gt,na.rm=TRUE)/mean(P.L),apply(P.L.gt.week/mean(P.L),1,mean,na.rm=TRUE),
                 mean(P.R.gt,na.rm=TRUE)/mean(P.R),apply(P.R.gt.week/mean(P.R),1,mean,na.rm=TRUE),
                 mean(P.theta.gt,na.rm=TRUE)/mean(P.theta),apply(P.theta.gt.week/mean(P.theta),1,mean,na.rm=TRUE)), 
               "Paper/Empirical analysis R/Nowcast results/2004_2017/now_cc_gt.xlsx")
  }
}