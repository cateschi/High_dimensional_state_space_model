# This R-script is party based on "An Introduction to State Space Models" by Marc Wildi.

# The function "LFS_GT" performs the estimation of the LFS model, with the multivariate auxiliary series of Google Trends. It requires the following arguments:
# par: initial values for the model's parameters ((9+r)x1 vector)g; r is the number of Google Trends' factors.
# y: (5+n)xT matrix of the unemployed labour force and the n Google trends (the first 5 series are the unemployed labour force) (T=185).
# opti: if TRUE, optimizes the function.
# k: Tx5 matrix of the standard errors of the GREG estimates.
# outofsample: if TRUE, computes the loglikelihood based on the out-of-sample forecast errors.
# parP10: large number for the diffuse initialization.
# nfactors: number of Google Trends' factors (r).
# nstates: number of state variables in the model.
# lambda: nxr vector of estimated, by PCA, factors' loadings.
# H: the nxn estimated covariance matrix, by PCA, of the idiosyncratic components of the Google Trends.
# ns.id: vector of length equal to the number of nonstationary idiosyncratic components. It should contain the ordered numbers corresponing to the nonstationary idiosyncratic components.

# The lines that are not commented here, are commented on /LFS_model.R.  

# Packages required to run the scripts:
library(magic)
library(ucminf)


is.scalar <- function(x) is.atomic(x) && length(x) == 1L


LFS_GT <- function(par,y,opti,k,outofsample,parP10,nfactors,nstates,lambda,H,ns.id){
    
  len <- length(y[1,])
  sigma_Ry <- par[1]
  sigma_omegay <- par[2]
  sigma_lambda <- par[3]
  sd_nu <- diag(exp(c(par[4], par[5], par[6], par[7], par[8])), 5,5)
  sigma_Rx <- log(1)      # variance of the Google Trends' factors' innovations (it is fixed).
  x10 <- rep(0,nstates)
  Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
  Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
  delta <- par[9]
  P10 <- diag(c(rep(parP10[1],17),c(1,rep((1-delta^2),4),1,rep((1-delta^2),3),1,rep((1-delta^2),3)),rep(parP10[1],nstates-30)),nstates,nstates)     
  Pttm1[[1]] <- P10
  xtt <- matrix(0,nstates,(len))
  xttm1 <- matrix(0,nstates,(len+1))
  xttm1[,1] <- x10
  R <- diag(1,nstates,nstates)
  D <- adiag(0, exp(sigma_Ry), exp(sigma_omegay)*diag(11), exp(sigma_lambda)*diag(4), sd_nu, diag(0,8,8), diag(rep(exp(sigma_Rx),nfactors)), diag(sqrt(diag(H)[ns.id])))
  for (j in 1:nfactors){
    R[(30+j),2] <- tanh(par[9+j])      # correlation between the LFS solpe's and the Google Trends' factors' innovations.
    R[2,(30+j)] <- tanh(par[9+j])
  }
  Q <- D%*%R%*%D
  
  # Build T (the transition matrix)::
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
  TyE <- cbind(TyE, rbind(c(0,0,0,0),diag(delta,nrow=4,ncol=4),matrix(0,8,4)))
  TyE <- cbind(TyE, rbind(matrix(0,5,4),diag(4),matrix(0,4,4)))
  Ty <- adiag(Tymu, Tyomega, Tylambda, TyE)
  Tx <- diag(1,ncol=(length(ns.id)+nfactors), nrow=(length(ns.id)+nfactors))      # transition matrix of the Google Trends' factors and the nonstationary idiosyncratic components.
  Tmatrix <- adiag(Ty, Tx)
  
  # initialization of log-likelihood
  logl <- 0
  
  # Start of KF recursions
  for (i in 1:len){
    
    # Bulild Z:
    Zy <- c(1,0)
    Zy <- rep(Zy,6)
    Zy <- c(Zy,1)
    Zy <- rbind(Zy,Zy,Zy,Zy,Zy)
    Zy <- cbind(Zy,rbind(c(0,0,0,0),diag(4)))
    if (is.na(k[i,])) {
      k[i,which(is.na(k[i,]))] <- 0
    }
    Zy <- cbind(Zy, diag(as.numeric(k[i,]), nrow=5, ncol=5), matrix(0, nrow=5, ncol=8))
    one.ns.id <- rep(0, nrow(H))
    one.ns.id[ns.id] <- 1       # one if element of vector corresponds to nonstationary idiosyncratic component
    Zx <- cbind(as.matrix(lambda),diag(one.ns.id))
    Zx <- Zx[,which(!apply(Zx,2,FUN = function(x){all(x == 0)}))]
    Z <- cbind(adiag(Zy,as.matrix(Zx)))
    ncol(Z)
    nrow(Z)
    
    W <- diag(1,length(y[,i]))     
    if (length(which(is.na(y[,i]))) > 0 && length(which(is.na(y[,i]))) < length(y[,i])){      
      W <- matrix(W[-which(is.na(y[,i])),], nrow=(length(y[,i])-length(which(is.na(y[,i])))), ncol=ncol(W))
      Z <- W%*%Z
      y[which(is.na(y[,i])),i] <- 0
    }
    
    if (length(which(is.na(y[,i]))) > 0 && length(which(is.na(y[,i]))) == length(y[,i])){  
      xtt[,i] <- xttm1[,i]
      Ptt[[i]] <- Pttm1[[i]]
      Pttm1[[i+1]] <- Tmatrix%*%Pttm1[[i]]%*%t(Tmatrix) + Q
      xttm1[,i+1] <- Tmatrix%*%xttm1[,i]
    } else {
      epshatoutofsample <- W%*%y[,i] - Z%*%xttm1[,i]
      diag.H <- diag(H)
      diag.H[ns.id] <- 0       # only variances of stationary idiosyncratic components.
      Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + W%*%adiag(diag(0,ncol(waves),ncol(waves)),diag(diag.H,nrow(y)-ncol(waves),nrow(y)-ncol(waves)))%*%t(W)
      #Fmatrix[1,1] <- ifelse(!is.na(epshatoutofsample[1,]), Fmatrix[1,1], parP10[1])
      if (is.scalar(Fmatrix) == TRUE){
        Fmatrix.inv = 1/Fmatrix
      } else {
        svdFmatrix <- svd(Fmatrix)
        Fmatrix.inv <- svdFmatrix$v%*%diag(1/svdFmatrix$d)%*%t(svdFmatrix$u)
      }
      Kg <- Tmatrix%*%Pttm1[[i]]%*%t(Z)%*%Fmatrix.inv 
      xtt[,i] <- xttm1[,i]+Pttm1[[i]]%*%t(Z)%*%Fmatrix.inv%*%epshatoutofsample 
      epshatinsample <- W%*%y[,i]-Z%*%xtt[,i] 
      Ptt[[i]] <- Pttm1[[i]]-Pttm1[[i]]%*%t(Z)%*%Fmatrix.inv%*%Z%*%Pttm1[[i]] 
      Pttm1[[i+1]] <- Tmatrix%*%Pttm1[[i]]%*%t(Tmatrix-Kg%*%Z)+Q 
      xttm1[,i+1] <- Tmatrix%*%xttm1[,i] + Kg%*%epshatoutofsample
    }
    
    # The optimization criterion:
    if (outofsample) {
      if (i <= (30-13) ){
        logl <- logl - nrow(y)/2*log(2*pi)
      } else if (i > (30-13) ){ 
        logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%Fmatrix.inv%*%epshatoutofsample
        if ((NaN %in% logl)==T){
          logl<- -P10[1]
        }
      }
    } else {
      if (i <= (30-13) ){
        logl <- logl - nrow(y)/2*log(2*pi)
      } else if (i > (30-13) ){ 
        logl <- logl - nrow(y)/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%Fmatrix.inv%*%epshatinsample
        if ((NaN %in% logl)==T){
          logl<- -P10[1]
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


init.val.GT <- c(log(2000),log(0.02),log(900),log(1.07),log(0.99*(1-0.21^2)),
                 log(1.01*(1-0.21^2)),log(1.13*(1-0.21^2)),log(1.06*(1-0.21^2)), 0.21, rep(0,nfactors))

objopt.GT <- ucminf(par=init.val.GT, LFS_GT,y,opti=T,k,outofsample=T,parP10=1000000000000,nfactors,nstates,lambda,H=H,
                        ns.id,hessian=2,control=list(grad="central", gradstep = c(1e-2, 1e-3), trace=T))

par.GT <- objopt.GT$par
                
obj <- LFS_GT(par=objopt.GT$par,y,opti=F,k,outofsample=T,parP10=1000000000000,nfactors,nstates,lambda,H=H,ns.id)
