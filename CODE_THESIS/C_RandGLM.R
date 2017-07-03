## Author: Robert Lange
## DATE: 06/2017
## BGSE Master Project- RandNLA for GLMS with Big Datasets

####################################################
####################################################
##Applying Fast LS Methods to Generalized LinModel##
####################################################
####################################################
#Idea: Use approximations in iterative weight least squares procedure!
#Presentation - work on weighted leverage scores!
#Canonical Link Function for Binomial Data
inv.link <- function(eta){ exp(eta)/(1+exp(eta))}

data.gen.logit <- function(N,M,beta,type){
  y <- matrix(0,N,1)
  X <- matrix(0,N,M)
  eta <- matrix(0,N,1)
  mu <- matrix(0,N,1)

  X <- data_lev_gen(N, M, type)

  eta <- X%*%beta
  eta
  mu <- inv.link(eta)
  mu[is.na(mu)] <- 1
  y <- runif(N) < mu
  data <- data.frame(y,X)
  data
}

logit.mean <- function(eta){exp(eta)/(1+exp(eta))}
logit.var <- function(mu){mu*(1-mu)}

##########################################################

iwls.logit <- function(y, X, max.iter, delta){

  beta <- rep(0,ncol(X))
  beta.prev <- beta
  beta.trace <- matrix(0, max.iter+1, length(beta))

  eta <- X%*%beta
  mu <- (y+0.5)/2
  beta.trace[1,] <- beta

  for(iter in 1:max.iter){
    cat(".")

    z <- eta + (y-mu)/logit.var(mu)
    w <- logit.var(mu)
    W <- diag(as.numeric(w))

    beta <- solve(t(X)%*%W%*%X, t(X)%*%W%*%z)
    eta <- X %*% beta
    mu <- logit.mean(eta)

    beta.trace[iter+1, ] <- beta

    if(sum((beta-beta.prev)**2) < delta){return(beta.trace[1:iter+1,])}
    beta.prev <- beta
    }
  cat("\n")
  beta.trace[1:iter+1,]
}

iwls.logit.cov <- function(y,X,beta){
  mu <- logit.mean(X%*%beta)
  w <- logit.var(mu)
  W <- diag(as.numeric(w))
  C <- solve(t(X)%*%W%*%X)
  C
}

##########################################################

slow.iwls.logit <- function(y,X,max.iter,eps, delta,type="weighted"){
  beta <- rep(0,ncol(X))
  beta.prev <- beta
  beta.trace <- matrix(0, max.iter+1, length(beta))

  eta <- X%*%beta
  mu <- (y+0.5)/2
  beta.trace[1,] <- beta

  for(iter in 1:max.iter){
    cat(".")

    z <- eta + (y-mu)/logit.var(mu)
    if((abs(sum(z))==Inf) | is.nan(sum(z))==1){
        beta.trace[1:iter+1,]
        break
    }
    w <- logit.var(mu)
    W <- diag(as.numeric(w))

    if (type=="unweighted"){
        beta <- slow.rSampling.GLM(X,W,z,type="unweighted",eps)
    }else if (type=="weighted"){
        beta <- slow.rSampling.GLM(X,W,z,type="weighted",eps)
    }else if (type=="influence"){
        beta <- slow.rSampling.GLM(X,W,z,type="influence",eps)
    }

    eta <- X %*% beta
    mu <- logit.mean(eta)

    beta.trace[iter+1, ] <- beta
    if(sum((beta - beta.prev)**2)<delta){return(beta.trace[1:iter+1,])}
    beta.prev <- beta
  }
  cat("\n")
  beta.trace[1:iter+1,]
}

slow.rSampling.GLM <- function(X,W,z,type,eps){
  #INPUT: Lin. sytem of eq./Euclidean Norm min problem (min||Ax-b||_2), sampling type ("SVD" or "QR"), approx. error - eps)
  #OUTPUT: Approximate solution to system based on random sampling: x_opt
  #NOTES: dimensionality (d) and desired approx. error determine number of sampled rows
  n <- dim(X)[1]
  d <- dim(X)[2]
  r <- round((d*log(d))/eps)

  sampled <- sampling.GLM(A=X,B=z,W=W,c=r,type)

  SX <- matrix(unlist(sampled[1]), ncol=d)
  Sz <- matrix(unlist(sampled[2]), ncol=1)
  SW <- matrix(unlist(sampled[3]), ncol=r)

  x_opt <- solve(t(SX)%*%SW%*%SX, t(SX)%*%SW%*%Sz)
  x_opt
}


sampling.GLM <- function(A,B,W,c,type="unweighted"){
  #INPUT: Matrices A, B, positive integer c, probability type on which sampling is based
  #OUTPUT: c sampled columns and rows
  #NOTES: QR/SVD used for LS
  m= dim(A)[1]
  p= dim(B)[2]
  n <- dim(A)[2] #or dim(B)[1]

  #sampling for LS/solving of Linear systems of Equations
  length <- seq(1,m)
  C <- matrix(NA,c,n)
  R <- c(NA)
  H <- matrix(0,c,c)
  if (type == "unweighted"){
      P <- sampling.prob.GLM(A,W, type="unweighted")
  } else if (type == "weighted"){
      P <- sampling.prob.GLM(A,W, type="weighted")
  } else if (type == "influence"){
      P <- sampling.prob.GLM(A,W=W,z=B, "influence")
  }

  for (i in 1:c){
    #sample in the easiest way - sample uniformly - p_i/k = 1/n
    i_c = sample(length,1,prob=P,replace=T)
    C[i,] = A[i_c,]/sqrt(c*P[i_c])
    R[i] = B[i_c]/sqrt(c*P[i_c])
    H[i,i] = W[i_c,i_c]
  }

  Out <- list("SA" = C, "SB" = R, "SW" = H)
  Out
}


sampling.prob.GLM <- function(A,W,z=0,type="unweighted"){
  #INPUT: Matrix A,probability type on which sampling is based (SVD/QR) or sLev
  #OUTPUT: c sampled columns and rows
  #NOTES: Both versions yield exactly the same probabilities; Use "thin" SVD - things that are zeroed-out by singular values are not included
    d <- dim(A)[2]
    w <- diag(W)
    #Probabilities based on SVD

    if (type == "unweighted"){
        U <- fast.svd(A)$u
        p = eucl_matrix_row_norm(U)^2/d
    } else if (type == "weighted"){
        U <- fast.svd(A)$u
        p = eucl_matrix_row_norm(w*U)^2/norm(W%*%U, type="F")
    } else if (type == "influence"){
        Xz <- cbind(A,z)
        #Probabilities based on SVD
        U <- fast.svd(Xz)$u
        p = eucl_matrix_row_norm(U)^2/norm(U, type="F")
    }
  return(p)
}
