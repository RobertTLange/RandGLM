## Author: Robert Lange
## DATE: 06/2017
## BGSE Master Project- RandNLA for GLMS with Big Datasets

#Description: This file provides slow LS approximation algorithms
#such as the ones introduced in the first few chapters of Mahoney (2016)
#Furthermore, I code the shrinked leveraging estimator (SLEV) and the
#Unweighted Leveraging Estimator (LEVUNW) first introduced by Ma, et al. (2015).

####################################################
####################################################
##Slow Sampling/Projection algorithms for LS##
####################################################
####################################################

fast.approx.samp.dist <- function(X, eps){
    # Only first approximation!!!
    # Second approx makes running time very slow!
    n <- dim(X)[1]
    d <- dim(X)[2]

    r1 <- round((d*log(n))/eps^2 * log((d*log(n))/eps^2))
    Pi1 <- matrix(rnorm(r1*n), r1, n)

    lev_approx = eucl_matrix_row_norm(X %*% ginv(Pi1 %*% X))^2
    samp_dist <- lev_approx/sum(lev_approx)
    samp_dist
}

sampling.prob.LS <- function(A,b=0,p.type,q=0){
  #INPUT: Matrix A,probability type on which sampling is based (SVD/QR) or sLev
  #OUTPUT: c sampled columns and rows
  #NOTES: Both versions yield exactly the same probabilities; Use "thin" SVD - things that are zeroed-out by singular values are not included
  d <- dim(A)[2]
  if (p.type=="SVD"){
    #Probabilities based on SVD
    U <- fast.svd(A)$u
    p = (1/d)*eucl_matrix_row_norm(U)^2
  }else if (p.type=="QR"){
    #Probabilities based on QR
    QR <- qr(A)
    Q <- qr.Q(QR)
    p = (1/d)*eucl_matrix_row_norm(Q)^2
  }else if (p.type=="sLev"){
    #Probabilities based on QR
    n <- dim(A)[1]
    U <- fast.svd(A)$u
    p = (1/d)*eucl_matrix_row_norm(U)^2
    p <- q*p + (1-q)*(1/n)
  }else if (p.type=="influence"){
    Ab <- cbind(A,b)
    U <- fast.svd(Ab)$u
    p = (1/d)*eucl_matrix_row_norm(U)^2
  }else if (p.type=="approx.lev"){
    p = fast.approx.samp.dist(A, eps=0.2)
  }else if (p.type=="approx.infl"){
    Ab <- cbind(A,b)
    p = fast.approx.samp.dist(Ab, eps=0.2)
  }
  return(p)
}

#sampling.prob.LS(A,b=b,p.type="SVD",q=0) - sampling.prob.LS(A,b=b,p.type,q=0)

##########################################################
sampling.LS <- function(A,B,c,p.type,q=0){
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
  P <- sampling.prob.LS(A,b=B,p.type,q)
  for (i in 1:c){
    #sample in the easiest way - sample uniformly - p_i/k = 1/n
    i_c = sample(length,1,prob=P,replace=T)
    C[i,] = A[i_c,]/sqrt(c*P[i_c])
    R[i] = B[i_c]/sqrt(c*P[i_c])
  }

  Out <- list("SA" = C, "SB" = R)
  Out
}

##########################################################

slow.rSampling.LS <- function(A,b,p.type,eps){
  #INPUT: Lin. sytem of eq./Euclidean Norm min problem (min||Ax-b||_2), sampling type ("SVD" or "QR"), approx. error - eps)
  #OUTPUT: Approximate solution to system based on random sampling: x_opt
  #NOTES: dimensionality (d) and desired approx. error determine number of sampled rows
  n <- dim(A)[1]
  d <- dim(A)[2]
  r <- round((d*log(d))/eps)
  sampled <- sampling.LS(A=A,B=b,c=r,p.type)
  SA <- matrix(unlist(sampled[1]), ncol=d)
  Sb <- matrix(unlist(sampled[2]), ncol=1)
  x_opt <- solve(t(SA)%*%SA,t(SA)%*%Sb)
  x_opt
}

#beta.est <- slow.rSampling.LS(A,b,p.type,eps)
#beta.base <- solve(t(A)%*%A)%*%t(A)%*%b

qual.of.approx <- function(A,b,beta.est, beta.base,eps){
    #Computes quality of approx. statistic from Drineas(2011)
    #0.8 of the probability mass should lie below 1
    EN.error <- eucl_vec_norm(beta.est - beta.base)
    U <- fast.svd(A)$u
    gamma <- eucl_vec_norm(U %*% t(U) %*% b)/eucl_vec_norm(b)
    denom <- sqrt(eps)*(kappa(A)*sqrt(gamma^(-2) - 1))* eucl_vec_norm(beta.base)
    quality <- EN.error/denom
    quality
}

#qual.of.approx(A,b,beta.est, beta.base,eps)

#beta.Sampling <- slow.rSampling.LS(A,b,"SVD",eps=0.01)

##########################################################

eps.JLT <- function(A,b,r,n,proj.type){
  #INPUT: Lin. sytem of eq./Euclidean Norm min problem (min||Ax-b||_2), projection type ("FM" or "Ach"), approx. error - eps)
  #OUTPUT: Matrix and vector projected to r-dim subspace
  #NOTES: dimensionality (d) and desired approx. error determine number of sampled rows
  #NOTES: FM - Frankl and Meahara (Gaussian), Ach - Achlioptas (binary coin toss)
  #k <- (9*log(n))/(eps^2 - eps^3)
  if (proj.type=="FM"){
    P.fill <- rnorm((n*r))
    P <- (1/sqrt(r))*matrix(P.fill,r,n)
  }else if (proj.type=="Ach"){
    P.fill <- sample(c(3/sqrt(r),0,-3/sqrt(r)),size=((n*r)),replace=T,prob=c(1/6,2/3,1/6))
    P <- matrix(P.fill,r,n)
  }
  A.proj <- P%*%A
  if (is.null(b) == F){
    b.proj <- S%*%H%*%D%*%b
    Out <- list("SA" = A.proj, "Sb" = b.proj)
    Out
  } else{
    A.proj
  }
}

#projected <- projection.LS(A,b,r,n,"FM")

##########################################################

slow.rProjection.LS <- function(A,b,eps,proj.type){
  #INPUT: Lin. sytem of eq./Euclidean Norm min problem (min||Ax-b||_2), projection type ("FM" or "Ach"), approx. error - eps)
  #OUTPUT: Approximate solution to system based on random projection on row space: x_opt
  #NOTES: dimensionality (d) and desired approx. error determine number of sampled rows
  #NOTES: FM - Frankl and Meahara (Gaussian), Ach - Achlioptas (binary coin toss)
  d <- dim(A)[2]
  n <- dim(A)[1]
  r <- round((d*log(d))/eps)
  projected <- eps.JLT(A,b,r,n,proj.type)
  SA <- matrix(unlist(projected[1]), ncol=d)
  Sb <- matrix(unlist(projected[2]), ncol=1)
  x_opt <- solve(t(SA)%*%SA,t(SA)%*%Sb)
  x_opt
}

#beta.Projection <- slow.rProjection.LS(A,b,eps=0.01,"FM")

##########################################################

sLev.est <- function(A,b,q,eps){
    n <- dim(A)[1]
    d <- dim(A)[2]
    r <- round((d*log(d))/eps)
    sampled <- sampling.LS(A=A,B=b,c=r,p.type="sLev",q)
    SA <- matrix(unlist(sampled[1]), ncol=d)
    Sb <- matrix(unlist(sampled[2]), ncol=1)
    x_opt <- solve(t(SA)%*%SA,t(SA)%*%Sb)
    x_opt
}

#beta.sLev <- sLev.est(A,b,q=0.4,eps=0.05)

##########################################################

Lev.unw <- function(A,B,c,p.type){
    m= dim(A)[1]
    p= dim(B)[2]
    n <- dim(A)[2] #or dim(B)[1]

    #sampling for LS/solving of Linear systems of Equations
    length <- seq(1,m)
    C <- matrix(NA,c,n)
    R <- c(NA)
    P <- sampling.prob.LS(A,p.type,q)
    for (i in 1:c){
      #sample in the easiest way - sample uniformly - p_i/k = 1/n
      i_c = sample(length,1,prob=P,replace=T)
      C[i,] = A[i_c,]
      R[i] = B[i_c]
    }
    Out <- list("SA" = C, "SB" = R)
    Out
}

Lev.unw.est <- function(A,b,p.type,eps){
  #INPUT: Lin. sytem of eq./Euclidean Norm min problem (min||Ax-b||_2), sampling type ("SVD" or "QR"), approx. error - eps)
  #OUTPUT: Approximate solution to system based on random sampling: x_opt
  #NOTES: dimensionality (d) and desired approx. error determine number of sampled rows
  n <- dim(A)[1]
  d <- dim(A)[2]
  r <- round((d*log(d))/eps)
  sampled <- Lev.unw(A=A,B=b,c=r,p.type)
  SA <- matrix(unlist(sampled[1]), ncol=d)
  Sb <- matrix(unlist(sampled[2]), ncol=1)
  x_opt <- solve(t(SA)%*%SA,t(SA)%*%Sb)
  x_opt
}

#beta.Lev.unw <- Lev.unw.est(A,b,p.type="SVD",eps=0.05)

##########################################################

# print(system.time(solve(t(A)%*%A,t(A)%*%b)))
# print(system.time(slow.rSampling.LS(A,b,"SVD",eps=0.01)))
# print(system.time(slow.rSampling.LS(A,b,"QR",eps=0.01)))
# print(system.time(slow.rProjection.LS(A,b,eps=0.01,"FM")))
# print(system.time(slow.rProjection.LS(A,b,eps=0.01,"Ach")))
