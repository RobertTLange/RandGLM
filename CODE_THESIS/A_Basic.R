## Author: Robert Lange
## DATE: 06/2017
## BGSE Master Project- RandNLA for GLMS with Big Datasets

# CONTENT: USEFUL BASIC FUNCTIONS FOR FURTHER ANALYSIS#

################################################################################
## Data Simulation Functions##
################################################################################
LEV <- function(X) {
    U <- fast.svd(X)$u
    leverage = eucl_matrix_row_norm(U)^2
    data.frame(leverage)
}

data_lev_gen <- function(n, d, type = "nearly uniform lev") {
    # Function generates data according to Ma et. al (2015)
    Sigma <- diag(rep(2,d))
    for (i in 1:d){
        for (j in 1:d){
            Sigma[i,j] <- 2*0.5^(abs(i-j))
        }
    }

    if (type == "nearly uniform lev") {
        X <- rmvnorm(n, mean=rep(0, nrow(Sigma)), sigma=Sigma)
    } else if (type == "moderately nonuniform lev") {
        X <- rmvt(n, sigma=Sigma, df=3)
    } else if (type == "very nonuniform lev") {
        X <- rmvt(n, sigma=Sigma, df=1)
    }
    X
}

data_sim_gauss_indep <- function(mu=0, sigma=1, N, M, cov.structure=0) {
  #Function simulates data from the univariate normal distribution
  #Assume that variables are iid
    X <- matrix(NA, N, M)
        if (cov.structure == 0){
            X[,] <- rnorm(N*M,0,1)
        } else{
            Sigma = diag(M)
            Sigma[lower.tri(Sigma)] = rep(cov.structure, M*(M-1)/2)
            Sigma[upper.tri(Sigma)] = Sigma[lower.tri(Sigma)]
            X = mvrnorm(N*M, mu=rep(0, M), Sigma)
        }
    return(X)
}

#data_sim_gauss_indep(N=2, M=3)

gen_cov_matrix_dense <- function(M){
  #Many dependence relations in the features
  #Creates dense PSD and symmetric matrix which is used as a covariance matrix in the following simulations
  A <- data_sim_gauss_indep(0,1,M,M)[,-1]
  #Construct symmetric matrix and ensure PSD by adding nI
  #Also: A <- A%*%t(A) - The first is significantly faster: O(n^2) compared to O(n^3)
  A <- 0.5*(A + t(A))
  A <- A + M*diag(1,M)
  return(A)
}

#m=1000; n=1000; p=500
#A <- data_sim_gauss_indep(N=m, M=n)
#B <- data_sim_gauss_indep(N=n, M=p)

################################################################################
##Matrix and Vector Norms##
################################################################################

eucl_vec_norm <- function(x){
  sqrt(sum(x^2))
}

eucl_matrix_row_norm <- function(A){
  rows <- dim(A)[1]
  row_norms <- c()
  for (i in 1:rows){
    row_norms[i] <- sqrt(t(A[i,])%*%A[i,])
  }
  row_norms
}

################################################################################
##Algo 1: Select Algorithm##
################################################################################

select <- function(n,alpha,beta){
  #INPUT: n - length of vector a; alpha,beta - parameters
  #OUTPUT: a.star - randomly sampled element from vector, i.star - index of element
  #NOTES: Sample vector from beta to ensure nonnegativity
  a <- rbeta(n,alpha,beta)
  D=0
  a.star=0
  i.star=0
  for (i in 1:length(a)){
    D = D + a[i]
    if (a[i]/D > runif(1,0,1)){
      #Sample from uniform to ensure that a[i] is selected with correct prob.
      i.star=i
      a.star=a[i]
    }
  }
  list(i.star, a.star)
}

#system.time(select(100000,2,3))

################################################################################
##Algo 2: Inner product Matrix Multiplication (3-loop mat mult algorithm)##
################################################################################

##View an element of AB as an inner product between row of A and an column of B
inner.MatMult <- function(A,B){
  #INPUT: Matrix A,B (#column A = #row B)
  #OUTPUT: Matrix product AB (mxn)(nxp)
  #NOTES: High-school algorithm for MatMatMult

  m= dim(A)[1]
  p= dim(B)[2]

  AB <- matrix(NA, m, p)

  #for (i in 1:m){
  # #Directly compute dot-products - faster -> vectorized
  #  for (j in 1:p){
  #    AB[i,j] = t(A[i,])%*%B[,j]
  #  }
  #}

  for (i in 1:m){
    for (j in 1:p){
      AB[i,j]=0
      for (k in 1:n){
        AB[i,j] = AB[i,j] + A[i,k]*B[k,j]
      }
    }
  }
  return(AB)
}

#Frobenius Norm - Spectral Norm: type="2"
#norm(inner.MatMult(A,B)-A%*%B, type = "F")


################################################################################
##MatrixMultiplication as sum of outer products algorithm##
################################################################################

outer.MatMult <- function(A,B){
  #INPUT: Matrix A,B (#column A = #row B)
  #OUTPUT: Matrix product AB (mxn)(nxp)
  #NOTES: Computes Matrix Matri Product by summing over Rank One Matrices

  m= dim(A)[1]
  p= dim(B)[2]
  n <- dim(A)[2] #or dim(B)[1]
  AB_outer = matrix(0,m,p)

  for(i in 1:n){
    AB_outer = AB_outer + A[,i]%*%t(B[i,])
    #assign(paste0("AB", i, sep = ""), AB_outer)
  }
  return(AB_outer)
}

#norm(outer.MatMult(A,B)-A%*%B, type = "F")

################################################################################
##Solve LS using Cholesky Decomposition##
################################################################################

my.chol <- function(A){
  #INPUT: Symmetric Positive Definite Matrix
  #OUTPUT: Lower Triangular Matrix derived from Cholesky Decomposition
  n <- dim(A)[1]
  L <- matrix(c(rep(0,n*n)),n,n)

  #Construct matrix by looping over elements in matrix
  for (i in 1:n){
    for (j in 1:i){
      s <- t(L[i,])%*%L[j,]
      if (i == j){
        L[i,j] = sqrt(A[i,i] - s)
      } else{
        L[i,j] = (1/L[j,j]) * (A[i,j] - s)
      }
    }
  }
  return(L)
}

my.forward.solve <- function(L, b){
  #INPUT: Lower Triangular Matrix L and vector of same dimension b
  #OUTPUT: Vector x that solves linear system of equations Lx=b
  n <- length(b)
  x <- rep(0,n)

  for (i in 1:n){
    s = t(L[i,])%*%x
    x[i] = (1/L[i,i])*(b[i] - s)
  }

  return(x)
}

my.back.solve <- function(U, b){
  #INPUT: Upper Triangular Matrix U and vector of same dimension b
  #OUTPUT: Vector x that solves linear system of equations Ux=b
  n <- length(b)
  x <- rep(0,n)

  for (i in n:1){
    s = t(U[i,])%*%x
    x[i] = (1/U[i,i])*(b[i] - s)
  }

  return(x)
}

exact.LS <- function(A, b){
  #INPUT: Symmetric Positive Definite Matrix and vector of same dimension b
  #OUTPUT: Vector x that solves linear system of equations Ax=b
  Gramm <- t(A)%*%A
  L <- my.chol(Gramm)
  y <- my.forward.solve(L,t(A)%*%b)
  x <- my.back.solve(t(L),y)
  return(x)
}

################################################################################
##Hadamard matrix generation functions
################################################################################

legendre.symbol <- function(a, p){
  if(isprime(p)){
    K <- as.bigz(a)^(as.bigz((p-1)/2))
    modulus(K) <- p
    K[K==p-1] <- -1
    k <- integer(length(K))
    k[K==-1] <- -1
    k[K==1] <- 1
    A <- as.bigz(a)
    modulus(A) <- p
    k[A==0] <- 0
    k
  } else stop("p not prime")
}

jacobsthal <- function(p){
  if(isprime(p)){
    Q=outer(1:p,1:p,function(i,j)legendre.symbol(i-j, p))
    return(Q)
  }else
  stop("p not prime")
}

paley <- function(n)
{
  if(isprime(n-1) && (n-1)%%4==3) {
                                        # Jacobsthal Matrix Q:
    Q=jacobsthal(n-1)
                                        # H =( 1 1^T)
                                        #    ( 1 Q-I)
    H=rbind(c(0,rep(1,n-1)),
            cbind(rep(-1,n-1), Q))+diag(1,n)
    return(H)
  }else stop("n-1 is not prime or n%%4 != 0")
}

paley2 <- function(n){

  if( n%%4==0) {
    N=n/2
    if(isprime(N-1) && (N-1)%%4==1) {
      #Jacobsthal Matrix Q:
      Q=jacobsthal(N-1);
      # H =( 1 1^T)
      #    ( 1 Q-E)

      S0=matrix(c(1,-1,-1,-1),2,2,byrow=TRUE)
      S1=matrix(c(1,1,1,-1),2,2,byrow=TRUE)
      S=rbind(c(0,rep(1,N-1)),
              cbind(rep(1,N-1), Q))

      H=(S==0)%x%S0 + S%x%S1

      return(H);
    }else stop("n/2-1 is not prime or n/2%%4 != 2")
  }else stop("n%%4 != 0")
}



walsh <- function(H)
{
  N=dim(H)[1];
  if(N!=dim(H)[2])
    stop("H not square")
  H2=sylvester(1)%x%H

  return(H2)
}

sylvester <- function(k)
{
  n=2^k;
  if(n==1)
    return(matrix(1,1,1));
  if(n==2)
    return(matrix(c(1, 1, 1, -1),2,2,byrow=TRUE));
  N=2
  H=sylvester(1)
  while(N<n){
    H=walsh(H)
    N=2*N
  }
  return(H)
}

hadamard.matrix <- function(n, verbose=FALSE, getnearest=TRUE)
{
  if(n==1)
    return(sylvester(1))
  if(n==2)
    return(sylvester(2))

  f <- factorize(n)
  if(all(f==2)){
    pow <- length(f)
    if(verbose)cat("Hadamard matrix of sylvester type\n")
    return(sylvester(pow))
  } else {
    if(n%%4==0){
      if(isprime(n-1)){
        if(verbose)cat("Hadamard matrix of paley I type\n")
        H=paley(n)
        return(H)
      }else{
        if(isprime(n/2-1) && (n/2-1)%%4==1) {
          if(verbose)cat("Hadamard matrix of paley II type\n")
          H=paley2(n)
          return(H)
        }else{
                                        # try n/2
          if(verbose)cat("try recursive call (walsh step)\n")
          H2=hadamard.matrix(n/2,getnearest=FALSE)
          if(!is.null(H2)){
            H=walsh(H2)
            return(H)
          }else{
            if(verbose) cat(paste("use stored matrix for size=",n,"\n"))
            hadamard.table(n)
          }
        }
      }
    }else{
      if(getnearest){
        warning(paste(n," != 4*k, get the smallest Hadamard matrix with at least dimension ",n))
        hadamard.matrix((n%/%4 +1)*4)
      }else{
        warning(paste(n," != 4*k"))
        return(NULL)
      }
    }
  }
}

#dim(hadamard.matrix(8192))

is.hadamard <- function(H){
  n <- dim(H)[1]
  all((H%*%t(H)-n*diag(1,n))==0)
}
