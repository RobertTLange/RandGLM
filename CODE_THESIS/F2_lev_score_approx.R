## Author: Robert Lange
## DATE: 06/2017
## BGSE Master Project- RandNLA for GLMS with Big Datasets

set.seed(101)

################################################################################

fast.approx.lev.scores <- function(X,eps){
  n <- dim(X)[1]
  d <- dim(X)[2]

  r1 <- round((d*log(n))/eps^2 * log((d*log(n))/eps^2))

  X.proj <- eps.FJLT(X,r1,n)
  R.proj <- qr.R(qr(X.proj))

  r2 <- round(log(n)/eps^2)
  XR <- X%*%solve(R.proj)
  XR.proj <- eps.JLT(XR,r2,dim(XR)[2])
  p = data.frame(eucl_matrix_row_norm(XR.proj)^2)
  p
}


eps.FJLT <- function(X,r,n){
  P <- Matrix(0,r,n,sparse=T)
  ind <- sample(seq(1,n,1), size=r, replace=T)

  for (i in 1:r){
    P[i,ind[i]] <- sqrt(n/r)
  }

  H <- hadamard.matrix(n)

  D <- Matrix(0,n,n)
  diag(D) <- sample(c(-1,1), size=n, replace=T)

  PH <- P%*%H
  epsFJLT <- PH%*%D
  X.proj <- epsFJLT%*%X
  X.proj
}

eps.JLT <- function(X,r,d){
  P.fill <- rnorm((d*r))
  P <- (1/sqrt(r))*matrix(P.fill,d,r)

  X.proj <- X%*%P
  X.proj
}

###################################################################################

approx.samp.dist <- function(X, eps, type){
    n <- dim(X)[1]
    d <- dim(X)[2]

    r1 <- round((d*log(n))/eps^2 * log((d*log(n))/eps^2))
    Pi1 <- matrix(rnorm(r1*n), r1, n)

    if (type == "hat"){
        lev_approx = eucl_matrix_row_norm(X %*% ginv(Pi1 %*% X))^2
        samp_dist <- lev_approx/sum(lev_approx)
    }
    else if (type == "tilde"){
        r2 <- round(log(n)/eps^2)
        Pi2 <- matrix(rnorm(r1*r2), r1, r2)
        lev_approx = eucl_matrix_row_norm(X %*% ginv(Pi1 %*% X) %*% Pi2)^2
        samp_dist <- lev_approx/sum(lev_approx)
    }

    samp_dist
}


stat.QoA <- function(X, eps, type){
    samp.dist.approx <- approx.samp.dist(X, eps, type)
    samp.dist.exact <- LEV(X)/dim(X)[2]
    stat <- abs(samp.dist.exact - samp.dist.approx)/(eps*samp.dist.exact)
    stat
}

plot.stat.QoA <- function(X, eps, type){
    statQoA <- stat.QoA(X, eps, type)
    plot <- ggplot(statQoA, aes(x=leverage)) +
      geom_histogram(col="blue", aes(y=..density..), fill="deepskyblue1") +
      geom_density(aes(y=..density..)) +
      ylab(TeX('Count/Density')) +
      #geom_vline(xintercept = 1) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(0, 0.5)
    if (type=="hat"){
        plot = plot + xlab(TeX('$\\frac{|p_i - \\hat{p}_i|}{\\epsilon p_i}$'))
    }
    else if (type=="tilde"){
        plot = plot + xlab(TeX('$\\frac{|p_i - \\tilde{p}_i|}{\\epsilon p_i}$'))
    }
    plot
}

X=X_GA

X_GA <- data_lev_gen(n, d, type="nearly uniform lev")
X_T1 <- data_lev_gen(n, d, type="moderately nonuniform lev")
X_T3 <- data_lev_gen(n, d, type="very nonuniform lev")

p1 <- plot.stat.QoA(X_GA, eps, type="hat") + ggtitle("a) Nearly Uniform \n  Multivariate Normal")
p2 <- plot.stat.QoA(X_GA, eps, type="tilde")
p3 <- plot.stat.QoA(X_T1, eps, type="hat") + ggtitle("b) Moderately Nonunif. \n Multivariate t (df=1)")
p4 <- plot.stat.QoA(X_T1, eps, type="tilde")
p5 <- plot.stat.QoA(X_T3, eps, type="hat") + ggtitle("c) Very Nonuniform \n Multivariate t (df=3)")
p6 <- plot.stat.QoA(X_T3, eps, type="tilde")

if(save.plots) pdf("F2_Lev_Score_Approx.pdf")
multiplot(p1, p2, p3, p4, p5, p6, cols=3)
if(save.plots) dev.off()
