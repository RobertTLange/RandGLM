## Author: Robert Lange
## DATE: 06/2017
## BGSE Master Project- RandNLA for GLMS with Big Datasets

####################################################
##Approximation Error of RandLS estimator         ##
####################################################

no_cores <- detectCores() - 1

cl<-makeCluster(no_cores, outfile="")
registerDoParallel(cl)

######################################################################
#LS Estimator Simulations for three different data-gen Processes
######################################################################

LS.GA <- foreach(i = 1:iter, .combine=rbind, .packages=c('corpcor', 'mvtnorm', 'MASS')) %dopar%{
    X <- data_lev_gen(n, d, type="nearly uniform lev")
    beta <- seq(1,d,1)
    y <- X%*%beta + rnorm(dim(X)[1])

    beta.Chol <- exact.LS(X,y)
    beta.Rsamp <- slow.rSampling.LS(X,y,"approx.lev",eps=eps)
    beta.infl <- slow.rSampling.LS(X,y,"approx.infl",eps=eps)

    exact.GA <- qual.of.approx(X, y, beta.Rsamp, beta, eps=eps)
    chol.GA <- qual.of.approx(X, y, beta.Rsamp, beta.Chol, eps=eps)
    chol.influence.GA <- qual.of.approx(X, y, beta.infl, beta.Chol, eps=eps)

    out <- cbind(exact.GA, chol.GA, chol.influence.GA)
}

print("RandLS for GA done")

LS.T1 <- foreach(i = 1:iter, .combine=rbind, .packages=c('corpcor', 'mvtnorm', 'MASS')) %dopar%{
    X <- data_lev_gen(n, d, type="moderately nonuniform lev")
    beta <- seq(1,d,1)
    y <- X%*%beta + rnorm(dim(X)[1])

    beta.Chol <- exact.LS(X,y)
    beta.Rsamp <- slow.rSampling.LS(X,y,"SVD",eps=eps)
    beta.infl <- slow.rSampling.LS(X,y,"influence",eps=eps)

    exact.T1 <- qual.of.approx(X, y, beta.Rsamp, beta, eps=eps)
    chol.T1 <- qual.of.approx(X, y, beta.Rsamp, beta.Chol, eps=eps)
    chol.influence.T1 <- qual.of.approx(X, y, beta.infl, beta.Chol, eps=eps)

    out <- cbind(exact.T1, chol.T1, chol.influence.T1)
}

print("RandLS for T1 done")

LS.T3 <- foreach(i = 1:iter, .combine=rbind, .packages=c('corpcor', 'mvtnorm', 'MASS')) %dopar%{
    X <- data_lev_gen(n, d, type="very nonuniform lev")
    beta <- seq(1,d,1)
    y <- X%*%beta + rnorm(dim(X)[1])

    beta.Chol <- exact.LS(X,y)
    beta.Rsamp <- slow.rSampling.LS(X,y,"SVD",eps=eps)
    beta.infl <- slow.rSampling.LS(X,y,"influence",eps=eps)

    exact.T3 <- qual.of.approx(X, y, beta.Rsamp, beta, eps=eps)
    chol.T3 <- qual.of.approx(X, y, beta.Rsamp, beta.Chol, eps=eps)
    chol.influence.T3 <- qual.of.approx(X, y, beta.infl, beta.Chol, eps=0.01)

    out <- cbind(exact.T3, chol.T3, chol.influence.T3)
}

print("RandLS for T3 done")

stopCluster(cl)

######################################################################
data.full <- data.frame(cbind(LS.GA, LS.T1, LS.T3))
if(save.sim){
    write.table(data.full, "sim_ls_results.txt", sep="\t")
}

plot_QoA_LS <- function(data, data.to.plot, type){
    plot <- ggplot(data, aes_string(x=data.to.plot)) +
      geom_histogram(col="blue", aes(y=..density..), fill="deepskyblue1") +
      geom_density(aes(y=..density..)) +
      geom_vline(xintercept = 1) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(0, 1.3)
     if (type == "true"){
         plot = plot +
                xlab(TeX('$\\frac{||\\tilde{\\beta}_{RS} - \\beta_0||_2}{\\sqrt{\\epsilon} \\sqrt{\\kappa(X)^{-2} -1} ||\\beta_0||_2}$'))
     } else if (type == "est"){
         plot = plot +
                xlab(TeX('$\\frac{||\\tilde{\\beta}_{RS} - \\beta_{LS}^{CHOL}||_2}{\\sqrt{\\epsilon} \\sqrt{\\kappa(X)^{-2} -1} ||\\beta_{LS}^{CHOL}||_2}$'))
     } else if (type == "influence"){
         plot = plot +
                xlab(TeX('$\\frac{||\\tilde{\\beta}_{Influence} - \\beta_{LS}^{CHOL}||_2}{\\sqrt{\\epsilon} \\sqrt{\\kappa(X)^{-2} -1} ||\\beta_{LS}^{CHOL}||_2}$'))
     }
    plot
}

# data.full <- data.frame(read.table("sim_ls_results.txt" , header=T))

RSAMP.exact.GA <- plot_QoA_LS(data=data.full, data.to.plot="exact.GA", type="true") +
                ggtitle('a) Nearly Uniform \n Multivariate Normal \n QoA (Lev) - True Vec') +
                ylab(TeX('Count/Density'))
RSAMP.exact.T1 <- plot_QoA_LS(data=data.full, data.to.plot="exact.T1", type="true") +
                ggtitle('b) Moderately Nonuniform \n Multivariate t (df=1) \n QoA (Lev) - True Vec') +
                theme(axis.title.y=element_blank())
RSAMP.exact.T3 <- plot_QoA_LS(data=data.full, data.to.plot="exact.T3", type="true") +
                ggtitle('c) Very Nonuniform \n Multivariate t (df=3) \n QoA (Lev) - True Vec') +
                theme(axis.title.y=element_blank())

RSAMP.chol.GA <- plot_QoA_LS(data=data.full, data.to.plot="chol.GA", type="est") +
                ggtitle('QoA (Lev) - LS') +
                ylab(TeX('Count/Density'))
RSAMP.chol.T1 <- plot_QoA_LS(data=data.full, data.to.plot="chol.T1", type="est") +
                ggtitle('QoA (Lev) - LS') +
                theme(axis.title.y=element_blank())
RSAMP.chol.T3 <- plot_QoA_LS(data=data.full, data.to.plot="chol.T3", type="est") +
                ggtitle('QoA (Lev) - LS') +
                theme(axis.title.y=element_blank())

Infl.chol.GA <- plot_QoA_LS(data=data.full, data.to.plot="chol.influence.GA", type="influence") +
                ggtitle('QoA (Influence) - LS') +
                ylab(TeX('Count/Density'))
Infl.chol.T1 <- plot_QoA_LS(data=data.full, data.to.plot="chol.influence.T1", type="influence") +
                ggtitle('QoA (Influence) - LS') +
                theme(axis.title.y=element_blank())
Infl.chol.T3 <- plot_QoA_LS(data=data.full, data.to.plot="chol.influence.T3", type="influence") +
                ggtitle('QoA (Influence) - LS') +
                theme(axis.title.y=element_blank())

if(save.plots) pdf("F3_LS_comparison.pdf")
multiplot(RSAMP.exact.GA, RSAMP.chol.GA, Infl.chol.GA, RSAMP.exact.T1,  RSAMP.chol.T1, Infl.chol.T1, RSAMP.exact.T3, RSAMP.chol.T3,  Infl.chol.T3, cols=3)
if(save.plots) dev.off()
