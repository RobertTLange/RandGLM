## Author: Robert Lange
## DATE: 06/2017
## BGSE Master Project- RandNLA for GLMS with Big Datasets

####################################################
##Approximation Error of RandGLM estimator        ##
####################################################

#####SET SEED for reproducibility####
set.seed(1234, kind = NULL, normal.kind = NULL)
####################################################

beta = c(0.1,0.2,0.3,0.4,0.5)

no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)

######################################################################
#GLM Estimator Simulations - Generate One Randomized Estimator Per Simulated Data Set!
######################################################################
GLM.GA <- foreach(i = 1:iter, .combine=rbind, .packages=c('corpcor', 'mvtnorm', 'MASS')) %dopar%{
    data <- data.gen.logit(n,d,beta,type="nearly uniform lev")
    y <- as.matrix(data[,1])
    X <- as.matrix(data[,1:d+1])

    #beta.iwls.trace <- iwls.logit(y,X,max.iter,delta)
    beta.iwls <- glm(y ~ X -1, binomial(link="logit"))$coefficients
    #iter.iwls.GA <- dim(beta.iwls.trace)[1]

    beta.RandIWLS.trace <- slow.iwls.logit(y,X,max.iter,eps,delta,type="unweighted")
    beta.sampling.unw <- beta.RandIWLS.trace[nrow(beta.RandIWLS.trace),]
    iter.sampling.unweighted.GA <- dim(beta.RandIWLS.trace)[1]

    beta.RandIWLS.trace <- slow.iwls.logit(y,X,max.iter,eps,delta,type="weighted")
    beta.sampling.w <- beta.RandIWLS.trace[nrow(beta.RandIWLS.trace),]
    iter.sampling.weighted.GA <- dim(beta.RandIWLS.trace)[1]

    beta.RandIWLS.trace <- slow.iwls.logit(y,X,max.iter,eps,delta,type="influence")
    beta.sampling.infl <- beta.RandIWLS.trace[nrow(beta.RandIWLS.trace),]
    iter.sampling.infl.GA <- dim(beta.RandIWLS.trace)[1]

    EN.iwls <- eucl_vec_norm(beta.iwls)
    EN.unweighted.GA <- eucl_vec_norm(beta.sampling.unw - beta.iwls)/EN.iwls
    EN.weighted.GA <- eucl_vec_norm(beta.sampling.w - beta.iwls)/EN.iwls
    EN.influence.GA <- eucl_vec_norm(beta.sampling.infl - beta.iwls)/EN.iwls

    out <- cbind(EN.unweighted.GA, EN.weighted.GA, EN.influence.GA)
}

print("RandGLM for GA done")

GLM.T1 <- foreach(i = 1:iter, .combine=rbind, .packages=c('corpcor', 'mvtnorm', 'MASS'))  %dopar%{
    data <- data.gen.logit(n,d,beta,type="moderately nonuniform lev")
    y <- as.matrix(data[,1])
    X <- as.matrix(data[,1:d+1])

    #beta.iwls.trace <- iwls.logit(y,X,max.iter,delta)
    beta.iwls <- glm(y ~ X -1, binomial(link="logit"))$coefficients
    #iter.iwls.T1 <- dim(beta.iwls.trace)[1]

    beta.RandIWLS.trace <- slow.iwls.logit(y,X,max.iter,eps,delta,type="unweighted")
    beta.sampling.unw <- beta.RandIWLS.trace[nrow(beta.RandIWLS.trace),]
    iter.sampling.unweighted.T1 <- dim(beta.RandIWLS.trace)[1]

    beta.RandIWLS.trace <- slow.iwls.logit(y,X,max.iter,eps,delta,type="weighted")
    beta.sampling.w <- beta.RandIWLS.trace[nrow(beta.RandIWLS.trace),]
    iter.sampling.weighted.T1 <- dim(beta.RandIWLS.trace)[1]

    beta.RandIWLS.trace <- slow.iwls.logit(y,X,max.iter,eps,delta,type="influence")
    beta.sampling.infl <- beta.RandIWLS.trace[nrow(beta.RandIWLS.trace),]
    iter.sampling.infl.T1 <- dim(beta.RandIWLS.trace)[1]

    EN.iwls <- eucl_vec_norm(beta.iwls)
    EN.unweighted.T1 <- eucl_vec_norm(beta.sampling.unw - beta.iwls)/EN.iwls
    EN.weighted.T1 <- eucl_vec_norm(beta.sampling.w - beta.iwls)/EN.iwls
    EN.influence.T1 <- eucl_vec_norm(beta.sampling.infl - beta.iwls)/EN.iwls

    out <- cbind(EN.unweighted.T1, EN.weighted.T1, EN.influence.T1)
}

print("RandGLM for T1 done")

GLM.T3 <- foreach(i = 1:iter, .combine=rbind, .packages=c('corpcor', 'mvtnorm'))  %dopar%{
    data <- data.gen.logit(n,d,beta,type="moderately nonuniform lev")
    y <- as.matrix(data[,1])
    X <- as.matrix(data[,1:d+1])

    #beta.iwls.trace <- iwls.logit(y,X,max.iter,delta)
    beta.iwls <- glm(y ~ X -1, binomial(link="logit"))$coefficients
    #iter.iwls.T3 <-  dim(beta.iwls.trace)[1]

    beta.RandIWLS.trace <- slow.iwls.logit(y,X,max.iter,eps,delta,type="unweighted")
    beta.sampling.unw <- beta.RandIWLS.trace[nrow(beta.RandIWLS.trace),]
    iter.sampling.unweighted.T3 <- dim(beta.RandIWLS.trace)[1]

    beta.RandIWLS.trace <- slow.iwls.logit(y,X,max.iter,eps,delta,type="weighted")
    beta.sampling.w <- beta.RandIWLS.trace[nrow(beta.RandIWLS.trace),]
    iter.sampling.weighted.T3 <- dim(beta.RandIWLS.trace)[1]

    beta.RandIWLS.trace <- slow.iwls.logit(y,X,max.iter,eps,delta,type="influence")
    beta.sampling.infl <- beta.RandIWLS.trace[nrow(beta.RandIWLS.trace),]
    iter.sampling.infl.T3 <- dim(beta.RandIWLS.trace)[1]

    EN.iwls <- eucl_vec_norm(beta.iwls)
    EN.unweighted.T3 <- eucl_vec_norm(beta.sampling.unw - beta.iwls)/EN.iwls
    EN.weighted.T3 <- eucl_vec_norm(beta.sampling.w - beta.iwls)/EN.iwls
    EN.influence.T3 <- eucl_vec_norm(beta.sampling.infl - beta.iwls)/EN.iwls

    out <- cbind(EN.unweighted.T3, EN.weighted.T3, EN.influence.T3)
}

print("RandGLM for T3 done")

stopCluster(cl)

######################################################################
# Outfile Simulation Results
data.full <- data.frame(cbind(GLM.GA, GLM.T1, GLM.T3))

if(save.sim){
    write.table(data.full, "sim_glm_results.txt", sep="\t")
}
######################################################################
# Plot and save Simulation Results
plot_QoA_GLM <- function(data, data.to.plot, type){
    plot <- ggplot(data, aes_string(x=data.to.plot)) +
      geom_histogram(col="blue", aes(y=..density..), fill="deepskyblue1") +
      geom_density(aes(y=..density..)) +
      ylab("Count/Density") +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(0, 1.5)
    if (type == "weighted"){
        plot = plot +
                xlab(TeX('$\\frac{||\\tilde{\\beta}^{weightedLev}_{RandGLM} - \\beta_{IWLS}||_2}{||\\beta_{IWLS}||_2}$'))
    } else if (type == "unweighted"){
        plot = plot +
                xlab(TeX('$\\frac{||\\tilde{\\beta}^{unweightedLev}_{RandGLM} - \\beta_{IWLS}||_2}{||\\beta_{IWLS}||_2}$'))
    } else if (type == "influence"){
        plot = plot +
               xlab(TeX('$\\frac{||\\tilde{\\beta}^{influence}_{RandGLM} - \\beta_{IWLS}||_2}{||\\beta_{IWLS}||_2}$'))
    }
    plot
}

GLM.GA_unw <- plot_QoA_GLM(data=data.full, data.to.plot="EN.unweighted.GA", type="unweighted") +
                ggtitle('a) Nearly Uniform \n Multivariate Normal \n Unw. Lev Scores')
GLM.T1_unw <- plot_QoA_GLM(data=data.full, data.to.plot="EN.unweighted.T1", type="unweighted") +
                ggtitle('b) Moderately Nonuniform \n Multivariate t (df=1) \n Unw. Lev Scores') +
                theme(axis.title.y=element_blank())
GLM.T3_unw <- plot_QoA_GLM(data=data.full, data.to.plot="EN.unweighted.T3", type="unweighted") +
                ggtitle('c) Very Nonuniform \n Multivariate t (df=3) \n Unw. Lev Scores') +
                theme(axis.title.y=element_blank())


GLM.GA_w <- plot_QoA_GLM(data=data.full, data.to.plot="EN.weighted.GA", type="weighted") +
                ggtitle('Weighted')
GLM.T1_w <- plot_QoA_GLM(data=data.full, data.to.plot="EN.weighted.T1", type="weighted") +
                ggtitle('Weighted') +
                theme(axis.title.y=element_blank())
GLM.T3_w <- plot_QoA_GLM(data=data.full, data.to.plot="EN.weighted.T3", type="weighted") +
                ggtitle('Weighted') +
                theme(axis.title.y=element_blank())

GLM.GA_infl <- plot_QoA_GLM(data=data.full, data.to.plot="EN.influence.GA", type="influence") +
                ggtitle('Influence')
GLM.T1_infl <- plot_QoA_GLM(data=data.full, data.to.plot="EN.influence.T1", type="influence") +
                ggtitle('Influence') +
                theme(axis.title.y=element_blank())
GLM.T3_infl <- plot_QoA_GLM(data=data.full, data.to.plot="EN.influence.T3", type="influence") +
                ggtitle('Influence') +
                theme(axis.title.y=element_blank())

if(save.plots) pdf("F4_GLM_comparison.pdf")
multiplot(GLM.GA_unw, GLM.GA_w, GLM.GA_infl, GLM.T1_unw, GLM.T1_w, GLM.T1_infl, GLM.T3_unw, GLM.T3_w, GLM.T3_infl , cols=3)
if(save.plots) dev.off()
