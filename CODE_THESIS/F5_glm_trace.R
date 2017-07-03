## Author: Robert Lange
## DATE: 06/2017
## BGSE Master Project- RandNLA for GLMS with Big Datasets

####################################################
set.seed(1234, kind = NULL, normal.kind = NULL)

#2^13 = 8192
beta = c(0.1,0.15,-0.2,0.25,-0.3)

legend.names = c("iteration", "IWLS", "Unweighted RandGLM", "Weighted RandGLM", "Influence RandGLM")

plot.trace <- function(data){
    newdf <- melt(data,'iteration')
    plot  <- ggplot(newdf,aes(x=iteration,y=value,group=variable,color=variable)) +
                geom_line() +
                ylab(TeX('$\\beta^{RandGLM}_{(k)}$')) +
                xlab("Iteration") +
                theme(plot.title = element_text(hjust = 0.5, lineheight=.8)) +
                ylim(-0.5, 0)
    plot
}
############################################################################
data <- data.gen.logit(n,d,beta,type="nearly uniform lev")
y <- as.matrix(data[,1])
X <- as.matrix(data[,1:d+1])

beta.iwls.trace.GA <- iwls.logit(y,X,max.iter,delta)[,5]
conv.iwls.GA <- length(beta.iwls.trace.GA)

beta.RandIWLS.trace.GA.unw <- slow.iwls.logit(y,X,max.iter,eps,delta,type="unweighted")[,5]
beta.RandIWLS.trace.GA.w <- slow.iwls.logit(y,X,max.iter,eps,delta,type="weighted")[,5]
beta.RandIWLS.trace.GA.infl <- slow.iwls.logit(y,X,max.iter,eps,delta,type="influence")[,5]

max.iwls.iter.GA = max(length(beta.iwls.trace.GA),
                       length(beta.RandIWLS.trace.GA.unw),
                       length(beta.RandIWLS.trace.GA.w))

if (length(beta.iwls.trace.GA) < max.iwls.iter.GA){
    beta.iwls.trace.GA <- c(beta.iwls.trace.GA,
                            rep(beta.iwls.trace.GA[length(beta.iwls.trace.GA)],
                            max.iwls.iter.GA-length(beta.iwls.trace.GA)))
}
if (length(beta.RandIWLS.trace.GA.unw) < max.iwls.iter.GA){
    beta.RandIWLS.trace.GA.unw <- c(beta.RandIWLS.trace.GA.unw,
                                    rep(beta.RandIWLS.trace.GA.unw[length(beta.RandIWLS.trace.GA.unw)],
                                    max.iwls.iter.GA-length(beta.RandIWLS.trace.GA.unw)))
}
if (length(beta.RandIWLS.trace.GA.w) < max.iwls.iter.GA){
    beta.RandIWLS.trace.GA.w <- c(beta.RandIWLS.trace.GA.w,
                                  rep(beta.RandIWLS.trace.GA.w[length(beta.RandIWLS.trace.GA.w)],
                                  max.iwls.iter.GA-length(beta.RandIWLS.trace.GA.w)))
}
if (length(beta.RandIWLS.trace.GA.infl) < max.iwls.iter.GA){
    beta.RandIWLS.trace.GA.infl <- c(beta.RandIWLS.trace.GA.infl,
                                     rep(beta.RandIWLS.trace.GA.infl[length(beta.RandIWLS.trace.GA.infl)],
                                     max.iwls.iter.GA-length(beta.RandIWLS.trace.GA.infl)))
}


data.to.plot.GA <- as.data.frame(cbind(seq(1, length(beta.iwls.trace.GA)), beta.iwls.trace.GA, beta.RandIWLS.trace.GA.unw, beta.RandIWLS.trace.GA.w, beta.RandIWLS.trace.GA.infl))
colnames(data.to.plot.GA) = legend.names
plot.GA <- plot.trace(data.to.plot.GA) +
            geom_vline(xintercept = conv.iwls.GA) +
            theme(legend.position="top", legend.title=element_blank()) +
            ggtitle("a) Nearly Uniform - Multivariate Normal")


############################################################################

data <- data.gen.logit(n,d,beta,type="moderately nonuniform lev")
y <- as.matrix(data[,1])
X <- as.matrix(data[,1:d+1])

beta.iwls.trace.T1 <- iwls.logit(y,X,max.iter,delta)[,5]
conv.iwls.T1 <- length(beta.iwls.trace.T1)

beta.RandIWLS.trace.T1.unw <- slow.iwls.logit(y,X,max.iter,eps,delta,type="unweighted")[,5]
beta.RandIWLS.trace.T1.w <- slow.iwls.logit(y,X,max.iter,eps,delta,type="weighted")[,5]
beta.RandIWLS.trace.T1.infl <- slow.iwls.logit(y,X,max.iter,eps,delta,type="influence")[,5]

max.iwls.iter.T1 = max(length(beta.iwls.trace.T1), length(beta.RandIWLS.trace.T1.unw), length(beta.RandIWLS.trace.T1.w))

if (length(beta.iwls.trace.T1) < max.iwls.iter.T1){
    beta.iwls.trace.T1 <- c(beta.iwls.trace.T1, rep(beta.iwls.trace.T1[length(beta.iwls.trace.T1)], max.iwls.iter.T1-length(beta.iwls.trace.T1)))
}
if (length(beta.RandIWLS.trace.T1.unw) < max.iwls.iter.T1){
    beta.RandIWLS.trace.T1.unw <- c(beta.RandIWLS.trace.T1.unw, rep(beta.RandIWLS.trace.T1.unw[length(beta.RandIWLS.trace.T1.unw)], max.iwls.iter.T1-length(beta.RandIWLS.trace.T1.unw)))
}
if (length(beta.RandIWLS.trace.T1.w) < max.iwls.iter.T1){
    beta.RandIWLS.trace.T1.w <- c(beta.RandIWLS.trace.T1.w, rep(beta.RandIWLS.trace.T1.w[length(beta.RandIWLS.trace.T1.w)], max.iwls.iter.T1-length(beta.RandIWLS.trace.T1.w)))
}
if (length(beta.RandIWLS.trace.T1.infl) < max.iwls.iter.T1){
    beta.RandIWLS.trace.T1.infl <- c(beta.RandIWLS.trace.T1.infl, rep(beta.RandIWLS.trace.T1.infl[length(beta.RandIWLS.trace.T1.infl)], max.iwls.iter.T1-length(beta.RandIWLS.trace.T1.infl)))
}

data.to.plot.T1 <- as.data.frame(cbind(seq(1, length(beta.iwls.trace.T1)), beta.iwls.trace.T1, beta.RandIWLS.trace.T1.unw, beta.RandIWLS.trace.T1.w, beta.RandIWLS.trace.T1.infl))
colnames(data.to.plot.T1) = legend.names
plot.T1 <- plot.trace(data.to.plot.T1) +
            geom_vline(xintercept = conv.iwls.T1) +
            theme(legend.position="top", legend.title=element_blank()) +
            ggtitle('b) Moderately Nonuniform - Multivariate t (df=1)')

############################################################################

data <- data.gen.logit(n,d,beta,type="nonuniform lev")
y <- as.matrix(data[,1])
X <- as.matrix(data[,1:d+1])

beta.iwls.trace.T3 <- iwls.logit(y,X,max.iter,delta)[,5]
conv.iwls.T3 <- length(beta.iwls.trace.T3)
beta.RandIWLS.trace.T3.unw <- slow.iwls.logit(y,X,max.iter,eps,delta,type="unweighted")[,5]
beta.RandIWLS.trace.T3.w <- slow.iwls.logit(y,X,max.iter,eps,delta,type="weighted")[,5]
beta.RandIWLS.trace.T3.infl <- slow.iwls.logit(y,X,max.iter,eps,delta,type="influence")[,5]

max.iwls.iter.T3 = max(conv.iwls.T3, length(beta.RandIWLS.trace.T3.unw), length(beta.RandIWLS.trace.T3.w))

if (length(beta.iwls.trace.T3) < max.iwls.iter.T3){
    beta.iwls.trace.T3 <- c(beta.iwls.trace.T3, rep(beta.iwls.trace.T3[length(beta.iwls.trace.T3)], max.iwls.iter.T3-length(beta.iwls.trace.T3)))
}
if (length(beta.RandIWLS.trace.T3.unw) < max.iwls.iter.T3){
    beta.RandIWLS.trace.T3.unw <- c(beta.RandIWLS.trace.T3.unw, rep(beta.RandIWLS.trace.T3.unw[length(beta.RandIWLS.trace.T3.unw)], max.iwls.iter.T3-length(beta.RandIWLS.trace.T3.unw)))
}
if (length(beta.RandIWLS.trace.T3.w) < max.iwls.iter.T3){
    beta.RandIWLS.trace.T3.w <- c(beta.RandIWLS.trace.T3.w, rep(beta.RandIWLS.trace.T3.w[length(beta.RandIWLS.trace.T3.w)], max.iwls.iter.T3-length(beta.RandIWLS.trace.T3.w)))
}
if (length(beta.RandIWLS.trace.T3.infl) < max.iwls.iter.T3){
    beta.RandIWLS.trace.T3.infl <- c(beta.RandIWLS.trace.T3.infl, rep(beta.RandIWLS.trace.T3.infl[length(beta.RandIWLS.trace.T3.infl)], max.iwls.iter.T3-length(beta.RandIWLS.trace.T3.infl)))
}

data.to.plot.T3 <- as.data.frame(cbind(seq(1, length(beta.iwls.trace.T3)), beta.iwls.trace.T3, beta.RandIWLS.trace.T3.unw, beta.RandIWLS.trace.T3.w, beta.RandIWLS.trace.T3.infl))
colnames(data.to.plot.T3) = legend.names
plot.T3 <- plot.trace(data.to.plot.T3) +
            geom_vline(xintercept = conv.iwls.T3) +
            theme(legend.position="top", legend.title=element_blank()) +
            ggtitle('c) Very Nonuniform - Multivariate t (df=2)')


if(save.plots) pdf("F5_Trace_GLM.pdf")
multiplot(plot.GA, plot.T1, plot.T3, cols=1)
if(save.plots) dev.off()
