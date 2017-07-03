## Author: Robert Lange
## DATE: 06/2017
## BGSE Master Project- RandNLA for GLMS with Big Datasets

# CODE VISUALIZES THE TRANSFORM DONE IN HDX - FJLT PROCESSING HDX ESSENTIALLY
# UNIFORMIZES THE LEVERAGE SCORES - GEOMETRICAL INTERPRETATION
set.seed(101)

HD_TRAFO <- function(X) {
    n <- nrow(X)
    # Create Hadamard Transform
    H <- hadamard.matrix(n)
    # Create spiky preprocessor
    D <- matrix(0, n, n)
    diag(D) <- sample(c(-1, 1), size = n, replace = T)

    HD <- H %*% D
    HDX <- HD %*% X
    HDX
}

plot_points <- function(X, HDX){
    data.points <- data.frame(cbind(X, HDX))

    plot.points <- ggplot(data.points) +
      geom_point(aes(x=X1, y=X2), colour="red", alpha = 0.25) +
      geom_point(aes(x=X3, y=X4), colour="blue",  alpha = 0.25) +
      ylab(TeX('$X_2$ and $HDX_2$')) +
      xlab(TeX('$X_1$ and $HDX_1$')) +
      theme(plot.title = element_text(hjust = 0.5), legend.position="none")
    plot.points
}

plot_lev <- function(lev, HDX_lev){
    data.lev <- data.frame(cbind(lev, HDX_lev))

    plot.lev <- ggplot(data.lev)+
      geom_point(aes(x=leverage, y=leverage.1), colour="black") +
      ylab(TeX('Leverage Scores $HDX$')) +
      xlab(TeX('Leverage Scores $X$')) +
      theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
      xlim(0, 1) +
      ylim(0, 0.1)
    plot.lev
}

plot_lev_dens <- function(lev, HDX_lev){
    data.lev <- data.frame(cbind(lev, HDX_lev))

    plot.lev.dens <- ggplot(data.lev)+
    geom_density(aes(x=leverage), colour="red", alpha = 0.25) +
    geom_density(aes(x=leverage.1), colour="blue", alpha = 0.25) +
    ylab("Kernel Density") +
    xlab("Leverage Score") +
    theme(plot.title = element_text(hjust = 0.5, lineheight=.8)) +
    xlim(0, 1) +
    ylim(0, 1)
    plot.lev.dens
}

X_GA <- data_lev_gen(n, d, type="nearly uniform lev")
lev_GA <- LEV(X_GA)
HDX_GA <- HD_TRAFO(X_GA)
lev_HDX_GA <- LEV(HDX_GA)

X_T1 <- data_lev_gen(n, d, type="moderately nonuniform lev")
lev_T1 <- LEV(X_T1)
HDX_T1 <- HD_TRAFO(X_T1)
lev_HDX_T1 <- LEV(HDX_T1)

X_T3 <- data_lev_gen(n, d, type="very nonuniform lev")
lev_T3 <- LEV(X_T3)
HDX_T3 <- HD_TRAFO(X_T3)
lev_HDX_T3 <- LEV(HDX_T3)

p1 <- plot_points(X=X_GA, HDX=HDX_GA) + ggtitle("a) Nearly Uniform \n  Multivariate Normal")
p3 <- plot_points(X=X_T1, HDX=HDX_T1) + ggtitle("b) Moderately Nonunif. \n Multivariate t (df=1)")
p5 <- plot_points(X=X_T3, HDX=HDX_T3) + ggtitle("c) Very Nonuniform \n Multivariate t (df=3)")

p2a <- plot_lev(lev=lev_GA, HDX=lev_HDX_GA) + ggtitle("a) Leverage Scores \n HDX and X")
p2b <- plot_lev_dens(lev=lev_GA, HDX=lev_HDX_GA) + ggtitle("a) Density Leverage \n X (Red) - HDX (Blue)")

p4a <- plot_lev(lev=lev_T1, HDX=lev_HDX_T1) + ggtitle("b) Leverage Scores \n HDX and X")
p4b <- plot_lev_dens(lev=lev_T1, HDX=lev_HDX_T1) + ggtitle("b) Density Leverage \n X (Red) - HDX (Blue)")

p6a <- plot_lev(lev=lev_T3, HDX=lev_HDX_T3) + ggtitle("c) Leverage Scores \n HDX and X")
p6b <- plot_lev_dens(lev=lev_T3, HDX=lev_HDX_T3) + ggtitle("c) Density Leverage \n X (Red) - HDX (Blue)")

if(save.plots) pdf("F1_HDX_Transform.pdf")
multiplot(p1, p3, p5, p2a, p4a, p6a, p2b, p4b, p6b, cols=3)
if(save.plots) dev.off()
