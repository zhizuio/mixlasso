#' mixlasso
#' @title Sub-function \code{tree.lasso()} for tree-lasso
#' @description
#' Main function for tree-lasso
#' 
#' @importFrom graphics axis box text mtext par image
#' @importFrom grDevices colorRampPalette dev.off grey
#' 
#' @export
tree.lasso <- function(x=x, y=y, lambda=10, tree.parm=NULL, t.idx=FALSE, num.nonpen=0, intercept=TRUE, mu=0.01, threshold=0, tol=1e-6, NoVar=50, y.mis=NULL, x.mis=NULL, t.glasso=FALSE, alpha=1, gamma=0, cov.proxy="FL", L=0, maxiter=10000, predict.re=FALSE){
  x <- data.matrix(x)
  y <- data.matrix(y)
  option <- c(maxiter, threshold, tol)
  names(option) <- c("maxiter", "threshold", "tol")
  
  if(is.null(t.idx)){
    fit <- accgrad(y, x, lambda, tree.parm$Tree, tree.parm$pre0$C, tree.parm$pre0$g_idx, tree.parm$pre0$TauNorm, mu, option, num_nonpen=num.nonpen, intercept=intercept, y_mis=y.mis, x_mis=x_mis)
    return(list(Beta=fit$Beta))
  }else{
    fit <- accgrad.re(y, x, lambda, tree.parm$Tree, t.idx, tree.parm$pre0$C, tree.parm$pre0$g_idx, tree.parm$pre0$TauNorm, mu, NoVar, option, num_nonpen=num.nonpen, intercept=intercept, y_mis=y.mis, x_mis=x.mis, t.glasso=t.glasso, alpha=alpha, gamma=gamma, cov.proxy=cov.proxy, L=L, predict.re=predict.re)
    return(list(Beta=fit$Beta, random.effects=fit$u.re, sigma=fit$sigma))
  }
  #return(list(Beta=fit$Beta))#, sigmaEU=fit$sigmaEU))
}
# function ipf.tree.lasso is the directly adapted tree-lasso method for two data sources
ipf.tree.lasso <- function(x, y, p=NULL, h=0.7, lambda=rep(10,2), tree.parm, num.nonpen=0, intercept=TRUE, mu=0.01, threshold=0){
  X <- data.matrix(x)
  Y <- data.matrix(y)
  m <- dim(Y)[2]
  option <- c(10000, threshold, 1e-6)
  names(option) <- c("maxiter", "threshold", "tol")
  
  pre0 <- pre.grad2(tree.parm$Tree, tree.parm$Tw, lambda)
  
  acc <- accgrad2(Y, X, lambda, tree.parm$Tree, pre0$C, pre0$g_idx, pre0$TauNorm, p, mu, option)
  Beta <- acc$Beta
  obj <- acc$obj
  time <- acc$time
  
  return(Beta)
}
# the function vertical.image.legend is orginally from the R package "aqfig"
vertical.image.legend <- function (zlim, col, legend.cex.axis=1) {
  starting.par.settings <- par(no.readonly = TRUE)
  #on.exit(par(starting.par.settings))
  mai <- par("mai")
  fin <- par("fin")
  x.legend.fig <- c(1 - (mai[4]/fin[1]), 1)
  y.legend.fig <- c(mai[1]/fin[2], 1 - (mai[3]/fin[2]))
  x.legend.plt <- c(x.legend.fig[1] + (0.18 * (x.legend.fig[2] - 
                                                 x.legend.fig[1])), x.legend.fig[2] - (0.6 * (x.legend.fig[2] - 
                                                                                                x.legend.fig[1])))
  y.legend.plt <- y.legend.fig
  cut.pts <- seq(zlim[1], zlim[2], length = length(col) + 1)
  z <- (cut.pts[1:length(col)] + cut.pts[2:(length(col) + 1)])/2
  par(new = TRUE, pty = "m", plt = c(x.legend.plt, y.legend.plt))
  # If z is not inincreasing, only two values
  if(all(diff(z) > 0)){
    image(x = 1.5, y = z, z = matrix(z, nrow = 1, ncol = length(col)), 
          col = col, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    axis(4, mgp = c(3, 0.2, 0), las = 2, cex.axis = legend.cex.axis, tcl = -0.1)
    box()
  }
  mfg.settings <- par()$mfg
  par(starting.par.settings)
  par(mfg = mfg.settings, new = FALSE)
  
}



