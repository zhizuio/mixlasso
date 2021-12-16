#' mixlasso
#' @title Structured penalized regression
#' @description
#' Function producing results of the structured penalized regression
#' 
#' @importFrom Matrix Diagonal bdiag
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats sd
#' @importFrom parallel detectCores makeCluster	stopCluster clusterEvalQ
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom penalized penalized
#' 
#' @param x,y \code{x} is the input design matrix; \code{y} is the input response matrix
#' @param x_test,y_test \code{x} is the input validated design matrix; \code{y} is the input validated response matrix
#' @param p the number of predictors from different data source.
#' @param foldid an vector of values for the cross-validation.
#' @param num.nonpen number of predictors forced to be estimated (i.e., nonpenalization).
#' @param method specify the the method to optimize its penalty parameters. The penalty parameters of \code{elastic-net}, \code{IPF-lasso}, \code{sIPF-elastic-net}, \code{IPF-elastic-net}, \code{IPF-tree-lasso} and \code{clogitLasso} are optimzed by the EPSGO algorithm. The penalty parameter of \code{lasso} and \code{tree-lasso} is optimzed by cross-validation. The default method is \code{IPF-lasso} 
#' @param lambda optional user-supplied \code{lambda} sequence; default is NULL, and \code{espsgo} chooses its own sequence except the tree-lasso methods.
#' @param bounds bounds for the interval-searching parameters
#' @param strata.surv stratification variable for the Cox survival model.
#' @param search.path save the visited points, default is \code{FALSE}.
#' @param EI.eps he convergence threshold for the expected improvement between fmin and the updated point 
#' @param fminlower minimal value for the function Q.func, default is 0.
#' @param threshold threshold for estimated coefficients of the tree-lasso methods.
#' @param N define the number of start points depending on the dimensionality of the parameter space.
#' @param min.iter the minimus iterations after the initial \code{N} iterations.
#' @param seed random seed.
#' @param parallel If \code{TRUE}, use parallel foreach to fit each fold except parallelizing each lambda for the tree-lasso methods. If \code{c(TRUE,TRUE)}, use parallel foreach to fit each fold and each lambda. 
#' @param verbose print the middle search information, default is \code{TRUE}.
#' ##param lib.loc a character vector describing the location of R library trees to search through, or NULL by default.
#' @return An object of list "\code{mixlasso}" is returned:
#'  \item{cvm}{the mean cross-validated error}  
#'  \item{cvm_cv}{the mean cross-validated error if providing external dataset "\code{x_test}" and "\code{y_test}". } 
#'  \item{alpha}{optimized \code{alpha}}  
#'  \item{lambda}{optimized \code{lambda}}  
#'  \item{pred}{the prediction of the responses}  
#'  \item{ipf}{optimzed penalty factors}  
#'  \item{Beta}{estimate of the coefficients}  
#'  \item{cv}{number of nonzero coefficients}  
#' 
#' @references Zhao, Z. & Zucknick, M. (2020). \emph{Stuctured penalized regression for drug sensitivity prediction.} JRSSC.
#' @export
mixlasso <- function(x, y, z=NULL, x_test=NULL, y_test=NULL, z_test=NULL, p=NA, foldid=NULL, num.nonpen=0, method="IPF-lasso", family="gaussian", nfolds=5,
                                y.mis=NULL, y.mis_test=NULL, x.mis=NULL, x.mis_test=NULL, tree.parm=NULL, cv.measure="mse", type.measure="deviance", type.min="lambda.min", standardize.response=FALSE,
                        lambda=NULL, bounds=NULL, bound.scale=NA, strata.surv=NULL, search.path=FALSE, EI.eps=0.01, fminlower=0, intercept=TRUE,
                        threshold=0, tol=1e-6, mu=0.01, NoVar=50, N=NULL, min.iter=20, seed=1234,parallel=FALSE, verbose=TRUE,t.idx=NULL,t.idx_test=NULL,
                        t.glasso=FALSE,alpha=1,gamma=0,
                        # the following parameters have problem to pass into tune.tree.re.interval
                        maxiter=10000,cov.proxy="FL",predict.re=FALSE,...){
 
  if((method!="lasso") & (method!="tree-lasso") & is.null(bounds)){
    if(method=="elastic-net"){
      bounds <- t(data.frame(alpha=c(0,1)))
    }else{
      # the default setting for bounds is only for two data sources
      if(method=="IPF-lasso" | method=="IPF-cox"){
        bounds <- t(data.frame(ipf=c(0.1,10)))
      }else{
        if(method=="sIPF-elastic-net"){
          bounds <- t(data.frame(alpha=c(0,1), ipf1 = c(0.1, 10)))
        }else{
          if(method=="IPF-elastic-net"){
            bounds <- t(data.frame(alpha1=c(0,1), alpha2=c(0,1), ipf1 = c(0.5, 8)))
          }else{
            if(method=="IPF-tree-lasso"){
              bounds <- t(data.frame(tree.ipf1 = c(0.1, 5)))
            }else{
              stop("Please give searching bounds for EPSGO algorithm!")
            }
          }
        } 
      }
    } 
    colnames(bounds) <- c("lower", "upper")
  }
  
  fit <- NULL
  #=============
  # IPF-lasso, sIPFEN
  #=============
  if((method=="IPF-lasso") | (method=="sIPF-elastic-net")){
    
    fit0 <- epsgo(Q.func = "tune.glmnet.interval", bounds = bounds, lambda=lambda, N=N,
                 parms.coding = "none", seed = seed, fminlower = fminlower, x = x, y = y, num.nonpen=num.nonpen, bound.scale=bound.scale,
                 family = family, foldid = foldid, nfolds = nfolds, type.min = type.min, p=p, intercept=intercept, standardize.response=standardize.response,
                 type.measure = type.measure, min.iter=min.iter, verbose=verbose, parallel=parallel, EI.eps=EI.eps, ...)
    sumint <- summary(fit0, verbose=F)
    mse_cv <- sumint$opt.error
    alpha <- sumint$opt.alpha
    lambda <- sumint$opt.lambda
    ipf <- sumint$opt.ipf
    
    # prediction
    if(length(p) == 1){
      adpen <- rep(1, p)
    }else{
      adpen <- c(rep(1,p[1]), rep(ipf,p[2]))
    }
    fit <- glmnet(x=x,y=y, family=family,
                          alpha=alpha,
                          offset = NULL,
                          lambda = seq(lambda*0.8,lambda*1.2,length=11),
                          penalty.factor=adpen,
                          intercept=intercept,
                          standardize.response=F)
    
    ypred <- cbind(rep(1,dim(x_test)[1]),x_test)%*%matrix(c(fit$a0[ceiling(11/2)], fit$beta[,ceiling(11/2)]),ncol=1)
    mse_val <- sum((y_test - ypred)^2)/prod(dim(y_test))
    vs <- sum(fit$beta[,ceiling(11/2)]!=0)
  }
  
  #=============
  # IPF-Cox
  #=============
  if(method=="IPF-cox"){
    
    fit0 <- epsgo(Q.func = "tune.surv.interval", strata.surv=strata.surv, bounds = bounds, lambda=lambda, alpha=alpha, N=N,
                 parms.coding = "none", seed = seed, fminlower = fminlower, x = x, y = y, num.nonpen=num.nonpen, bound.scale=bound.scale,
                 family = "cox", foldid = foldid, nfolds = nfolds, type.min = type.min, p=p, intercept=intercept, standardize.response=standardize.response,
                 type.measure = type.measure, min.iter=min.iter, verbose=verbose, parallel=parallel, EI.eps=EI.eps, ...)
    sumint <- summary(fit0, verbose=F)
    mse_cv <- sumint$opt.error
    alpha <- sumint$opt.alpha
    lambda <- sumint$opt.lambda
    ipf <- sumint$opt.ipf
    
    # fit the optimal model
    adpen <- c(rep(0,num.nonpen), rep(1, sum(p)))
    if(length(p) > 1){
      for(i in 1:(length(p)-1)){
        adpen[num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]] <- ipf[i]
      }
    }
    
    if(is.null(strata.surv)){
      y2 <- y
    }else{
      y2 <- stratifySurv(y, strata.surv)
    }
    fit <- glmnet(x=x,y=y2, family="cox",
                   alpha=alpha,
                   offset = NULL,
                   lambda = seq(lambda*0.8,lambda*1.2,length=11),
                   penalty.factor=adpen,
                   intercept=intercept)
    
    ypred <- NULL
    mse_val <- NULL
    vs <- NULL
  }
  
  #==================
  # Tree-lasso
  #==================
  
  if(method=="tree-lasso"){
    if( sum(is.na(y)) | sum(is.na(x)) | sum(is.na(y_test)) | sum(is.na(x_test)) )
      stop("All data sets have to be complete without missing values!")
    if(is.null(lambda)){
      cat("Warning: Please provide a proper lambda sequence!")
      lambda <- seq(2,5,length=10)
      # fun.lambda <- function(x,y){
      #   sum(matrix(x,nrow=1)%*%(y-matrix(rep(apply(y,2,mean),each=dim(y)[1]),ncol=dim(y)[2])))
      # }
      # lambda.max <- sqrt(max(sapply(split(x, rep(1:ncol(x), each=nrow(x))), fun.lambda, y=y)))
      # lambda <- seq(lambda.max*0.1, lambda.max, length=5)
    }
    
    cvm0 <- numeric(length(lambda))
    cv5<-function(xx,la) {sum((y[foldid==xx,] - cbind(rep(1,sum(foldid==xx)),x[foldid==xx,]) %*% tree.lasso(x=x[!foldid==xx,], y=y[!foldid==xx,],lambda=la,tree.parm=tree.parm,num.nonpen=num.nonpen, threshold=threshold, maxiter=maxiter, tol=tol, mu=mu, t.idx=t.idx)$Beta)^2)/(sum(foldid==xx)*(dim(y)[2]))}
    la.seq<-function(la) {mean(sapply(1:max(foldid), cv5, la))}
    la.xx <- cbind(rep(lambda,times=max(foldid)), rep(1:max(foldid),each=length(lambda)))
    cvm0 <- numeric(nrow(la.xx))
    
    if(parallel){
      cores <- length(lambda) * max(foldid)
      cl <- makeCluster(cores)
      #clusterEvalQ(cl, library(mixlasso))
      registerDoParallel(cl)
      cvm0[1:cores] <- foreach(i = 1:cores, .combine=c, .packages= c('base','Matrix','MASS')) %dopar%{
        cv5(la.xx[i,2], la.xx[i,1])
      }
      stopCluster(cl)
    }else{
      #lambda <- matrix(lambda,ncol=1)
      #cvm0<-apply(lambda, 1, la.seq)
      for(i in 1:nrow(la.xx)){
        cvm0[i] <- cv5(la.xx[i,2], la.xx[i,1])
      }
    }
    
    cvm0 <- rowMeans(matrix(cvm0,ncol=max(foldid)))
    mse_cv0 <- min(cvm0)
    alpha <- 1
    lambda <- lambda[which.min(cvm0)]
    
    # prediction
    Beta <- tree.lasso(x=x, y=y, lambda=lambda, tree.parm=tree.parm, num.nonpen=num.nonpen, threshold=threshold, maxiter=maxiter, tol=tol, mu=mu, t.idx=t.idx)$Beta
    mse_cv <- sum((y - cbind(rep(1,dim(x)[1]),x)%*%Beta)^2)/prod(dim(y))
    
    ypred <- cbind(rep(1,dim(x_test)[1]),x_test)%*%Beta
    mse_val <- sum((y_test - ypred)^2)/prod(dim(y_test))
    vs <- sum(Beta[-1,]!=0)
    ipf <- 1
    
  }
  
  #==================
  # IPF-tree-lasso
  #==================
  
  if(method=="IPF-tree-lasso"){
    if( sum(is.na(y)) | sum(is.na(x)) | sum(is.na(y_test)) | sum(is.na(x_test)) )
      stop("All data sets have to be complete without missing values!")
    if(is.null(lambda)){
      warning("Please provide a proper lambda sequence!")
      lambda <- seq(2,25,length=10)
      # fun.lambda <- function(x,y){
      #   sum(matrix(x,nrow=1)%*%(y-matrix(rep(apply(y,2,mean),each=dim(y)[1]),ncol=dim(y)[2])))
      # }
      # lambda.max <- sqrt(max(sapply(split(x, rep(1:ncol(x), each=nrow(x))), fun.lambda, y=y)))
      # lambda <- seq(lambda.max*0.1, lambda.max, length=5)
    }
    fit0 <- epsgo(Q.func = "tune.tree.interval", bounds = bounds, lambda=lambda, N=N,
                 parms.coding = "none", seed = seed, fminlower = fminlower, x = x, y = y, 
                 intercept=intercept, foldid = foldid, p=p, standardize.response=F, num.nonpen=num.nonpen,
                 min.iter=min.iter, verbose=verbose, parallel=parallel, EI.eps=EI.eps, threshold=threshold, ...)
    sumint <- summary(fit, verbose=F)
    
    ipf <- sumint$opt.ipf
    mse_cv <- sumint$opt.error
    alpha <- 1
    lambda <- sumint$opt.lambda
    
    # prediction
    Xtemp <- x
    for(i in 1:(length(p)-1)) Xtemp[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]] <- Xtemp[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]]/ipf[i]
    Beta <- tree.lasso(x=Xtemp, y=y, lambda=lambda, tree.parm=tree.parm, num.nonpen=num.nonpen, threshold=threshold, maxiter=maxiter)$Beta
    pf.r <- c(1, ipf)
    for(s in 1:length(p)) Beta[num.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s]+1,] <- 1/pf.r[s]*Beta[num.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s]+1,]  
    ypred <- data.matrix(cbind(rep(1,dim(x_test)[1]),x_test))%*%Beta
    mse_val <- sum((y_test - ypred)^2)/prod(dim(y_test))
    vs <- sum(Beta[-1,]!=0)
  }
  
  #==================
  # IPF-tree-lasso-re
  #==================
  if(method=="IPF-tree-lasso-re"){
    if( sum(is.na(x)) | sum(is.na(x_test)) )
      stop("X data have to be complete without missing values!")
    
    if(is.null(y.mis)){
      y.mis <- matrix(0, nrow=dim(y)[1], ncol=dim(y)[2])
    }else{
      y[y.mis==1] <- 0
    }
    if(is.null(y.mis_test)){
      y.mis_test <- matrix(0, nrow=dim(y_test)[1], ncol=dim(y_test)[2])
    }else{
      y_test[y.mis_test==1] <- 0
    }
    
    if(!is.null(x.mis)) x[x.mis==1] <- 0
    if(!is.null(x.mis_test)) x_test[x.mis_test==1] <- 0
    # 
    # if(is.null(lambda)){
    #   warning("Please provide a proper lambda sequence!")
    #   lambda <- seq(2,5,length=10)
    # }
    fit0 <- epsgo(Q.func = "tune.tree.re.interval", bounds = bounds, bound.scale=bound.scale, alpha=alpha, N=N, tree.parm=tree.parm,
                 parms.coding = "none", seed = seed, fminlower = fminlower, x = x, y = y, z = z,
                 intercept=intercept, foldid = foldid, cv.measure=cv.measure, p=p, standardize.response=F, num.nonpen=num.nonpen,
                 min.iter=min.iter, verbose=verbose, parallel=parallel, EI.eps=EI.eps, threshold=threshold, 
                 tol=tol,mu=mu,NoVar=NoVar,y.mis=y.mis,x.mis=x.mis,t.idx=t.idx,t.glasso=t.glasso,maxiter=maxiter,cov.proxy=cov.proxy,predict.re=predict.re, ...)
    sumint <- summary(fit0, verbose=F)
    
    mse_cv0 <- sumint$opt.error
    lambda <- sumint$opt.lambda
    gamma <- sumint$opt.gamma
    alpha <- sumint$opt.alpha
    ipf <- as.numeric(sumint$opt.ipf)
    
    Xtemp <- x
    if(length(p) > 1){
      for(i in 1:(length(p)-1)) Xtemp[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]] <- Xtemp[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]]/ipf[i]
      #x <- Xtemp
    }
    
    # prediction
    fit <- tree.lasso(x=Xtemp, y=y, lambda=lambda, tree.parm=tree.parm, num.nonpen=num.nonpen,threshold=threshold,maxiter=maxiter,tol=tol,mu=mu,NoVar=NoVar,y.mis=y.mis,x.mis=x.mis,t.idx=t.idx,t.glasso=t.glasso,alpha=alpha,gamma=gamma,predict.re=predict.re)
    
    if(length(p) > 1){
      pf.r <- c(1, ipf)
      for(s in 1:length(p)) fit$Beta[num.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s]+1,] <- 1/pf.r[s]*fit$Beta[num.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s]+1,]  
    }
    
    if(is.null(t.idx)){
      #if(cv.measure == "mse")
        mse_cv <- sum(((y - crossprod(t(cbind(rep(1,dim(x)[1]),x)), fit$Beta))[y.mis!=1])^2)/sum(y.mis!=1)
      #if(cv.measure == "spearman")
      #  mse_cv <- mean(cor(y, crossprod(t(cbind(rep(1,dim(x)[1]),x)), fit$Beta), method="spearman" )[1,])
      ypred <- crossprod(t(cbind(rep(1,dim(x_test)[1]),x_test)), fit$Beta)
    }else{
      #if(cv.measure == "mse")
        #mse_cv <- sum((y - Matrix::crossprod(t(x.star), fit$Beta) - z %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y)[2]) )[y.mis!=1]^2)/sum(y.mis!=1)#prod(dim(y))
      #if(cv.measure == "spearman")
      #  mse_cv <- mean(cor(y, data.matrix(Matrix::crossprod(t(x.star), fit$Beta) + z %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y)[2])), method="spearman" )[1,] )
        ypred0 <- matrix(NA, nrow=nrow(y), ncol=ncol(y))
        for(i in unique(t.idx)){
          if(intercept){
            #mse_cv <- mse_cv + sum((y[t.idx==i,] - Matrix::crossprod(t(cbind(rep(1,sum(t.idx==i)),x[t.idx==i,])), fit$Beta[(i-1)*(sum(p)+1)+1:(sum(p)+1),])
            #                        - z[t.idx==i,] %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y)[2]) )[y.mis[t.idx==i,]!=1]^2)
            x.tmp <- cbind(rep(1,sum(t.idx==i)),matrix(x[t.idx==i,], ncol=ncol(x)))
            ypred0[t.idx==i,] <- data.matrix(Matrix::crossprod(t(x.tmp), fit$Beta[(i-1)*(sum(p)+1)+1:(sum(p)+1),])
                                    + z[t.idx==i,] %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y)[2]))
            
          }else{
            #mse_cv <- mse_cv + sum((y[t.idx==i,] - Matrix::crossprod(t(x[t.idx==i,]), fit$Beta[(i-1)*sum(p)+1:sum(p),])
            #                        - z[t.idx==i,] %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y)[2]) )[y.mis[t.idx==i,]!=1]^2)
            x.tmp <- matrix(x[t.idx==i,], ncol=ncol(x))
            ypred0[t.idx==i,] <- data.matrix(Matrix::crossprod(t(x.tmp), fit$Beta[(i-1)*sum(p)+1:sum(p),])
                                    + z[t.idx==i,] %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y)[2]))
          }
        }
        #mse_cv <- mse_cv/sum(y.mis!=1)
        mse_cv <- sum(((y - ypred0)[y.mis==0])^2)/sum(y.mis==0)
        mse_cv_complete <- sum((y[rowSums(y.mis)==0,] - ypred0[rowSums(y.mis)==0,])^2)/(sum(rowSums(y.mis)==0)*ncol(y))
        
        ypred <- matrix(NA, nrow=nrow(y_test), ncol=ncol(y_test))
        for(i in unique(t.idx_test)){
          if(intercept){
            x.tmp <- cbind(rep(1,sum(t.idx_test==i)),matrix(x_test[t.idx_test==i,], ncol=ncol(x)))
            ypred[t.idx_test==i,] <- data.matrix(Matrix::crossprod(t(x.tmp), fit$Beta[(i-1)*(sum(p)+1)+1:(sum(p)+1),])
                                     + z_test[t.idx_test==i,] %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y_test)[2]))
          }else{
            x.tmp <- matrix(x_test[t.idx_test==i,], ncol=ncol(x))
            ypred[t.idx_test==i,] <- data.matrix(Matrix::crossprod(t(x.tmp), fit$Beta[(i-1)*sum(p)+1:sum(p),]) 
                                     + z_test[t.idx_test==i,] %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y_test)[2]))
          }
        }
    }
    #if(cv.measure == "mse")
      mse_val <- sum(((y_test - ypred)[y.mis_test==0])^2)/sum(y.mis_test==0)#prod(dim(y_test))
      mse_val_complete <- sum((y_test[rowSums(y.mis_test)==0,] - ypred[rowSums(y.mis_test)==0,])^2)/(sum(rowSums(y.mis_test)==0)*ncol(y_test))
      
    #if(cv.measure == "spearman")
    #  mse_val <- mean( cor(y_test, ypred, method="spearman" )[1,] )
    vs <- sum(fit$Beta!=0)
  }
  
    if(method == "tree-lasso" ){
      return(list(cvm=mse_val, cvm_cv=mse_cv, cvm_cv0=mse_cv0, gamma=gamma, alpha=alpha, lambda=lambda, pred=ypred, ipf=ipf, Beta=Beta, cv=vs))
    }else{
      if(method == "IPF-tree-lasso-re") {
        return(list(cvm=mse_val, cvm_cv=mse_cv, cvm_cv0=mse_cv0, cvm_complete=mse_val_complete, cvm_cv_complete=mse_cv_complete, gamma=gamma, alpha=alpha, lambda=lambda, pred=ypred, ipf=ipf, Beta=fit$Beta, cv=vs, random.effects=fit$random.effects))
      }else{
        return(list(cvm=mse_val, cvm_cv=mse_cv, alpha=alpha, lambda=lambda, pred=ypred, ipf=ipf, Beta=c(fit$a0[ceiling(11/2)], fit$beta[,ceiling(11/2)]), cv=vs, fit=fit, search.fit=fit0))
      }
    }
}
