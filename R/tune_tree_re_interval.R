#' mixlasso
#' @title Wrapper function for tree-lasso objects.
#' @description
#' Wrapper function for tree-lasso objects used by epsgo function. This function is mainly used within the function \code{epsgo}.
#' 
#' @importFrom parallel detectCores makeCluster	stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' 
#' @param parms tuning parameter alpha for the tree-lasso object.
#' @param x,y \code{x} is a matrix where each row refers to a sample a each column refers to a gene; \code{y} is a factor which includes the class for each sample
#' @param lambda A user supplied lambda sequence. 
#' @param nfolds number of cross-validation's folds, default 5.
#' @param foldid an optional vector of values between 1 and nfold identifying what fold each observation is in. If supplied, nfold can be missing.
#' @param num.nonpen number of predictors forced to be estimated (i.e., nonpenalization).
#' @param seed random seed
#' @param intercept should  intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE}).
#' @param standardize.response standardization for the response variables. Default: \code{TRUE}.
#' @param p the number of predictors from different data source.
#' @param verbose print the middle search information, default is \code{TRUE}.
#' @param parallel  If \code{TRUE}, use parallel foreach to fit each lambda. If \code{c(TRUE,TRUE)}, use parallel foreach to fit each lambda and each fold. 
#' @param search.path save the visited points, default is \code{FALSE}.
#' @param threshold threshold for estimated coefficients of the tree-lasso models.
#' @export
tune.tree.re.interval<-function(parms, x, y, z=z,
                             #lambda = NULL, 
                             alpha = 1,
                             nfolds = 5,
                             foldid=NULL,
                             cv.measure=cv.measure,
                             tree.parm=tree.parm,
                             num.nonpen = 0,
                             seed=12345, 
                             intercept=TRUE,
                             standardize.response=FALSE,
                             p=NULL,
                             verbose=FALSE,
                             parallel=FALSE,
                             search.path=FALSE,
                             threshold=threshold,
                             tol=tol,
                             mu=mu,
			                       NoVar=NoVar,
                             y.mis=y.mis,
                             x.mis=x.mis,
                             t.idx=t.idx,
                             t.glasso=t.glasso,
                             maxiter=10000,
                             cov.proxy="FL",
                             L=0,
			                       predict.re=predict.re,
                             ...){
  
  # 1. decode the parameters ############################################################
  if(!identical(grep("lambda", names(parms)), integer(0))){
    lambda <- parms[grep("lambda", names(parms))]
    if(verbose) print(paste("lambda=",paste(as.character(lambda),collapse=",")))
  # }else{
  #   stop("Please give searching bounds for lambda!")
  }
  if(!identical(grep("gamma", names(parms)), integer(0))){
    gamma <- parms[grep("gamma", names(parms))]
    if(verbose) print(paste("gamma=",paste(as.character(gamma),collapse=",")))
  }else{
    gamma <- 0
    names(gamma) <- "gamma"
  }
  if(!identical(grep("ipf", names(parms)), integer(0))){
    ipf <- parms[grep("ipf", names(parms))]
    if(verbose) print(paste("IPF=",paste(as.character(ipf),collapse=":")))
  }else{
    ipf<-NA
  }
  
  if(t.glasso & (alpha[1]==1)){
    alpha <- matrix(seq(.1, .9, length=5),ncol=1)
  }
  # else{
  #   alpha <- 0
  # }
  if( (length(alpha)>1) & (!is.matrix(alpha)) ) alpha <- matrix(alpha, ncol=1)
  opt.alpha <- alpha[1]
  
  if(standardize.response) y<-scale(y)[,]
  
  #browser()
  #=========
  # using augmented data
  #=========
  #if(is.null(lambda)) stop("No given lambda sequence!")
  #lambda <- matrix(lambda,ncol=1)
  if(length(p) > 1){
    x2 <- x
    for(i in 1:(length(p)-1)) 
      x2[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]] <- x2[,num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]]/parms[i]
    x <- x2
  }
  
  #tree.parm <- tree.parms(y=y)
  
  #cvm0 <- numeric(length(lambda))
  cvm0 <- numeric(length(alpha))
  cvk<-function(xx,al){#,alpha=1){
    #browser()
    if(is.null(t.idx)){
      fit <- tree.lasso(x=x[!foldid==xx,], y=y[!foldid==xx,],lambda=lambda,tree.parm=tree.parm,num.nonpen=num.nonpen, 
                        threshold=threshold,tol=tol,mu=mu,NoVar=NoVar,y.mis=y.mis[!foldid==xx,],x.mis=x.mis[!foldid==xx,])
      cvm0 <- sum((y[foldid==xx,] - crossprod(t(cbind(rep(1,sum(foldid==xx)),x[foldid==xx,])), fit$Beta))[y.mis[foldid==xx,]!=1]^2)/sum(y.mis[foldid==xx,]!=1)#(sum(foldid==xx)*(dim(y)[2]))
      if(cv.measure == "aic")
        cvm0 <- cvm0 + 2*(sum(fit$Beta!=0)-(num.nonpen+1)*ncol(y))
    }else{
      if( L[1]==0 ){
        fit <- tree.lasso(x=x[!foldid==xx,], y=y[!foldid==xx,],lambda=lambda,tree.parm=tree.parm,num.nonpen=num.nonpen, 
                          threshold=threshold,tol=tol,mu=mu,NoVar=NoVar,y.mis=y.mis[!foldid==xx,],x.mis=x.mis[!foldid==xx,],t.idx=t.idx[!foldid==xx],t.glasso=t.glasso,alpha=al,gamma=gamma,cov.proxy=cov.proxy,L=0,maxiter=maxiter,predict.re=predict.re)
      }else{
        fit <- tree.lasso(x=x[!foldid==xx,], y=y[!foldid==xx,],lambda=lambda,tree.parm=tree.parm,num.nonpen=num.nonpen, 
                          threshold=threshold,tol=tol,mu=mu,NoVar=NoVar,y.mis=y.mis[!foldid==xx,],x.mis=x.mis[!foldid==xx,],t.idx=t.idx[!foldid==xx],t.glasso=t.glasso,alpha=al,gamma=gamma,cov.proxy=cov.proxy,L=L[,xx],maxiter=maxiter,predict.re=predict.re)
      }
      
      cvm0 <- 0
      if(intercept){
        for(ii in unique(t.idx)){
          y.tmp <- matrix(y[foldid==xx & t.idx==ii,], ncol=ncol(y))
          x.tmp <- matrix(x[foldid==xx & t.idx==ii,], ncol=ncol(x))
          cvm0 <- cvm0 + sum((y.tmp - Matrix::crossprod(t(cbind(rep(1,sum(foldid==xx & t.idx==ii)),x.tmp)), fit$Beta[(ii-1)*(sum(p)+1)+1:(sum(p)+1),])
                              - z[foldid==xx & t.idx==ii,] %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y)[2]) )[y.mis[foldid==xx & t.idx==ii,]!=1]^2)
        }
        if(cv.measure == "aic")
          cvm0 <- cvm0 + 2*(sum(fit$Beta!=0)-(num.nonpen+1)*ncol(y))
      }else{
        for(ii in unique(t.idx)){
          y.tmp <- matrix(y[foldid==xx & t.idx==ii,], ncol=ncol(y))
          x.tmp <- matrix(x[foldid==xx & t.idx==ii,], ncol=ncol(x))
          cvm0 <- cvm0 + sum((y.tmp - Matrix::crossprod(t(x.tmp), fit$Beta[(ii-1)*sum(p)+1:sum(p),])
                              - z[foldid==xx & t.idx==ii,] %*% matrix(fit$random.effects,ncol=1) %*% matrix(1,nrow=1,ncol=dim(y)[2]) )[y.mis[foldid==xx & t.idx==ii,]!=1]^2)
        }
        if(cv.measure == "aic")
          cvm0 <- cvm0 + 2*(sum(fit$Beta!=0)-num.nonpen*ncol(y))
      }
      ##if( sum(abs(fit$Beta)>threshold) <= max(t.idx)*(ifelse(intercept,1,0)+num.nonpen)*ncol(y) )
      ##  cvm0 <- 99999999#Inf
      cvm0 <- cvm0/sum(y.mis[foldid==xx,]!=1)
    }
    
    return(cvm0)
  }
  al.seq<-function(al) {mean(sapply(1:max(foldid), cvk, al))}
  #al.xx <- cbind(rep(lambda,times=max(foldid)), rep(1:max(foldid),each=length(lambda)))
  al.xx <- cbind(rep(alpha,times=max(foldid)), rep(1:max(foldid),each=length(alpha)))
  #browser()
  if(sum(parallel)==2){
    cvm0 <- rep(0, nrow(al.xx))
    cores <- length(cvm0)
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    cvm0[1:cores] <- foreach(i = 1:cores, .combine='c', .packages= c('base','MASS','Matrix')) %dopar%{
        cvk(al.xx[i,2], al.xx[i,1])#, alpha)
    }
    cvm0 <- colMeans( matrix(cvm0, ncol=length(alpha), byrow=TRUE) )
    if( length(alpha)>1 )
      opt.alpha <- alpha[which.min(cvm0)]
    q.val <- min(cvm0)
    stopCluster(cl)
  }else{
    if(sum(parallel)==1){
      #cores <- length(lambda)
      cores <- max(foldid)
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      if(length(alpha)>1){
        cvm0 <- rep(0,length(alpha))
        for(idx.al in 1:length(alpha)){
          cvm0[idx.al] <- mean( foreach(i = 1:cores, .combine='c', .packages= c('base','MASS','Matrix')) %dopar%{cvk(i, alpha[idx.al])} )
        }
        opt.alpha <- alpha[which.min(cvm0)]
        q.val <- min(cvm0)
      }else{
        q.val <- mean( foreach(i = 1:cores, .combine='c', .packages= c('base','MASS','Matrix')) %dopar%{cvk(i, alpha)} )
      }
      stopCluster(cl)
    }else{
      for(i in 1:nrow(al.xx))
        cvm0[i] <-  cvk(al.xx[i,2], al.xx[i,1])
      
      cvm0 <- rowMeans(matrix(cvm0,ncol=max(foldid)))
      if( length(alpha)>1 )
        opt.alpha <- alpha[which.min(cvm0)]
      q.val <- min(cvm0)
    }
  }
  
  #opt.gamma <- gamma
  #opt.lambda <- lambda
  
  #=========
  if(!search.path){
    ret<-list(q.val=q.val, model=list(lambda=lambda, gamma=gamma, alpha=opt.alpha, ipf=ipf, nfolds=nfolds))
  }else{
    ret<-list(q.val=q.val, model=list(lambda=lambda, gamma=gamma, alpha=opt.alpha, ipf=ipf, nfolds=nfolds, search.cvm=c(lambda,cvm0)))
  }
  return(ret)
}



