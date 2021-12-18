#' mixlasso
#' @title Wrapper function for conditional logistic lasso objects.
#' @description
#' Wrapper function for conditional logistic lasso objects optized based on the full likelihood by the EPSGO algorithm. This function is mainly used within the function \code{epsgo}. 
#' 
#' @importFrom glmnet stratifySurv
#' @importFrom survival Surv
#' @importFrom foreach foreach %dopar%
#' @importFrom stats as.formula
#' @importFrom parallel detectCores makeCluster	stopCluster
#' @importFrom doParallel registerDoParallel
#' 
#' @param x,y input matrix.
#' @param family response type.
#' @param lambda optional user-supplied \code{lambda} sequence; default is NULL, and \code{espsgo} chooses its own sequence.
#' @param intercept should  intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE}).
#' @param strata.surv stratification variable for the Cox survival model.
#' @param foldid an vector of values for the cross-validation.
#' @param standardize.response standardization for the response variables. Default: \code{TRUE}.
#' @param p the number of predictors from different data source.
#' @param verbose print the middle search information, default is \code{TRUE}.
#' @param seed random seed
#' @param parallel If \code{TRUE}, use parallel foreach to fit each fold except parallelizing each lambda for the tree-lasso methods. If \code{c(TRUE,TRUE)}, use parallel foreach to fit each fold and each lambda.
#' @return An object with list "\code{tune.clogit.interval}"
#' \item{q.val}{the minimum MSE (or minus likelihood for the Cox model) through the cross-validation}
#' \item{model}{some model related quantities:
#' \itemize{\code{alpha  }}{the optimzed alpha}
#' \itemize{\code{lambda }}{the optimzed (first) lambda}
#' \itemize{\code{ipf    }}{the optimzed penalty factors}
#' \itemize{\code{p      }}{a vector of the numbers of features from multiple data sources}
#' \itemize{\code{nfolds }}{number of folds used for the cross-validation}
#' \itemize{\code{cvreg  }}{the cross-validation results}
#' }
#' @export
tune.surv.interval<-function(parms, x=x, y=y,
                                         x_test=NULL,
                                         y_test=NULL,
                                         family = "cox",
                                         strata.surv=NULL,
                                         intercept=TRUE,
                                         alpha = 0,
                                         lambda = NULL,
                                         num.nonpen = 0,
                                         nfolds = 3,
                                         foldid =NULL,
                                         seed=12345, 
                                         standardize.response=FALSE,
                                         p=NULL,
                                         type.measure="deviance", 
                                         type.min="lambda.min",
                                         parallel=FALSE,
                                         verbose=FALSE,
                                         search.path=FALSE,
                                         ...){
  
  # 1. decode the parameters ############################################################
  #x <- data.matrix(x)
  #y <- data.matrix(y)
  
  if(!identical(grep("alpha", names(parms)), integer(0))){
    alpha<-parms[grep("alpha", names(parms))]
    if(verbose) print(paste("alpha=",paste(as.character(alpha),collapse=",")))
  }#else{
  #  alpha<-1
  #}
  if(!identical(grep("ipf", names(parms)), integer(0))){
    ipf<-parms[grep("ipf", names(parms))]
    if(verbose) print(paste("IPF=",paste(as.character(ipf),collapse=":")))
  }else{
    ipf<-NA
  }
  #browser()
  adpen0 <- c(rep(0,num.nonpen), rep(1, sum(p)))
  if(length(p) > 1){
    for(i in 1:(length(p)-1)){
      adpen0[num.nonpen+(cumsum(p)[i]+1):cumsum(p)[i+1]] <- ipf[i]
    }
  }
  
  if(is.null(strata.surv)){
    y2 <- y
  }else{
    y2 <- glmnet::stratifySurv(y, strata.surv)
  }
  
  #  2. find optimal lambda for given alpha and penalty factors ###########################
  # find optimal lambda for given alpha
  set.seed(seed)
  if(!is.null(foldid)){
    nfolds <- max(foldid)
  }else{
    foldid <- rep(sample(rep(seq(nfolds),length(table(strata.surv))/nfolds)), each=2)
  }
  
  if(length(alpha)==1){
    
    if(is.null(lambda)){
      lambda <- glmnet(x=x,y=y2, family=family,alpha=alpha,nlambda=20,intercept=intercept,offset=offset,penalty.factor=adpen0,standardize.response=standardize.response)$lambda
    }
    cv<-cv.glmnet(x=x,y=y2, family=family,
                  alpha=alpha,
                  lambda=lambda,
                  offset = NULL,
                  type.measure =type.measure,
                  nfolds = nfolds,
                  foldid = foldid,
                  penalty.factor=adpen0,
                  intercept=intercept,
                  standardize.response=standardize.response, 
                  parallel=parallel)
    
    opt.lambda<-ifelse(type.min=="lambda.min", cv$lambda.min, cv$lambda.1se)
    
    #if((opt.lambda<=min(cv$lambda)) | (opt.lambda>=max(cv$lambda))) warning("Searching lambda reaches the boundary value of the given lambda sequence!")
    lambda<-cv$lambda
    if(type.measure=="auc") cv$cvm <- 1-cv$cvm
    cvm0<-cv$cvm
    q.val<-min(cv$cvm)
  }
  if(length(alpha)>1){
    stop("To be developed...")
  }
  
  if(!search.path){
    ret<-list(q.val=q.val, model=list(alpha=alpha, lambda=opt.lambda, ipf=ipf, p=p, nfolds=nfolds, cvreg=cv))
  }else{
    ret<-list(q.val=q.val, model=list(alpha=alpha, lambda=opt.lambda, ipf=ipf, p=p, nfolds=nfolds, cvreg=cv, search.cvm=c(lambda,cvm0)))
  }
  
  return(ret)
}
