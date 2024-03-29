#' mixlasso
#' @title Wrapper function for glmnet objects.
#' @description
#' Wrapper function for glmnet objects used by epsgo function. This function is mainly used within the function \code{epsgo}. See the \code{R} package \pkg{c060} for details.
#' 
#' @importFrom Matrix Diagonal bdiag
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom parallel detectCores makeCluster	stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' 
#' @export
tune.glmnet.interval<-function(parms, x, y,
                                         weights, 
                                         x_test=NULL,
                                         y_test=NULL,
                                         offset = NULL, 
                                         lambda = NULL,
                                         num.nonpen = 0,
                                         nfolds = 10,
                                         type.measure = "mse",
                                         seed=12345, 
                                         grouped = TRUE, 
                                         type.min="lambda.min",
                                         family="gaussian",
                                         foldid=NULL,
                                         intercept=TRUE,
                                         standardize.response=FALSE,
                                         p=NULL,
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
  }else{
    alpha<-1
  }
  if(!identical(grep("ipf", names(parms)), integer(0))){
    ipf<-parms[grep("ipf", names(parms))]
    if(verbose) print(paste("IPF=",paste(as.character(ipf),collapse=":")))
  }else{
    ipf<-NA
  }
  
  #  2. find optimal lambda for given alpha and penalty factors ###########################
  # find optimal lambda for given alpha
  set.seed(seed)
  if(is.null(foldid)) foldid <- sample(rep(seq(nfolds),length=dim(x)[1]))
    
    if(length(alpha)==1){
      if(is.null(lambda))
        if(length(p) == 1){
          adpen0 <- rep(1, p)
        }else{
          adpen0 <- c(rep(1,p[1]), rep(ipf,p[2]))
        }
        lambda <- glmnet(x=x,y=y, family=family,alpha=alpha,nlambda=20,intercept=intercept,offset=offset,penalty.factor=adpen0,standardize.response=standardize.response)$lambda
        cv<-cv.glmnet(x=x,y=y, family=family,
                  alpha=alpha,
                  lambda=lambda,
                  offset = NULL,
                  type.measure =type.measure,
                  foldid = foldid,
                  grouped = grouped,
                  penalty.factor=adpen0,
                  intercept=intercept,
                  standardize.response=standardize.response, 
                  parallel=parallel)
       
      opt.lambda<-ifelse(type.min=="lambda.min", cv$lambda.min, cv$lambda.1se)
      
      if((opt.lambda<=min(cv$lambda)) | (opt.lambda>=max(cv$lambda))) warning("Searching lambda reaches the boundary value of the given lambda sequence!")
      lambda<-cv$lambda
      if(type.measure=="auc") cv$cvm <- 1-cv$cvm
      cvm0<-cv$cvm
      q.val<-min(cv$cvm)
    }
    if(length(alpha)>1){
      if(is.null(lambda)){
        #======
        # Method 1 to give lambda sequence
        #======
        lambda <- NULL
        for(s in 1:length(p)) lambda<-c(lambda, glmnet(x=x_star,y=y_star,family=family,alpha=alpha[s],nlambda=20,intercept=FALSE,penalty.factor=adpen0)$lambda)
        lambda <- sort(lambda, decreasing=T)
        lambda <- lambda[floor(seq(1,length(lambda),length=15))]
        #======
        # Method 2 to give lambda sequence as sIPFEN
        #======
        #lambda <- glmnet(x=x_star,y=y_star, family=family,alpha=mean(alpha),nlambda=20,intercept=FALSE,offset=offset,penalty.factor=adpen0,standardize.response=standardize.response)$lambda
      }
      lambda <- matrix(lambda,ncol=1)
      
      alpha[alpha==0] <- 0.0001 
      # implement k-fold cross validation and search the lambda sequence
      cvHeteAlpha<-function(xx,la){
        y2 <- as.vector(rbind(y[!foldid==xx,], matrix(0, nrow=sum(p),ncol=m)))
        sd.y2 <- sd(y2)*sqrt((length(y2)-1)/length(y2))
        y2 <- y2/sd.y2
         
        adpen00 <- adpen00*(dim(x_star)[2])/sum(adpen00)
        lambda2 <- la[1]*c(1,ipf)*m*dim(x)[1]
        # penalty parameters ratios
        pf.r <- c(1, ipf)
        x2 <- aug <- NULL
        # generate the transformed X
        for(s in 1:length(p)){
          if(s==1){
            x2 <- pf.r[s]/alpha[s]*x[,num.nonpen+(ifelse(s==1,1,cumsum(p)[s-1]+1)):cumsum(p)[s]]
            aug <- pf.r[s]/alpha[s]*sqrt(1/2*lambda2[s]*(1-alpha[s]))*diag(p[s])
          } 
          if(s>1){
            x2 <- cbind(x2, pf.r[s]/alpha[s]*x[,num.nonpen+(ifelse(s==1,1,cumsum(p)[s-1]+1)):cumsum(p)[s]])
            aug <- bdiag(aug, pf.r[s]/alpha[s]*sqrt(1/2*lambda2[s]*(1-alpha[s]))*diag(p[s]))
          } 
        }
        if(num.nonpen>0){
          x2 <- cbind(x[,c(1:num.nonpen)], x2)
          aug <- cbind(matrix(0,ncol=num.nonpen,nrow=sum(p)), aug)
        }
        if(intercept){
          org.const <- c(rep(1,dim(x2)[1])[!foldid==xx],rep(0,sum(p)))
        }else{
          org.const <- NULL
        }
        x2_star <- kronecker(Diagonal(m), cbind(org.const,rbind(x2[!foldid==xx,], aug))/sd.y2)
        lam_star <- lambda2[1]/(prod(dim(y))+sd.y2^2*sum(dim(x))*m*2)
        la2<- glmnet(x2_star,y2,family=family,lambda=seq(lam_star*0.8,lam_star*1.2,length=11),maxit=1000000,intercept=FALSE,penalty.factor=adpen00)
        
        if(sum(foldid==xx)==0){
          idx_test <- rep(TRUE, dim(x)[1])
        }else{
          idx_test <- c(foldid==xx)
        } 
        Beta0 <- matrix(la2$beta[,ceiling(11/2)], ncol=dim(y)[2])
        for(s in 1:length(p)) Beta0[n.int.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s],] <- pf.r[s]/alpha[s]*Beta0[n.int.nonpen+ifelse(s==1,1,cumsum(p)[s-1]+1):cumsum(p)[s],]
        res <- sum((y[idx_test,] - cbind(rep(ifelse(intercept,1,NULL),sum(idx_test)),x[idx_test,]) %*% Beta0)^2)/(sum(idx_test)*(dim(y)[2]))
        return(res)
      }
      laSeqHeteAlpha<-function(la) {mean(sapply(1:max(foldid), cvHeteAlpha, la))}
      
      # parallelise the lambda sequence in Abel (HPC)
      cvm0 <- numeric(length(lambda))
      if(sum(parallel)==1){
        cl <- makeCluster(min(length(lambda),15))
        registerDoParallel(cl)
        cvm0[1:min(length(lambda),15)] <- foreach(i = 1:min(length(lambda),15), .combine=c, .packages= c('base','Matrix')) %dopar%{
          adpen00 <- adpen00
          laSeqHeteAlpha(lambda[i])
        }
        stopCluster(cl)
        
        if(length(lambda)>15){
          cvm0[(15+1):length(lambda)]<-apply(matrix(lambda[(15+1):length(lambda)],ncol=1), 1, laSeqHeteAlpha)
        }
      }
      if(sum(parallel)!=1){
        cvm0<-apply(lambda, 1, laSeqHeteAlpha)
      }
      
      opt.lambda<-lambda[which.min(cvm0)]#*max(adpen00)
      # q.val= mean cross-validated error over the folds
      q.val<-min(cvm0)
      
      # cv below is useless currently
      cv<-cv.glmnet(x=x_star,y=y_star, family=family,
                    alpha=1,
                    offset = offset,
                    type.measure =type.measure,
                    foldid = foldid2,
                    grouped = grouped,
                    penalty.factor=adpen0,
                    intercept=FALSE,
                    parallel=parallel)

      # prediction with all dataset
      if(is.null(x_test)){
        x_test <- x
        y_test <- y
      }
    }

  if(!search.path){
    ret<-list(q.val=q.val, model=list(alpha=alpha, lambda=opt.lambda, ipf=ipf, p=p, nfolds=nfolds, cvreg=cv))
  }else{
    ret<-list(q.val=q.val, model=list(alpha=alpha, lambda=opt.lambda, ipf=ipf, p=p, nfolds=nfolds, cvreg=cv, search.cvm=c(lambda,cvm0)))
  }
  
  return(ret)
}



