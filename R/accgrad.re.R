#' mixlasso
#' @title Sub-function \code{accgrad()} for tree-lasso
#' @importFrom RSpectra eigs
#' @description
#' A version of \code{accgrad()} with random effects for clustered samples
#' 
#' @export
accgrad.re <- function(y, x, lambda, Tree, t.idx, C, g_idx, TauNorm,  mu, NoVar, option, num_nonpen=0, intercept=TRUE, y_mis=NULL, x_mis=NULL, t.glasso=FALSE, alpha=1, gamma=0, cov.proxy="FL", L=0, predict.re=FALSE){
  
  Y <- data.matrix(y)
  X <- data.matrix(x)
  p <- dim(X)[2]
  K <- dim(Y)[2]
  
  #L <- eigen(XX)$values[1] + lambda^2 * TauNorm/mu
  C <- C * lambda #* alpha
  
  for(i in unique(t.idx)){
    if(i == 1) Z <- matrix(1,ncol=1,nrow=sum(t.idx==1))
    if(i > 1) Z <- bdiag(Z, matrix(1,ncol=1,nrow=sum(t.idx==i)))
  }
  
  if(intercept){
    if(num_nonpen>0){
      nonvs.idx <- rep(1:(num_nonpen+1), length(unique(t.idx))) + rep((1:length(unique(t.idx))-1)*(1+p),each=num_nonpen+1)
    }else{
      nonvs.idx <- rep(1, length(unique(t.idx))) + (1:length(unique(t.idx))-1)*(1+p)
    }
    
  }else{
    if(num_nonpen>0){
      nonvs.idx <- rep(1:num_nonpen, length(unique(t.idx))) + rep((1:length(unique(t.idx))-1)*p,each=num_nonpen)
    }else{
      nonvs.idx <- -c(1:(p*length(unique(t.idx))))
    }
  }
  
  # sigma.eps <- sum(Y.diff.mean^2)/(n-max(t.idx))/K
  # sigma.re <- (sum(Y.mean^2)/K - sum(1/as.vector(table(t.idx)))*sigma.eps)/max(t.idx)
  if(predict.re){
    if(cov.proxy=="FL"){
      V.full <- chol2inv(chol(diag(dim(Y)[1])+log(dim(Y)[1])*Matrix::tcrossprod(Z,Z)))
      #V.full <- diag(dim(Y)[1])+log(dim(Y)[1])*Matrix::tcrossprod(Z,Z)
    }else{
      if(cov.proxy=="BCG"){
        V.full <- chol2inv(chol(diag(dim(Y)[1])+2/(3*length(unique(t.idx)))*Matrix::tcrossprod(Z,Z)))
        #V.full <- diag(dim(Y)[1])+2/(3*length(unique(t.idx)))*Matrix::tcrossprod(Z,Z)
      }else{
        stop("Please specify correct argument 'cov.proxy'!")
      }
    }
    V.full <- data.matrix(V.full)
  }else{
    V.full <- diag(dim(Y)[1])
  }
  XX <- XY <- rep(list(NA), length(unique(t.idx)))
  
  if( length(L) == 1 ) L <- rep(0, length(unique(t.idx)))
  
  for(i in unique(t.idx)){
    #XX[[i]] <- as( Matrix::crossprod( t(Matrix::crossprod(X[t.idx==i,], V.full[t.idx==i,t.idx==i])), X[t.idx==i,]), "dgCMatrix")
    #XY[[i]] <- as( Matrix::crossprod( t(Matrix::crossprod(X[t.idx==i,], V.full[t.idx==i,t.idx==i])), Y[t.idx==i,]), "dgCMatrix")
    
    #L[i] <- as.numeric(RSpectra::eigs(XX[[i]], 1, which="LR", opts = list(retvec = FALSE))$values) + (lambda)^2 * TauNorm/mu + (lambda*(1-alpha))^2 * ifelse(t.glasso,Tglasso$TauNorm,0)/mu#/10
    #L[I] <- as.numeric(RSpectra::eigs(XX[[i]], 1, which="LR", opts = list(retvec = FALSE))$values) + (lambda)^2 * TauNorm/mu + (lambda*(1-alpha))^2 * ifelse(t.glasso,sum(p),0)/mu#/10
    #L[i] <- (lambda)^2 * TauNorm/mu + (lambda*(1-alpha))^2 * ifelse(t.glasso,sum(p),0)/mu
    if( L[i]==0 ){
      XX <- as( Matrix::crossprod( t(Matrix::crossprod(X[t.idx==i,], V.full[t.idx==i,t.idx==i])), X[t.idx==i,]), "dgCMatrix")
      L[i] <- as.numeric(RSpectra::eigs(XX, 1, which="LR", opts = list(retvec = FALSE))$values)# + (lambda)^2 * TauNorm/mu + (lambda*(1-alpha))^2 * ifelse(t.glasso,Tglasso$TauNorm,0)/mu#/10
    }
    #L[i] <- L[i] + (lambda)^2 * TauNorm/mu + (1-alpha) * gamma * sqrt(length(unique(t.idx))) * ifelse(t.glasso,sum(p),0)/mu
    L[i] <- L[i] + (lambda)^2 * TauNorm/mu + ((1-alpha) * gamma)^2 * length(unique(t.idx))/mu
    
  }
  
  #browser()
  #mixlassoObj <- list(V.full=V.full, C=C, g_idx=g_idx, t.glasso=t.glasso, t.idx=t.idx, X=X, Y=Y, intercept=intercept, num_nonpen=num_nonpen, L=L, lambda=lambda, option=option, mu=mu, NoVar=NoVar, gamma=gamma, y_mis=y_mis, alpha=alpha)
  #save(mixlassoObj, file="mixlassoObj.RData")
  bx <- mixlassoLoop(V.full, C, g_idx, t.glasso, t.idx, X, Y, intercept, num_nonpen, L, lambda, option, mu, NoVar, gamma, y_mis, alpha ) 
  
  # obj <- 0.
  # for(iter in 1:option["maxiter"]){
  #   
  #   #if(sigma.re < 0) browser()#sigma.re <- 1
  #   ## run one step gradient
  #   theta_new <- 2./(iter+2.)
  #   
  #   obj_new <- 0.
  #   onegradLoop <- onegradLoop( t.idx, X, Y, XX, XY, bx, intercept, num_nonpen, RC_tree, L, lambda, RC_Tglasso, theta, theta_new , option, alpha, l1_alpha, y_mis ) 
  #   
  #   bx_new <- onegradLoop$bx_new
  #   obj_new <- onegradLoop$obj_new
  #   #sigma.eps.new <- sigma.re.new <- 0
  #   #sigma.re=2; sigma.eps=1;
  #   # for(i in 1:max(t.idx)){
  #   #   
  #   #   #V <- matrix(sigma.re, nrow=sum(t.idx==i), ncol=sum(t.idx==i))
  #   #   #diag(V) <- sigma.re + sigma.eps
  #   #   
  #   #   #V <- matrix(log(n), nrow=sum(t.idx==i), ncol=sum(t.idx==i))
  #   #   #diag(V) <- diag(V) + 1
  #   #   #V <- chol2inv(chol(V))
  #   #   
  #   #   if(intercept){
  #   #     p_i <- (i-1)*(1+p)+1:(1+p)
  #   #   }else{
  #   #     p_i <- (i-1)*p+1:p
  #   #   }
  #   #   
  #   #   if(t.glasso){
  #   #     #onegrad.obj <- onegrad(X=X[t.idx==i,], Y=Y[t.idx==i,], XX=XX[[i]], XY=XY[[i]], bx=bx[p_i,], intercept=intercept, num_nonpen=num_nonpen,
  #   #     #                       RC_tree=RC_tree[(i-1)*p+1:p,], mu=mu, L=L[i], lambda=lambda, RC_Tglasso=RC_Tglasso[(i-1)*(p-num_nonpen)+1:(p-num_nonpen),], 
  #   #     #                       theta=theta, theta_new=theta_new, option=option, alpha=alpha, l1_alpha=l1_alpha, y_mis=y_mis[t.idx==i,])
  #   #     onegrad.obj <- onegrad(X=X[t.idx==i,], Y=Y[t.idx==i,], XX=XX[[i]], XY=XY[[i]], bx=bx[p_i,], intercept=intercept, num_nonpen=num_nonpen,
  #   #                            RC_tree=RC_tree[(i-1)*p+1:p,], L=L[i], lambda=lambda, RC_Tglasso=RC_Tglasso[(i-1)*(p-num_nonpen)+1:(p-num_nonpen),], 
  #   #                            theta=theta, theta_new=theta_new, option=option, alpha=alpha, l1_alpha=l1_alpha, y_mis=y_mis[t.idx==i,])
  #   #   }else{
  #   #     onegrad.obj <- onegrad(X=X[t.idx==i,], Y=Y[t.idx==i,], XX=XX[[i]], XY=XY[[i]], bx=bx[p_i,], intercept=intercept, num_nonpen=num_nonpen,
  #   #                            RC_tree=RC_tree[(i-1)*p+1:p,], L=L[i], lambda=lambda, RC_Tglasso=RC_Tglasso, 
  #   #                            theta=theta, theta_new=theta_new, option=option, alpha=alpha, l1_alpha=l1_alpha, y_mis=y_mis[t.idx==i,])
  #   #   }
  #   #   
  #   #   bx_new[p_i,] <- onegrad.obj$bx_new
  #   #   obj_new <- obj_new + onegrad.obj$obj_new
  #   #   
  #   #   # ## calculate the general variance of residuals and inner group variance
  #   #   # if(intercept){
  #   #   #   sigma.eps.new <- sigma.eps.new + sum((Y.diff.mean[t.idx==i,] - crossprod(t(X.diff.mean[t.idx==i,]), onegrad.obj$bx_new[-1,]))^2)
  #   #   #   sigma.re.new <- sigma.re.new + sum((Y.mean[i,] - onegrad.obj$bx_new[1,] - crossprod(matrix(X.mean[i,],ncol=1), onegrad.obj$bx_new[-1,]))^2)
  #   #   #   
  #   #   # #   epsilon[t.idx==i,] <- Y.diff.mean[t.idx==i,] - crossprod(t(X.diff.mean[t.idx==i,]), onegrad.obj$bx_new[-1,])
  #   #   # #   ub[t.idx==i,] <- matrix(Y.mean[i,] - onegrad.obj$bx_new[1,] - crossprod(matrix(X.mean[i,],ncol=1), onegrad.obj$bx_new[-1,]),nrow=sum(t.idx==i),ncol=K,byrow=T)*sqrt(sum(t.idx==i))
  #   #   # }else{
  #   #   #   sigma.eps.new <- sigma.eps.new + sum((Y.diff.mean[t.idx==i,] - crossprod(t(X.diff.mean[t.idx==i,]), onegrad.obj$bx_new))^2)
  #   #   #   sigma.re.new <- sigma.re.new + sum((Y.mean[i,] - crossprod(matrix(X.mean[i,],ncol=1), onegrad.obj$bx_new))^2)
  #   #   # }
  #   # }
  #   #sigma.eps <- sigma.eps.new/(dim(Y)[1]-max(t.idx))/K
  #   #sigma.re <- (sigma.re.new/K - sum(1/as.vector(table(t.idx)))*sigma.eps)/max(t.idx)
  #   
  #   # sigma.eps <- sigma.re <- 0
  #   # for(j in 1:K){
  #   #   sigma.eps <- sigma.eps + sum(epsilon[,j]^2)/(n-max(t.idx)-sum(bx_new[-1,j]!=0))
  #   #   sigma.re <- sigma.re + (as.vector(crossprod(t(crossprod(matrix(ub[,j],ncol=1), P)), matrix(ub[,j],ncol=1))) - (n-sum(bx_new[-1,j]!=0))*sigma.eps)/(n-max(t.idx))
  #   # }
  #   # sigma.eps <- sigma.eps/K
  #   # sigma.re <- sigma.re/K
  #   
  #   #compute grad(f(w_k)) for the tree-lasso penalty
  #   RC_tree <- as( Matrix::crossprod( shrink(Matrix(Matrix::tcrossprod(C,bx_new[-nonvs.idx,])/mu,sparse=T), g_idx), C), "dgCMatrix")
  #   
  #   obj_new <- obj_new + cal2norm(Matrix(Matrix::tcrossprod(C,bx_new[-nonvs.idx,]),sparse=T), g_idx)
  #   
  #   # update grad(f(w_k)) for the tissue-group-lasso penalty
  #   if(t.glasso){
  #     obj_new <- obj_new + cal2norm(Matrix(Matrix::crossprod(t(Tglasso$C),bx_new[-nonvs.idx,]),sparse=T), Tglasso$g_idx)
  # 
  #     R.Tglasso <- shrink(Matrix(Matrix::crossprod(t(Tglasso$C),bx_new[-nonvs.idx,])/mu,sparse=T), Tglasso$g_idx)
  #     RC_Tglasso <- Matrix::crossprod(Tglasso$C, R.Tglasso)
  #   }
  #   #cat("iter=",iter,"; L=", L, "; obj_new-obj=",obj_new-obj, "; obj_new=",obj_new,"; onegrad.obj$obj_new=",onegrad.obj$obj_new,"\n", sep="")
  #   #time0 <- proc.time() - ptm
  #   #time[iter] <- time0$user
  #   # if(iter %% 50 == 0)
  #      #cat("iter=",iter,"; L=", L, "; sigma.eps=", sigma.eps, "; sigma.re=", sigma.re, "; obj_new-obj=",obj_new-obj,"\n", sep="")
  #      cat("iter=",iter,"; L=", L, "; obj_new-obj=",obj_new-obj,"\n", sep="")
  #   # 
  #   
  #   if((iter>10) && (abs(obj_new-obj)/abs(obj)<option["tol"])) break
  #   theta <- theta_new
  #   bx <- bx_new
  #   obj <- obj_new
  #   
  #   if(iter==option["maxiter"])
  #     warning("Stopped at the maximum number of iterations!")
  # }
  
  bx[abs(bx) < option["threshold"]] <- 0
  
  # ## calculate the general variance of residuals and inner group variance
  sigma.eps <- NULL#0
  # for(i in 1:max(t.idx)){
  #   if(intercept){
  #     sigma.eps <- sigma.eps + sum((Y.diff.mean[t.idx==i,] - Matrix::crossprod(t(X.diff.mean[t.idx==i,]), bx[(i-1)*(1+p)+1:p,]))^2)
  #   }else{
  #     sigma.eps <- sigma.eps + sum((Y.diff.mean[t.idx==i,] - Matrix::crossprod(t(X.diff.mean[t.idx==i,]), bx[(i-1)*p+1:p,]))^2)
  #   }
  # }
  # sigma.eps <- sigma.eps/(dim(Y)[1]-max(t.idx))/K
  
  # u.re <- rep(0,max(t.idx)^2)
  # for(i in 1:max(t.idx)){
  #   u.first <- K*data.matrix(Matrix::crossprod(Z[t.idx == i,],Z[t.idx == i,]))
  #   diag(u.first) <- diag(u.first) + 1/log(dim(Y)[1])
  #   u.first <- chol2inv(chol(u.first))
  #   if(intercept){
  #     u.last <- rowSums(Y[t.idx == i,] - data.matrix(Matrix::crossprod(t(cbind(rep(1,sum(t.idx == i)),X[t.idx == i,])), bx[(i-1)*(1+p)+1:(1+p),])))
  #   }else{
  #     u.last <- rowSums(Y[t.idx == i,] - data.matrix(Matrix::crossprod(t(X[t.idx == i,]), bx[(i-1)*p+1:p,])))
  #   }
  #   u.re[max(t.idx)*(i-1)+1:max(t.idx)] <- Matrix::crossprod( t(Matrix::tcrossprod(u.first,Z[t.idx==i,])), matrix(u.last, ncol=1) )
  # }
  #browser()
  u.re <- rep(0,length(unique(t.idx)))
  if(predict.re){
    for(i in unique(t.idx)){
      u.first <- K * Matrix::crossprod(Z[t.idx == i,i])
      #diag(u.first) <- diag(u.first) + 1/log(dim(Y)[1])
      #u.first <- chol2inv(chol(u.first))
      u.first <- (u.first + 1/log(dim(Y)[1]))^(-1)
      if(intercept){
        u.last <- rowSums(Y[t.idx == i,] - data.matrix(Matrix::crossprod(t(cbind(rep(1,sum(t.idx == i)),X[t.idx == i,])), bx[(i-1)*(1+p)+1:(1+p),])))
      }else{
        u.last <- rowSums(Y[t.idx == i,] - data.matrix(Matrix::crossprod(t(X[t.idx == i,]), bx[(i-1)*p+1:p,])))
      }
      #u.re[max(t.idx)*(i-1)+1:max(t.idx)] <- Matrix::crossprod( t(Matrix::tcrossprod(u.first,Z[t.idx==i,])), matrix(u.last, ncol=1) )
      u.re[i] <- u.first * Z[t.idx == i,i] %*% matrix(u.last, ncol=1)
    }
  }
  
  #Beta <- bx
  #obj <- obj[1:iter]
  #time <- time[1:iter]
  #return(list(Beta=Beta, obj=obj, time=time, iter=iter))
  return(list(Beta=bx, u.re=u.re, sigma=sigma.eps))
}
