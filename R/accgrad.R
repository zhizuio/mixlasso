#' mixlasso
#' @title Sub-function \code{accgrad()} for tree-lasso
#' @description
#' Based on the matlab code from http://www.cs.cmu.edu/~sssykim/softwares/softwares.html
#' 
#' @export
accgrad <- function(y, x, lambda, Tree, C, g_idx, TauNorm,  mu, option, num_nonpen=0, intercept=TRUE, y_mis=NULL, x_mis=NULL){
  
  #Y Centered Matrix: N by K
  #X Centered Matrix: N by p
  #lam: lambda
  #Tree: sparse matrix: group info. rows: number of group, cols: number of tasks
  #Tw: n_group by 1: weight for each group
  #C: note \sum_|g| by K
  #g_idx: n_group by 2, group index
  #L1, Lipschitz cond
  #TauNorm: \|\Tau\|_1,2^2 
  #mu: mu in nesterov paper
  #maxiter
  
  Y <- data.matrix(y)
  X <- data.matrix(x)
  
  C <- C * lambda 
  #p <- dim(X)[2]
  #K <- dim(Tree)[2]
  
  #treeLassoObj <- list(C=C, g_idx=g_idx, X=X, Y=Y, intercept=intercept, TauNorm=TauNorm, num_nonpen=num_nonpen, lambda=lambda, option=option, mu=mu, y_mis=y_mis)
  #save(treeLassoObj, file="treeLassoObj.RData")
  #sourceCpp("/Users/zhiz/Downloads/IPFStructPenalty/mixlasso/src/treeLassoLoop.cpp")
  #load("treeLassoObj0.RData")
  #attach(treeLassoObj)
  #browser()
  bx <- treeLassoLoop(X, Y, C, g_idx, TauNorm, intercept, num_nonpen, lambda, option, mu) 
  
  # XX <- crossprod(X)
  # XY <- crossprod(X, Y)
  # 
  # L <- eigen(XX)$values[1] + lambda^2 * TauNorm/mu
  # C <- C * lambda
  # if(intercept){
  #   bw0 <- matrix(0, ncol=dim(Y)[2])
  # }else{
  #   b0 <- bw0 <- matrix(0, nrow=0, ncol=dim(Y)[2])
  # }
  # 
  # #We can also use "bw0 <- matrix(0, ncol=K)" since "dim(Y)[2]=K"
  # if(num.nonpen==0){
  #   bw1 <- matrix(0, nrow=p, ncol=dim(Y)[2])
  #   bx <- rbind(bw0, bw1)
  # }else{
  #   bW0 <- matrix(0, nrow=num.nonpen, ncol=dim(Y)[2])
  #   bw1 <- matrix(0, nrow=p-num.nonpen, ncol=dim(Y)[2])
  #   bx <- rbind(bw0, bW0, bw1)
  # }
  # 
  # theta <- 1
  # obj <- 0
  # #ptm <- proc.time()
  # for(iter in 1:option["maxiter"]){
  #   #compute grad(f(w_k))
  #   R <- shrink(data.matrix(Matrix::tcrossprod(C,bw1))/mu, g_idx, Atranspose=FALSE, 0)
  #   
  #   if(num.nonpen==0){
  #     if(intercept) grad_bw0 <- matrix(1, ncol=dim(X)[1]) %*% (matrix(1, nrow=dim(X)[1])%*%bw0 + crossprod(t(X),bw1) - Y)
  #     grad_bw <- crossprod(t(XX), bw1) - XY + crossprod(R, C)
  #   }else{
  #     if(!is.null(bw0)) grad_bw0 <- matrix(1, ncol=dim(X)[1]) %*% (matrix(1, nrow=dim(X)[1])%*%bw0 + crossprod(t(X),rbind(bW0,bw1)) - Y)
  #     grad_bW0 <- crossprod(X[,1:num.nonpen], matrix(1, nrow=dim(X)[1])%*%bw0 + crossprod(t(X),rbind(bW0,bw1)) - Y)
  #     grad_bw <- crossprod(X[,-c(1:num.nonpen)], matrix(1, nrow=dim(X)[1])%*%bw0 + crossprod(t(X),rbind(bW0,bw1)) - Y) + crossprod(R, C)
  #   }
  #   
  #   if(intercept) b0 <- bw0 - 1/L * grad_bw0
  #   if(num.nonpen>0) B0 <- bW0 - 1/L * grad_bW0
  #   bv <- bw1 - 1/L * grad_bw
  #   #browser()
  #   b <- abs(bv)-lambda/L
  #   b[b<0] <- 0
  #   b <- data.matrix(sign(bv) * b)
  #   bx_new <- data.matrix(rbind(b0, b))
  #   if(num.nonpen>0) bx_new <- data.matrix(rbind(b0, B0, b))
  #   
  #   if(intercept){
  #     obj_new <- sum((Y-crossprod(t(cbind(rep(1,dim(X)[1]),X)),bx_new))[y.mis!=1]^2)/2/nrow(Y) + cal2norm(Matrix::tcrossprod(C,b), g_idx, Atranspose=FALSE, TgCB_T=0)
  #   }else{
  #     obj_new <- sum((Y-crossprod(t(X),bx_new))[y.mis!=1]^2)/2/nrow(Y) + cal2norm(Matrix::tcrossprod(C,b), g_idx, Atranspose=FALSE, TgCB_T=0)
  #   }
  #   
  #   theta_new <- 2/(iter+2)
  #   
  #   bw <- bx_new + (1-theta)/theta * theta_new * (bx_new-bx)
  #   
  #   #bw[abs(bw) < option["threshold"]] <- 0
  #   
  #   if(intercept) bw0 <- bw[1,]
  #   if(num.nonpen>0){
  #     bW0 <- bw[dim(b0)[1]+1:num.nonpen,]
  #     bw1 <- bw[-(1:(dim(b0)[1]+num.nonpen)),]
  #   }else{
  #     if(intercept){
  #       bw1 <- bw[-1,]
  #     }else{
  #       bw1 <- bw
  #     }
  #   }
  #   
  #   #time0 <- proc.time() - ptm
  #   #time[iter] <- time0$user
  #   if( iter %% 10 == 0)
  #     message(paste("iter=",iter,"; L=", L, "; (obj_new-obj)/obbj=",abs(obj_new-obj)/abs(obj),"\n", sep=""))
  #   
  #   if(obj!=0)
  #     if((iter>10) && (abs(obj_new-obj)/abs(obj)<option["tol"])) break
  #   theta <- theta_new
  #   bx <- bx_new
  #   obj <- obj_new
  #}
  
  bx[abs(bx) < option["threshold"]] <- 0
  Beta <- bx
  #obj <- obj[1:iter]
  #time <- time[1:iter]
  #return(list(Beta=Beta, obj=obj, time=time, iter=iter))
  return(list(Beta=Beta, sigmaEU=NULL))
}
# function accgrad2 is for updating coefficients of different data sources separately with different penalty factors
accgrad2 <- function(y, x, lambda, Tree, C, g_idx, TauNorm, p, mu, option, num.nonpen=0, intercept=TRUE){
  Tree <- data.matrix(Tree)
  Y <- data.matrix(y)
  X <- data.matrix(x)
  X <- cbind(rep(1, dim(X)[1]), X)
  X1 <- data.matrix(x[,num.nonpen+1:p[1]])
  X2 <- data.matrix(x[,-c(1:(num.nonpen+p[1]))])
  C <- data.matrix(C)
  
  # to be optimized for more data sources 
  XX1 <- crossprod(X1)
  XX2 <- crossprod(X2)
  X12 <- crossprod(X1, X2)
  X21 <- crossprod(X2, X1)
  XY1 <- crossprod(X1, Y)
  XY2 <- crossprod(X2, Y)
  
  L <- eigen(crossprod(X))$values[1] + TauNorm/mu
  
  J1 <- dim(X1)[2]
  J2 <- dim(X2)[2]
  K <- dim(Y)[2]
  
  C1 <- C * lambda[1]
  C2 <- C * lambda[2]
  if(intercept){
    bw0 <- matrix(0, ncol=dim(Y)[2])
  }else{
    b0 <- bw0 <- matrix(0, nrow=0, ncol=dim(Y)[2])
  }
  if(num.nonpen==0){
    bw1 <- bx_new1 <- matrix(0, nrow=J1, ncol=K)
    bw2 <- bx_new2 <- matrix(0, nrow=J2, ncol=K)
    bx <- rbind(bw0, bw1, bw2)
  }else{
    bW0 <- matrix(0, nrow=num.nonpen, ncol=K)
    bw1 <- bx_new1 <- matrix(0, nrow=J1, ncol=K)
    bw2 <- bx_new2 <- matrix(0, nrow=J2, ncol=K)
    bx <- rbind(bw0, bW0, bw1, bw2)
  }
  if(intercept) bw0 <- matrix(0, ncol=dim(Y)[2])
  theta <- 1
  obj <- 0
  #ptm <- proc.time()
  for(iter in 1:option["maxiter"]){
    #compute grad(f(w_k))
    R1 <- shrink(data.matrix(Matrix::tcrossprod(C1,bw1))/mu, g_idx, Atranspose=FALSE, 0)
    R2 <- shrink(data.matrix(Matrix::tcrossprod(C2,bw2))/mu, g_idx, Atranspose=FALSE, 0)
    
    if(intercept) grad_bw0 <- matrix(1, ncol=dim(X)[1]) %*% (matrix(1, nrow=dim(X)[1])%*%bw0 + X1%*%bw1 + X2%*%bw2 - Y)
    grad_bw1 <- XX1 %*% bw1 + X12 %*% bw2 - XY1 + crossprod(R1,C1)
    grad_bw2 <- XX2 %*% bw2 + X21 %*% bw1 - XY2 + crossprod(R2,C2)     
    
    if(intercept) b0 <- bw0 - 1/L * grad_bw0
    bv1 <- bw1 - 1/L * grad_bw1
    bv2 <- bw2 - 1/L * grad_bw2
    
    b1 <- abs(bv1)-lambda[1]/L
    b1[b1<0] <- 0
    bx_new1 <- data.matrix(sign(bv1) * b1)
    
    b2 <- abs(bv2)-lambda[2]/L
    b2[b2<0] <- 0
    bx_new2 <- data.matrix(sign(bv2) * b2)
    
    bx_new <- rbind(b0, bx_new1, bx_new2)
    obj_new <- sum((Y-X%*%bx_new)^2)/2 + cal2norm(Matrix::tcrossprod(C1,bx_new1), g_idx, Atranspose=FALSE, TgCB_T=0) 
               + cal2norm(Matrix::tcrossprod(C2,bx_new2), g_idx, Atranspose=FALSE, TgCB_T=0)
    
    theta_new <- 2/(iter+2)
    
    bw <- bx_new + (1-theta)/theta * theta_new * (bx_new-bx)
    if(intercept) bw0 <- bw[1,]
    
    bw1 <- bw[dim(b0)[1]+1:J1,]
    bw2 <- bw[-c(1:(J1+dim(b0)[1])),]
    
    if((iter>10) && (abs(obj_new-obj)/abs(obj)<option["tol"])) break
    theta <- theta_new
    bx <- bx_new
    obj <- obj_new
  }
  
  bx[abs(bx) < option["threshold"]] <- 0
  Beta <- bx
  #obj <- obj[1:iter]
  #time <- time[1:iter]
  
  return(list(Beta=Beta))
}
