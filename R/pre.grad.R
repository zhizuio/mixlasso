#' mixlasso
#' @title Sub-function \code{pre.grad()} for tree-lasso
#' @description
#' Based on the matlab code from http://www.cs.cmu.edu/~sssykim/softwares/softwares.html
#' 
#' @importFrom Matrix sparseMatrix
#' 
#' @export
pre.grad <- function(Tree, Tw){#, t.glasso=FALSE){
  
  V <- dim(Tree)[1]
  K <- dim(Tree)[2]
  
  sum_col_T <- as.integer(apply(Tree, 1, sum))
  SV <- sum(sum_col_T)
  csum <- as.integer(cumsum(sum_col_T))
  #g_idx <- cbind(c(1,csum[1:(length(csum)-1)]+1), csum, sum_col_T)
  if(length(csum)!=1){
    g_idx <- cbind(c(as.integer(1),csum[1:(length(csum)-1)]+1), csum, sum_col_T)
  }else{
    g_idx <- cbind(as.integer(1), csum, sum_col_T)
  }
  
  J <- rep(0, SV)
  W <- rep(0, SV)
  for(v in 1:V){
    J[g_idx[v,1]:g_idx[v,2]] = which(Tree[v,] == 1)
    W[g_idx[v,1]:g_idx[v,2]] = Tw[v]
  }
  
  C <- sparseMatrix(i=1:SV, j=J, x=W, dims=c(SV, K))
  
  # if(t.glasso){
  #   TauNorm <- matrix(Tw[1:V], nrow=V,ncol=V) * Tree[,1:V]
  #   TauNorm <- max(apply(TauNorm^2, 2, sum))
  # }else{
    TauNorm <- matrix(Tw, nrow=length(Tw),ncol=K) * Tree
    TauNorm <- max(apply(TauNorm^2, 2, sum))
  #}
  
  return(list(C=C, g_idx=g_idx, TauNorm=TauNorm))
}
# function accgrad2 calculates the norm term in Lipschitz constant of the irectly adapted tree-lasso method for two data sources
pre.grad2 <- function(Tree, Tw, lambda){
  
  if(is.null(dim(Tree))) Tree <- matrix(Tree, nrow=1)
  
  V <- dim(Tree)[1]
  K <- dim(Tree)[2]
  
  sum_col_T <- apply(Tree, 1, sum)
  SV <- sum(sum_col_T)
  csum <- cumsum(sum_col_T)
  if(length(csum)!=1){
    g_idx <- cbind(c(1,csum[1:(length(csum)-1)]+1), csum, sum_col_T)
  }else{
    g_idx <- cbind(1, csum, sum_col_T)
  }
  
  J <- rep(0, SV)
  W <- rep(0, SV)
  for(v in 1:V){
    J[g_idx[v,1]:g_idx[v,2]] = which(Tree[v,] == 1)
    W[g_idx[v,1]:g_idx[v,2]] = Tw[v]
  }
  
  C <- sparseMatrix(i=1:SV, j=J, x=W, dims=c(SV, K))
  
  TauNorm0 <- matrix(Tw, ncol=1)
  for(r in 2:K) TauNorm0 <- cbind(TauNorm0, Tw)
  TauNorm <- TauNorm0 * Tree
  TauNorm <- rbind(lambda[1]*TauNorm, lambda[2]*TauNorm)
  TauNorm <- max(apply(TauNorm^2, 2, sum))
  
  return(list(C=C, g_idx=g_idx, TauNorm=TauNorm))
}