#' mixlasso
#' @title Sub-function \code{convNd2T()} for tree-lasso
#' @description
#' Based on the matlab code from http://www.cs.cmu.edu/~sssykim/softwares/softwares.html
#' 
#' @importFrom Matrix sparseMatrix
#' 
#' @export
convNd2T <- function(Nd, w, w_max){
  # Nd : node list
  # w : a vector of weights for internal nodes
  # Tree : VxK matrix
  #	V is the number of leaf nodes and internal nodes
  #	K is the number of tasks
  #	Element (v,k) is set to 1 if task k has a membership to
  #	the cluster represented by node v. Otherwise, it's 0.
  # Tw : V vector
  
  #===========================
  find_leaves <- function(Nd, ch, K, Jt, w, Tw){
    
    for(ii in 1:length(ch)){
      if(Nd[ch[ii], 2] > K){
        leaves0 <- find_leaves(Nd, which(Nd[,1] == Nd[ch[ii], 2]), K, Jt, w, Tw) 
        Jt <- leaves0$Jt
        Tw <- leaves0$Tw
      }else
        Jt <- c(Jt, Nd[ch[ii], 2])
    }
    
    Tw[Nd[ch, 2]] <- Tw[Nd[ch, 2]] * w
    
    return(list(Jt=Jt, Tw=Tw))
  }
  #===========================
  
  # of leaf nodes
  K <- Nd[1,1] - 1
  #V = Nd(size(Nd,1),1);
  #V = Nd(size(Nd,1),1)-1;		# without the root
  if(sum(w < w_max)<1){
    V <- 1 + K
  }else{
    ind0 <- which(w < w_max)    # only the internal nodes with w<w_max
    V <- ind0[length(ind0)] + K 
  }
  
  # for leaf nodes
  I <- 1:K
  J <- 1:K
  
  Tw <- rep(1, V)
  
  # for internal nodes
  for(i in (K+1):V){
    Jt <- NULL
    
    Tw[i] <- Tw[i] * (1 - w[i-K])
    leaves0 <- find_leaves(Nd, which(Nd[,1] == i), K, Jt, w[i-K], Tw)
    Jt <- leaves0$Jt
    Tw <- leaves0$Tw
    
    I <- c(I, rep(1,length(Jt)) * i)
    J <- c(J, Jt)
  }
  
  Tree <- sparseMatrix(i=I, j=J, x=rep(1, length(I)), dims=c(V, K))
  
  return(list(Tree=Tree, Tw=Tw))
  
}