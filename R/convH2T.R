#' mixlasso
#' @title Sub-function \code{covH2T()} for tree-lasso
#' @description
#' Based on the matlab code from http://www.cs.cmu.edu/~sssykim/softwares/softwares.html
#' 
#' @export
convH2T <- function(H, w_max){
  K <- dim(H)[1] + 1
  Nd <- cbind(rep((K+1):(2*K-1), each = 2), as.vector(t(H[,1:2])))
  W_norm <- H[,3]/max(H[,3])  
  conv0 <- convNd2T(Nd, W_norm, w_max)  
  
  return(conv0)  
}