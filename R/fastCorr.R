#' mixlasso
#' @title Sub-function \code{fastCorr()} for tree-lasso
#' @description
#' Based on the matlab code from http://www.cs.cmu.edu/~sssykim/softwares/softwares.html
#' 
#' @export
fastCorr <- function(A){
  C <- crossprod(scale(A))/(dim(A)[1]-1)
  return(C)
}