#' mixlasso
#' @title Sub-function \code{tree.parms()} for tree-lasso
#' @description
#' Tree structure from hierarchical cluster analysis for tree-lasso
#' 
#' @importFrom stats dist hclust
#' 
#' @export
tree.parms <- function(y=y, h=0.7, impute=FALSE){
  if(sum(is.na(y))>0){
    if(impute){
      y <- fill.USVT(y)$X
    }else{
      y0 <- y
      m0 <- ncol(y0)
      eps <- 0.05
      # r1.na is better to be not smaller than r2.na
      r1.na <- 0.3
      r2.na <- 0.2
      k <- 1
      while(sum(is.na(y0))>0){
        r1.na <- r1.na - eps/k
        r2.na <- r1.na - eps/k
        k <- k + 1
        ## select drugs with <30% (decreasing with k) missing data overall cell lines
        na.y <- apply(y0, 2, function(xx) sum(is.na(xx))/length(xx))
        while(sum(na.y<r1.na)<m0){
          y0 <- y0[,-c(which(na.y>=r1.na))]
          m0 <- sum(na.y<r1.na)
          na.y <- apply(y0, 2, function(xx) sum(is.na(xx))/length(xx))
        }
        
        ## select cell lines with treatment of at least 80% (increasing with k) drugs
        na.y0 <- apply(y0, 1, function(xx) sum(is.na(xx))/length(xx))
        while(sum(na.y0<r2.na)<(dim(y0)[1])){
          y0 <- y0[na.y0<r2.na,]
          na.y0 <- apply(y0, 1, function(xx) sum(is.na(xx))/length(xx))
        }
        num.na <- sum(is.na(y0))
        #message("#{NA}=", num.na, "\n", "r1.na =", r1.na, ", r2.na =", r2.na, "\n")
      }
      y <- y0
    }
    
  }
  # pdf("ctrp.dendrogram.pdf", width = 18, height=4)
  # tree <- hclust(dist(t(y)), "ave")
  # k <- 4
  # cl_members <- cutree(tree = tree, k = k)
  # plot(x = tree, labels =  row.names(tree), cex = 0.5)
  # rect.hclust(tree = tree, k = k, which = 1:k, border = 2:(k+1), cluster = cl_members)
  # dev.off()
    
  m <- dim(y)[2]
  myDist0 <- 1 - abs(fastCorr(y))
  myDist <- myDist0[lower.tri(myDist0)]
  a0 <- dist(t(y))
  a0[1:length(a0)] <- myDist
  # hierarchical clustering for multivariate responses
  myCluster0 <- hclust(a0, method = "complete")
  myCluster <- cbind(ifelse(myCluster0$merge < 0, - myCluster0$merge, myCluster0$merge + m), myCluster0$height)
  
  conv0 <- convH2T(myCluster, h)
  Tree <- conv0$Tree
  if(is.null(dim(Tree)))
    Tree <- matrix(Tree, nrow=1)
  Tw <- conv0$Tw
  idx <- c(apply(Tree,1,sum) == 1)
  Tree <- Tree[!idx,]
  if(is.null(dim(Tree)))
    Tree <- matrix(Tree, nrow=1)
  Tw <- Tw[!idx]
  
  return(list(Tree=Tree, Tw=Tw, hc=myCluster0, y.colnames=colnames(y)))
}
